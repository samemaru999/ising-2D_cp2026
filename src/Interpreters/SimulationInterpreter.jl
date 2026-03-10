"""
    Interpreters/SimulationInterpreter.jl - シミュレーション実行 [Impure]

eDSLコマンドを実際の数値シミュレーションに変換するインタープリター

FDDのPure/Impure分離により、このモジュールのみがODEソルバー等の
副作用を持つ処理を実行する
"""

using DifferentialEquations
using LinearAlgebra
using Printf
using Statistics  # median for threshold calculation

# =============================================================================
# 安定周期解（SPS）計算
# =============================================================================

"""
    interpret(cmd::ComputeSPS) -> ComputationResult{StablePeriodicSolution}

安定周期解計算コマンドを実行

1. 初期緩和シミュレーション
2. ピーク検出による周期推定
3. 1周期分のサンプリング
"""
function interpret(cmd::ComputeSPS)
    try
        params = cmd.params
        Ntheta = cmd.Ntheta
        t_relax = cmd.t_relax
        N = params.N

        # 初期条件（乱数またはデフォルト）
        X0 = create_initial_state(N)

        # ODE問題を構築
        function dynamics!(dX, X, p, t)
            dX .= network_full_dynamics(X, params)
        end

        # 緩和シミュレーション
        prob_relax = ODEProblem(dynamics!, X0, (0.0, t_relax))
        sol_relax = solve(prob_relax, Tsit5(); abstol=1e-8, reltol=1e-8, saveat=0.1)

        if sol_relax.retcode != :Success
            return Failure{StablePeriodicSolution}("Relaxation simulation failed: $(sol_relax.retcode)")
        end

        # 緩和後の状態を取得
        X_relaxed = sol_relax.u[end]

        # 周期推定のための追加シミュレーション
        T_estimate = estimate_period(X_relaxed, params)

        # 1周期分を高精度でサンプリング（緩和前の初期値を記録）
        sps = sample_periodic_orbit(X_relaxed, params, T_estimate, Ntheta; X0_initial=X0)

        return Success(sps)

    catch e
        return Failure{StablePeriodicSolution}("SPS computation failed: $e",
            Dict("exception" => e, "stacktrace" => catch_backtrace()))
    end
end

"""
    create_initial_state(N::Int) -> Vector{Float64}

初期状態を生成
"""
function create_initial_state(N::Int)
    X0 = zeros(Ds * N)
    for i in 1:N
        # FHN振動子の典型的な初期値付近
        X0[2i-1] = 1.0 #0.0 + 0.1 * randn()   # u
        X0[2i] = 1.0#-0.5 + 0.1 * randn()    # v
    end
    return X0
end

"""
    estimate_period(X0, params) -> Float64

ゼロクロス検出により周期を推定（Cコード準拠）

v[0]（第1振動子の膜電位）が0を下から上にクロスするタイミングを検出
"""
function estimate_period(X0::Vector{Float64}, params::NetworkParams)
    N = params.N

    # 十分な時間シミュレーション（RMAX相当: 200秒）
    function dynamics!(dX, X, p, t)
        dX .= network_full_dynamics(X, params)
    end

    t_search = 200.0  # Cコード: RMAX * DT = 200000 * 0.001 = 200秒
    dt_search = 0.001  # Cコード: DT = 0.001
    prob = ODEProblem(dynamics!, X0, (0.0, t_search))
    sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8, saveat=dt_search)

    # 第1振動子のv成分（膜電位）
    v1_series = [u[2] for u in sol.u]
    t_series = sol.t

    # ゼロクロス検出（Cコード準拠: v[0]が0を下から上にクロス）
    origin_threshold = 0.0  # Cコード: OV = 0.0
    crossings = Int[]
    for i in 2:length(v1_series)
        if v1_series[i-1] < origin_threshold && v1_series[i] > origin_threshold
            push!(crossings, i)
        end
    end

    if length(crossings) < 2
        @warn "Could not detect enough zero-crossings, falling back to peak detection"
        # フォールバック: ピーク検出
        peaks = Int[]
        for i in 2:length(v1_series)-1
            if v1_series[i] > v1_series[i-1] && v1_series[i] > v1_series[i+1]
                push!(peaks, i)
            end
        end
        if length(peaks) >= 2
            periods = [t_series[peaks[i]] - t_series[peaks[i-1]] for i in 2:length(peaks)]
            return mean(periods)
        end
        @warn "Could not detect period, using default"
        return 50.0
    end

    # 連続するゼロクロス間の時間差の平均
    periods = Float64[]
    for i in 2:length(crossings)
        push!(periods, t_series[crossings[i]] - t_series[crossings[i-1]])
    end

    return mean(periods)
end

"""
    sample_periodic_orbit(X0, params, T, Ntheta) -> StablePeriodicSolution

周期軌道を等間隔位相でサンプリング
"""
function sample_periodic_orbit(X0::Vector{Float64}, params::NetworkParams,
    T::Float64, Ntheta::Int;
    X0_initial::Vector{Float64}=zeros(Ds * params.N))
    N = params.N

    function dynamics!(dX, X, p, t)
        dX .= network_full_dynamics(X, params)
    end

    # 1周期分をNtheta点でサンプリング
    saveat = range(0, T, length=Ntheta)
    prob = ODEProblem(dynamics!, X0, (0.0, T))
    sol = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-10, saveat=saveat)

    # 結果をマトリックスに格納
    Xs = zeros(Ds * N, Ntheta)
    for (k, u) in enumerate(sol.u)
        Xs[:, k] = u
    end

    # 周期性の確認と補正
    endpoint_error = norm(Xs[:, 1] - Xs[:, end])
    if endpoint_error > 1e-3
        @warn "Periodic orbit may not be closed (endpoint error = $endpoint_error)"
    end

    # 厳密な周期性のため終点を始点と一致させる
    Xs[:, end] = Xs[:, 1]

    omega = 2π / T

    return StablePeriodicSolution(N, Ntheta, T, omega, Xs, params, X0_initial)
end

"""
    simulate_one_period(X0, params, T; dt=0.001) -> (t_series, X_series)

周期Tの間、細かい時間刻みでシミュレーションして軌道を返す。
Nθによるサンプリングではなく、実際の時間積分を行う。

# Returns
- `t_series`: 時間配列
- `X_series`: 状態ベクトルの配列 (各時刻での状態)
"""
function simulate_one_period(X0::Vector{Float64}, params::NetworkParams, T::Float64;
    dt::Float64=0.001)
    N = params.N

    function dynamics!(dX, X, p, t)
        dX .= network_full_dynamics(X, params)
    end

    # 周期T分を細かい時間刻みでシミュレーション
    n_steps = Int(ceil(T / dt))
    saveat = range(0, T, length=n_steps+1)

    prob = ODEProblem(dynamics!, X0, (0.0, T))
    sol = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-10, saveat=saveat)

    t_series = sol.t
    X_series = sol.u

    return t_series, X_series
end

"""
    sample_periodic_orbit_from_simulation(X0, params, T, Ntheta; dt=0.001) -> StablePeriodicSolution

周期Tの時間分シミュレーションを行い、その結果をNtheta点に補間してSPSを作成。
simulate_one_periodで得た軌道データから等間隔位相でサンプリング。
"""
function sample_periodic_orbit_from_simulation(X0::Vector{Float64}, params::NetworkParams,
    T::Float64, Ntheta::Int;
    dt::Float64=0.001,
    X0_initial::Vector{Float64}=zeros(Ds * params.N))
    N = params.N

    # 周期分シミュレーション
    t_series, X_series = simulate_one_period(X0, params, T; dt=dt)

    # Ntheta点に補間
    theta_times = range(0, T, length=Ntheta)
    Xs = zeros(Ds * N, Ntheta)

    for (k, t_target) in enumerate(theta_times)
        # 最も近い時刻のインデックスを見つける
        idx = argmin(abs.(t_series .- t_target))
        Xs[:, k] = X_series[idx]
    end

    # 周期性の確認と補正
    endpoint_error = norm(Xs[:, 1] - Xs[:, end])
    if endpoint_error > 1e-3
        @warn "Periodic orbit may not be closed (endpoint error = $endpoint_error)"
    end

    # 厳密な周期性のため終点を始点と一致させる
    Xs[:, end] = Xs[:, 1]

    omega = 2π / T

    return StablePeriodicSolution(N, Ntheta, T, omega, Xs, params, X0_initial)
end

# =============================================================================
# 位相感受関数（PSF）計算
# =============================================================================

"""
    interpret(cmd::ComputePSF) -> ComputationResult{PhaseSensitivityFunction}

位相感受関数計算コマンドを実行

FDD: 結合タイプ（CouplingType）に基づいて多重ディスパッチで計算方法を選択
- DiffusiveCoupling: PSF.mスタイルの簡略化計算
- GeneralCoupling: 完全なヤコビアン計算

随伴方程式を解いて Q(θ) を計算
"""
function interpret(cmd::ComputePSF)
    try
        sps = cmd.sps
        n_iterations = cmd.n_iterations
        coupling_type = cmd.coupling_type
        params = sps.params

        # 多重ディスパッチでPSF計算（Pure関数を呼び出し）
        # Logic/AdjointEquation.jl の compute_psf を使用
        psf = compute_psf(sps, params, coupling_type, n_iterations)

        return Success(psf)

    catch e
        return Failure{PhaseSensitivityFunction}("PSF computation failed: $e",
            Dict("exception" => e, "stacktrace" => catch_backtrace()))
    end
end

# =============================================================================
# 従来のPSF計算（後方互換性のため残す）
# =============================================================================

"""
    iterate_adjoint_equation(Q, sps, params) -> Matrix{Float64}

随伴方程式 (A1) を1周期分解く（従来実装、GeneralCoupling用）

随伴方程式:
    ω dQ_i/dθ = -J_i(θ)^T Q_i(θ) - Σ_j M_ij(θ)^T Q_i(θ) - Σ_j N_ij(θ) Q_j(θ)

ここで:
    J_i = ∂f_i/∂x_i        (固有ダイナミクスのヤコビアン)
    M_ij = ∂g_ij/∂x_i      (結合の自己項ヤコビアン)
    N_ij = ∂g_ij/∂x_j      (結合の交差項ヤコビアン)

集団振動の場合、これらは完全ヤコビアン J_full に統合される:
    ω dQ/dθ = -J_full(θ)^T Q(θ)

位相θでの後退積分:
    dQ/dθ = -(1/ω) J_full^T Q
    Q(θ - Δθ) = Q(θ) + (Δθ/ω) J_full^T Q(θ)

Δθ = 2π/(Ntheta-1), ω = 2π/T より Δθ/ω = T/(Ntheta-1) = dt
"""
function iterate_adjoint_equation(Q::Matrix{Float64}, sps::StablePeriodicSolution,
    params::NetworkParams)
    N = sps.N
    Ntheta = sps.Ntheta
    T = sps.T
    omega = sps.omega
    Xs = sps.Xs

    Q_new = zeros(size(Q))

    # 位相ステップ Δθ = 2π/(Ntheta-1)
    # 時間ステップ dt = T/(Ntheta-1)
    # 関係: Δθ = ω * dt, したがって Δθ/ω = dt
    dtheta = 2π / (Ntheta - 1)
    dt = dtheta / omega  # = T / (Ntheta - 1)

    # 逆時間方向に積分（θ: 2π → 0）
    # 初期条件を最終点に設定
    Q_new[:, end] = Q[:, end]

    for k in Ntheta:-1:2
        X_k = Xs[:, k]

        # 完全ヤコビアン J_full = J_i + M_ii + (off-diagonal coupling)
        # network_full_jacobian は以下を含む:
        #   - 対角ブロック: FHN単体ヤコビアン + 結合の自己項
        #   - 非対角ブロック: 結合の交差項 N_ij
        J_full_k = network_full_jacobian(X_k, params)

        # 随伴方程式の後退Euler積分:
        # Q(θ - Δθ) = Q(θ) + dt * J_full^T * Q(θ)
        # (dt = Δθ/ω を使用)
        Q_new[:, k-1] = Q_new[:, k] + dt * (J_full_k' * Q_new[:, k])
    end

    # 周期境界条件: Q(0) = Q(2π)
    # （収束後は自動的に満たされるはずだが、数値的に強制）
    Q_new[:, end] = Q_new[:, 1]

    return Q_new
end

"""
    normalize_psf(Q, sps, params) -> Matrix{Float64}

PSFの正規化条件を適用: Q(θ) · F(X(θ)) = ω
"""
function normalize_psf(Q::Matrix{Float64}, sps::StablePeriodicSolution,
    params::NetworkParams)
    Ntheta = sps.Ntheta
    omega = sps.omega
    Xs = sps.Xs

    Q_normalized = copy(Q)

    for k in 1:Ntheta
        X_k = Xs[:, k]
        F_k = network_full_dynamics(X_k, params)

        current_product = dot(Q[:, k], F_k)

        if abs(current_product) > 1e-10
            scale = omega / current_product
            Q_normalized[:, k] = Q[:, k] * scale
        end
    end

    return Q_normalized
end

# =============================================================================
# 3ネットワーク結合系シミュレーション
# =============================================================================

"""
    interpret(cmd::SimulateFullDynamics) -> ComputationResult{PhaseTimeSeries}

元力学系の3ネットワーク結合シミュレーションを実行
"""
function interpret(cmd::SimulateFullDynamics)
    try
        three_params = cmd.three_network_params
        initial_states = cmd.initial_states
        t_start = cmd.t_start
        t_end = cmd.t_end

        params = three_params.network
        C = three_params.C
        epsilon = three_params.epsilon
        N = params.N

        # 初期状態を1つのベクトルに結合（3N*Ds次元）
        X0 = vcat(initial_states...)

        # ODE問題を構築
        function dynamics!(dX, X, p, t)
            # 状態を3つのネットワークに分割
            X_A = X[1:Ds*N]
            X_B = X[Ds*N+1:2*Ds*N]
            X_C = X[2*Ds*N+1:3*Ds*N]

            state = ThreeNetworkState(X_A, X_B, X_C)
            dX_A, dX_B, dX_C = compute_three_network_dynamics(state, three_params)

            dX[1:Ds*N] = dX_A
            dX[Ds*N+1:2*Ds*N] = dX_B
            dX[2*Ds*N+1:3*Ds*N] = dX_C
        end

        # シミュレーション実行
        prob = ODEProblem(dynamics!, X0, (t_start, t_end))
        sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8, saveat=1.0)

        if sol.retcode != :Success
            return Failure{PhaseTimeSeries}("Simulation failed: $(sol.retcode)")
        end

        # 結果を位相時系列に変換（位相検出は別途必要）
        time = collect(sol.t)
        n_steps = length(time)

        # プレースホルダー（実際の位相検出はSPSが必要）
        phases = zeros(N_NET, n_steps)
        phase_diffs = zeros(2, n_steps)

        pts = PhaseTimeSeries(time, phases, phase_diffs)
        return Success(pts)

    catch e
        return Failure{PhaseTimeSeries}("Full dynamics simulation failed: $e",
            Dict("exception" => e))
    end
end

# =============================================================================
# 位相差シミュレーション（縮約方程式）
# =============================================================================

"""
    interpret(cmd::SimulatePhaseDynamics) -> ComputationResult{PhaseTimeSeries}

縮約位相方程式のシミュレーションを実行
"""
function interpret(cmd::SimulatePhaseDynamics)
    try
        total_pcf = cmd.total_pcf
        phi1_0 = cmd.initial_phi1
        phi2_0 = cmd.initial_phi2
        epsilon = cmd.epsilon
        t_start = cmd.t_start
        t_end = cmd.t_end

        # ODE問題を構築
        function phase_dynamics!(dphi, phi, p, t)
            dphi1, dphi2 = reduced_phase_dynamics(phi[1], phi[2], total_pcf, epsilon)
            dphi[1] = dphi1
            dphi[2] = dphi2
        end

        phi0 = [phi1_0, phi2_0]
        prob = ODEProblem(phase_dynamics!, phi0, (t_start, t_end))
        sol = solve(prob, Tsit5(); abstol=1e-8, reltol=1e-8, saveat=1.0)

        if sol.retcode != :Success
            return Failure{PhaseTimeSeries}("Phase dynamics simulation failed: $(sol.retcode)")
        end

        time = collect(sol.t)
        n_steps = length(time)

        # 位相時系列を構築
        # φ₃ = -(φ₁ + φ₂) の関係から3つのネットワークの位相を復元
        phases = zeros(N_NET, n_steps)
        phase_diffs = zeros(2, n_steps)

        for (t_idx, u) in enumerate(sol.u)
            phi1, phi2 = u
            phi3 = -(phi1 + phi2)

            # 基準位相を0として相対位相を記録
            phases[1, t_idx] = 0.0
            phases[2, t_idx] = phi1
            phases[3, t_idx] = phi1 + phi2

            phase_diffs[1, t_idx] = normalize_phase(phi1)
            phase_diffs[2, t_idx] = normalize_phase(phi2)
        end

        pts = PhaseTimeSeries(time, phases, phase_diffs)
        return Success(pts)

    catch e
        return Failure{PhaseTimeSeries}("Phase dynamics simulation failed: $e",
            Dict("exception" => e))
    end
end

# =============================================================================
# ユーティリティ
# =============================================================================

"""
    mean(x::Vector{Float64}) -> Float64

平均を計算
"""
function mean(x::Vector{Float64})
    return sum(x) / length(x)
end

"""
    run_relaxation_and_extract_sps(params, t_relax, Ntheta) -> StablePeriodicSolution

緩和シミュレーションと周期解抽出のワンショット実行
"""
function run_relaxation_and_extract_sps(params::NetworkParams, t_relax::Float64, Ntheta::Int)
    cmd = ComputeSPS(params, Ntheta, t_relax)
    result = interpret(cmd)

    if is_success(result)
        return unwrap(result)
    else
        throw(ErrorException("Failed to compute SPS: $(result.error)"))
    end
end

"""
    run_psf_computation(sps, n_iterations) -> PhaseSensitivityFunction

PSF計算のワンショット実行
"""
function run_psf_computation(sps::StablePeriodicSolution, n_iterations::Int)
    cmd = ComputePSF(sps, n_iterations)
    result = interpret(cmd)

    if is_success(result)
        return unwrap(result)
    else
        throw(ErrorException("Failed to compute PSF: $(result.error)"))
    end
end

