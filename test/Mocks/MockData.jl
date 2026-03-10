"""
    Mocks/MockData.jl - テスト用モックデータ

計算済みのSPS、PSF、PCFデータを提供
（実際のODEシミュレーションなしでテスト可能）
"""

using HONetSync
using Random

# =============================================================================
# Mock SPS (N=4, Ntheta=11 の小さなテスト用)
# =============================================================================

"""
    create_mock_sps(; N, Ntheta) -> StablePeriodicSolution

N=4, Ntheta=11 のモック安定周期解を生成

単純な正弦波で近似した周期軌道を返す
"""
function create_mock_sps(; N::Int=4, Ntheta::Int=11)
    T = 50.0
    omega = 2π / T

    # 単純な正弦波で近似した周期軌道
    Xs = zeros(Ds * N, Ntheta)
    for k in 1:Ntheta
        theta = 2π * (k - 1) / (Ntheta - 1)
        for i in 1:N
            phase_offset = (i - 1) * π / N
            Xs[2i-1, k] = 0.5 * sin(theta + phase_offset)      # u
            Xs[2i, k] = 1.5 * cos(theta + phase_offset) - 0.5  # v
        end
    end
    Xs[:, end] = Xs[:, 1]  # 周期性

    params = NetworkParams(N)
    return StablePeriodicSolution(N, Ntheta, T, omega, Xs, params)
end

# =============================================================================
# Mock PSF
# =============================================================================

"""
    create_mock_psf(sps) -> PhaseSensitivityFunction

モック位相感受関数を生成

実際の随伴方程式解ではなく、単純化した関数を返す
"""
function create_mock_psf(sps::StablePeriodicSolution)
    N = sps.N
    Ntheta = sps.Ntheta
    T = sps.T
    omega = sps.omega

    # 単純化したPSF（実際の随伴方程式解ではない）
    Q = zeros(Ds * N, Ntheta)
    for k in 1:Ntheta
        theta = 2π * (k - 1) / (Ntheta - 1)
        for i in 1:N
            Q[2i-1, k] = 0.1 * cos(theta)
            Q[2i, k] = omega * sin(theta)  # 正規化条件を近似
        end
    end
    Q[:, end] = Q[:, 1]

    return PhaseSensitivityFunction(N, Ntheta, Q, T, omega)
end

# =============================================================================
# Mock Individual PCFs
# =============================================================================

"""
    create_mock_individual_pcfs(N, Ntheta) -> Vector{IndividualPCF}

モック個別位相結合関数を生成

全ての (i, j, k) 組み合わせに対するγ_ijkを返す
"""
function create_mock_individual_pcfs(N::Int, Ntheta::Int)
    pcfs = IndividualPCF[]

    for i in 1:N
        for j in 1:N
            for k in 1:N
                # 単純な周期関数で近似
                gamma = zeros(Ntheta, Ntheta)
                for phi_idx in 1:Ntheta
                    for psi_idx in 1:Ntheta
                        phi = 2π * (phi_idx - 1) / (Ntheta - 1)
                        psi = 2π * (psi_idx - 1) / (Ntheta - 1)
                        # i==jのとき強い結合、それ以外は弱い結合
                        strength = (i == j) ? 1.0 : 0.3
                        gamma[phi_idx, psi_idx] = 0.01 * strength * sin(phi - psi + k * π / N)
                    end
                end
                push!(pcfs, IndividualPCF(gamma, i, j, k))
            end
        end
    end

    return pcfs
end

# =============================================================================
# Mock Total PCF
# =============================================================================

"""
    create_mock_total_pcf(; N, Ntheta, Gamma1, Gamma2) -> TotalPCF

モック全体位相結合関数を生成

指定したΓ₁, Γ₂の導関数値を持つように構成
"""
function create_mock_total_pcf(; N::Int=4, Ntheta::Int=11,
                                 Gamma1::Float64=0.1, Gamma2::Float64=0.05)
    Gamma = zeros(Ntheta, Ntheta)
    dphi = 2π / (Ntheta - 1)

    for phi_idx in 1:Ntheta
        for psi_idx in 1:Ntheta
            phi = dphi * (phi_idx - 1)
            psi = dphi * (psi_idx - 1)
            # 指定したΓ₁, Γ₂を持つように構成
            Gamma[phi_idx, psi_idx] = Gamma1 * phi + Gamma2 * psi +
                                       0.01 * sin(phi) * sin(psi)
        end
    end
    Gamma[1, 1] = 0.0  # 同期点での値は0

    C = create_uniform_coupling_tensor(N)
    return TotalPCF(N, Ntheta, Gamma, C, Gamma1, Gamma2)
end

"""
    create_mock_total_pcf_stable(; N, Ntheta) -> TotalPCF

安定な同期を示すモックTotalPCFを生成

Λ < 0 となるようにΓ₁ + Γ₂ > 0 を設定
"""
function create_mock_total_pcf_stable(; N::Int=4, Ntheta::Int=11)
    return create_mock_total_pcf(N=N, Ntheta=Ntheta, Gamma1=0.15, Gamma2=0.1)
end

"""
    create_mock_total_pcf_unstable(; N, Ntheta) -> TotalPCF

不安定な同期を示すモックTotalPCFを生成

Λ > 0 となるようにΓ₁ + Γ₂ < 0 を設定
"""
function create_mock_total_pcf_unstable(; N::Int=4, Ntheta::Int=11)
    return create_mock_total_pcf(N=N, Ntheta=Ntheta, Gamma1=-0.1, Gamma2=-0.05)
end

# =============================================================================
# Mock Optimization Result
# =============================================================================

"""
    create_mock_optimization_result(; N, opt_type) -> OptimizationResult

モック最適化結果を生成
"""
function create_mock_optimization_result(; N::Int=4, opt_type::Symbol=:linear)
    C_opt = randn(N, N, N) * 0.1
    target = opt_type == :linear ? -0.1 : 0.5
    achieved = target * (1 + 0.01 * randn())
    frobenius_norm = sqrt(sum(C_opt .^ 2))

    return OptimizationResult(N, C_opt, opt_type, target, achieved, frobenius_norm)
end

# =============================================================================
# Mock Phase Time Series
# =============================================================================

"""
    create_mock_phase_time_series(; t_end, dt, decay_rate) -> PhaseTimeSeries

モック位相時系列を生成

指数減衰する位相差の軌道を返す
"""
function create_mock_phase_time_series(; t_end::Float64=1000.0,
                                         dt::Float64=1.0,
                                         initial_phi1::Float64=0.5,
                                         initial_phi2::Float64=-0.3,
                                         decay_rate::Float64=0.01)
    time = collect(0.0:dt:t_end)
    n_steps = length(time)

    phases = zeros(N_NET, n_steps)
    phase_diffs = zeros(2, n_steps)

    for (t_idx, t) in enumerate(time)
        # 指数減衰
        phi1 = initial_phi1 * exp(-decay_rate * t)
        phi2 = initial_phi2 * exp(-decay_rate * t)

        phase_diffs[1, t_idx] = phi1
        phase_diffs[2, t_idx] = phi2

        # 絶対位相（ネットワークAを基準）
        phases[1, t_idx] = 0.0
        phases[2, t_idx] = phi1
        phases[3, t_idx] = phi1 + phi2
    end

    return PhaseTimeSeries(time, phases, phase_diffs)
end

# =============================================================================
# Deterministic Mock Data (再現可能なテスト用)
# =============================================================================

"""
    create_deterministic_mock_data(seed) -> NamedTuple

決定論的なシード付きモックデータセットを生成

テストの再現性を保証するために使用
"""
function create_deterministic_mock_data(seed::Int=42)
    Random.seed!(seed)

    N = 4
    Ntheta = 11

    sps = create_mock_sps(N=N, Ntheta=Ntheta)
    psf = create_mock_psf(sps)
    individual_pcfs = create_mock_individual_pcfs(N, Ntheta)
    total_pcf = create_mock_total_pcf(N=N, Ntheta=Ntheta)
    optimization_result = create_mock_optimization_result(N=N)
    phase_time_series = create_mock_phase_time_series()

    return (
        N = N,
        Ntheta = Ntheta,
        sps = sps,
        psf = psf,
        individual_pcfs = individual_pcfs,
        total_pcf = total_pcf,
        optimization_result = optimization_result,
        phase_time_series = phase_time_series
    )
end

# =============================================================================
# Edge Case Mock Data
# =============================================================================

"""
    create_edge_case_mock_data() -> NamedTuple

エッジケーステスト用のモックデータを生成

- 非常に小さい/大きい値
- ゼロに近い値
- 境界条件
"""
function create_edge_case_mock_data()
    N = 2  # 最小構成
    Ntheta = 5  # 最小分割

    # ゼロに近い値を持つSPS
    Xs_near_zero = fill(1e-10, Ds * N, Ntheta)
    Xs_near_zero[:, end] = Xs_near_zero[:, 1]

    # 大きな値を持つPCF
    Gamma_large = fill(100.0, Ntheta, Ntheta)
    Gamma_large[1, 1] = 0.0

    return (
        N = N,
        Ntheta = Ntheta,
        near_zero_Xs = Xs_near_zero,
        large_Gamma = Gamma_large
    )
end

