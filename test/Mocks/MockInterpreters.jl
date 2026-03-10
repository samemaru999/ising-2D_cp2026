"""
    Mocks/MockInterpreters.jl - インタープリターのモック実装

ODEシミュレーションを実行せずに即座に結果を返す
テストの高速化と決定論的な結果のために使用
"""

using HONetSync

# MockData.jl の関数を使用
include("MockData.jl")

# =============================================================================
# Mock Interpreter Functions
# =============================================================================

"""
    mock_interpret_sps(cmd::ComputeSPS) -> ComputationResult{StablePeriodicSolution}

モックSPS計算 - ODEを解かずに即座に事前計算データを返す
"""
function mock_interpret_sps(cmd::ComputeSPS)
    try
        sps = create_mock_sps(N=cmd.params.N, Ntheta=cmd.Ntheta)
        return Success(sps)
    catch e
        return Failure{StablePeriodicSolution}("Mock SPS failed: $e")
    end
end

"""
    mock_interpret_psf(cmd::ComputePSF) -> ComputationResult{PhaseSensitivityFunction}

モックPSF計算 - 随伴方程式を解かずに即座に結果を返す
"""
function mock_interpret_psf(cmd::ComputePSF)
    try
        psf = create_mock_psf(cmd.sps)
        return Success(psf)
    catch e
        return Failure{PhaseSensitivityFunction}("Mock PSF failed: $e")
    end
end

"""
    mock_interpret_individual_pcf(cmd::ComputeIndividualPCF) -> ComputationResult{IndividualPCF}

モック個別PCF計算
"""
function mock_interpret_individual_pcf(cmd::ComputeIndividualPCF)
    try
        Ntheta = cmd.sps.Ntheta
        gamma = zeros(Ntheta, Ntheta)

        for phi_idx in 1:Ntheta
            for psi_idx in 1:Ntheta
                phi = 2π * (phi_idx - 1) / (Ntheta - 1)
                psi = 2π * (psi_idx - 1) / (Ntheta - 1)
                gamma[phi_idx, psi_idx] = 0.01 * sin(phi - psi + cmd.k * 0.1)
            end
        end

        pcf = IndividualPCF(gamma, cmd.i, cmd.j, cmd.k)
        return Success(pcf)
    catch e
        return Failure{IndividualPCF}("Mock IndividualPCF failed: $e")
    end
end

"""
    mock_interpret_total_pcf(cmd::ComputeTotalPCF) -> ComputationResult{TotalPCF}

モックTotalPCF計算
"""
function mock_interpret_total_pcf(cmd::ComputeTotalPCF)
    try
        total_pcf = create_mock_total_pcf(N=cmd.N, Ntheta=cmd.Ntheta)
        return Success(total_pcf)
    catch e
        return Failure{TotalPCF}("Mock TotalPCF failed: $e")
    end
end

"""
    mock_interpret_optimize_linear(cmd::OptimizeLinearStability) -> ComputationResult{OptimizationResult}

モック線形安定性最適化
"""
function mock_interpret_optimize_linear(cmd::OptimizeLinearStability)
    try
        N = cmd.N
        C_opt = randn(N, N, N) * 0.1
        achieved = cmd.target_q * 0.98  # 目標の98%達成
        frobenius_norm = sqrt(sum(C_opt .^ 2))

        result = OptimizationResult(N, C_opt, :linear, cmd.target_q, achieved, frobenius_norm)
        return Success(result)
    catch e
        return Failure{OptimizationResult}("Mock linear optimization failed: $e")
    end
end

"""
    mock_interpret_optimize_rotation(cmd::OptimizeRotation) -> ComputationResult{OptimizationResult}

モック回転特性最適化
"""
function mock_interpret_optimize_rotation(cmd::OptimizeRotation)
    try
        N = cmd.N
        C_opt = randn(N, N, N) * 0.1
        achieved = cmd.target_mu * 0.98
        frobenius_norm = sqrt(sum(C_opt .^ 2))

        result = OptimizationResult(N, C_opt, :rotation, cmd.target_mu, achieved, frobenius_norm)
        return Success(result)
    catch e
        return Failure{OptimizationResult}("Mock rotation optimization failed: $e")
    end
end

"""
    mock_interpret_phase_dynamics(cmd::SimulatePhaseDynamics) -> ComputationResult{PhaseTimeSeries}

モック位相シミュレーション - 指数減衰する軌道を返す
"""
function mock_interpret_phase_dynamics(cmd::SimulatePhaseDynamics)
    try
        t_steps = 100
        time = collect(range(cmd.t_start, cmd.t_end, length=t_steps))

        # 安定性に基づいて減衰率を決定
        Lambda = compute_linear_stability(cmd.total_pcf.dGamma_dphi, cmd.total_pcf.dGamma_dpsi)
        decay_rate = Lambda < 0 ? abs(Lambda) * cmd.epsilon : -Lambda * cmd.epsilon

        phases = zeros(N_NET, t_steps)
        phase_diffs = zeros(2, t_steps)

        for (t_idx, t) in enumerate(time)
            dt = t - cmd.t_start
            phase_diffs[1, t_idx] = cmd.initial_phi1 * exp(decay_rate * dt)
            phase_diffs[2, t_idx] = cmd.initial_phi2 * exp(decay_rate * dt)
            phases[2, t_idx] = phase_diffs[1, t_idx]
            phases[3, t_idx] = phase_diffs[1, t_idx] + phase_diffs[2, t_idx]
        end

        pts = PhaseTimeSeries(time, phases, phase_diffs)
        return Success(pts)
    catch e
        return Failure{PhaseTimeSeries}("Mock phase dynamics failed: $e")
    end
end

"""
    mock_interpret_full_dynamics(cmd::SimulateFullDynamics) -> ComputationResult{PhaseTimeSeries}

モック元力学系シミュレーション
"""
function mock_interpret_full_dynamics(cmd::SimulateFullDynamics)
    try
        pts = create_mock_phase_time_series(t_end=cmd.t_end - cmd.t_start)
        return Success(pts)
    catch e
        return Failure{PhaseTimeSeries}("Mock full dynamics failed: $e")
    end
end

# =============================================================================
# Mock Interpreter Dispatcher
# =============================================================================

"""
    mock_interpret(cmd::ComputationCommand) -> ComputationResult

汎用モックインタープリター - コマンド型に応じて適切なモック関数を呼び出す
"""
function mock_interpret(cmd::ComputeSPS)
    return mock_interpret_sps(cmd)
end

function mock_interpret(cmd::ComputePSF)
    return mock_interpret_psf(cmd)
end

function mock_interpret(cmd::ComputeIndividualPCF)
    return mock_interpret_individual_pcf(cmd)
end

function mock_interpret(cmd::ComputeTotalPCF)
    return mock_interpret_total_pcf(cmd)
end

function mock_interpret(cmd::OptimizeLinearStability)
    return mock_interpret_optimize_linear(cmd)
end

function mock_interpret(cmd::OptimizeRotation)
    return mock_interpret_optimize_rotation(cmd)
end

function mock_interpret(cmd::SimulatePhaseDynamics)
    return mock_interpret_phase_dynamics(cmd)
end

function mock_interpret(cmd::SimulateFullDynamics)
    return mock_interpret_full_dynamics(cmd)
end

# =============================================================================
# Mock Service Factory
# =============================================================================

"""
    MockSimulationService

モックシミュレーションサービス構造体
"""
struct MockSimulationService
    interpret_sps::Function
    interpret_psf::Function
    interpret_pcf::Function
    interpret_optimization::Function
    interpret_simulation::Function
end

"""
    create_mock_simulation_service() -> MockSimulationService

モックシミュレーションサービスを作成（依存性注入用）
"""
function create_mock_simulation_service()
    return MockSimulationService(
        mock_interpret_sps,
        mock_interpret_psf,
        mock_interpret_total_pcf,
        mock_interpret_optimize_linear,
        mock_interpret_phase_dynamics
    )
end

# =============================================================================
# Failure Injection (テスト用)
# =============================================================================

"""
    create_failing_interpreter(fail_at::Symbol) -> Function

指定したコマンドで失敗するインタープリターを作成

テストで失敗ケースを検証するために使用
"""
function create_failing_interpreter(fail_at::Symbol)
    function failing_interpret(cmd::ComputationCommand)
        cmd_type = Symbol(typeof(cmd).name.name)
        if cmd_type == fail_at
            return Failure{Any}("Intentional failure at $fail_at")
        else
            return mock_interpret(cmd)
        end
    end
    return failing_interpret
end

"""
    create_slow_interpreter(delay_seconds::Float64) -> Function

意図的に遅いインタープリターを作成

タイムアウトやキャンセル処理のテストに使用
"""
function create_slow_interpreter(delay_seconds::Float64)
    function slow_interpret(cmd::ComputationCommand)
        sleep(delay_seconds)
        return mock_interpret(cmd)
    end
    return slow_interpret
end

