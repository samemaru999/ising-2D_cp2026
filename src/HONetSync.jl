"""
    HONetSync.jl - Higher-Order Network Synchronization

高次結合を持つ集団振動ネットワーク同期最適化シミュレーション

FDD (Functional Declarative Design) アーキテクチャに基づく実装:
- Domain Layer: 純粋なドメイン型とeDSL
- Logic Layer: 純粋なビジネスロジック
- Interpreters Layer: eDSLの実行（Impure）
- Runtime Layer: IO操作（Impure）

主要機能:
- F-SIM: 安定周期解（SPS）計算
- F-PSF: 位相感受関数（PSF）計算
- F-PCF: 位相結合関数（PCF）計算
- F-OPT: 高次結合テンソル最適化
- F-PHASE: 位相ダイナミクスシミュレーション

Reference: Kuroiwa & Kato, NOLTA2025
"""
module HONetSync

# =============================================================================
# Domain Layer [Pure]
# =============================================================================

# 定数定義
include("Domain/Constants.jl")

# ドメイン型定義
include("Domain/Types.jl")

# eDSL定義
include("Domain/DSL.jl")

# バリデーション
include("Domain/Validation.jl")

# =============================================================================
# Logic Layer [Pure]
# =============================================================================

# FHN振動子力学
include("Logic/FHN.jl")

# 結合関数
include("Logic/Coupling.jl")

# 位相縮約理論
include("Logic/PhaseReduction.jl")

# 最適化アルゴリズム
include("Logic/Optimization.jl")

# 随伴方程式（分離ヤコビアン実装）
include("Logic/AdjointEquation.jl")

# =============================================================================
# Interpreters Layer [Impure]
# =============================================================================

# シミュレーションインタープリター
include("Interpreters/SimulationInterpreter.jl")

# 位相計算インタープリター
include("Interpreters/PhaseInterpreter.jl")

# =============================================================================
# Runtime Layer [Impure]
# =============================================================================

# IO操作
include("Runtime/IO/IOHandler.jl")

# =============================================================================
# Exports
# =============================================================================

# --- 定数 ---
export N_OSC_DEFAULT, N_NET, Ds, Nθ_DEFAULT
export TREL, RMAX, REP, dt_DEFAULT, ε_DEFAULT
export δ_DEFAULT, a_DEFAULT, b_DEFAULT
export I_PACEMAKER, I_DEFAULT
export INTRA_K_4, INTRA_FULL_K10, INTRA_TREE_K10
export HO_INTER_C_UNI_4, HO_LOPT_C_4, HO_ROPT_C_4

# --- ドメイン型 ---
export OscillatorState, NetworkStateVector
export CouplingType, DiffusiveCoupling, GeneralCoupling
export TopologyType, FullCoupling, TreeCoupling
export FHNParams, NetworkParams, ThreeNetworkParams
export StablePeriodicSolution, PhaseSensitivityFunction
export IndividualPCF, TotalPCF
export OptimizationResult, PhaseTimeSeries

# --- eDSLコマンド ---
export ComputationCommand
export ComputeSPS, ComputePSF
export ComputeIndividualPCF, ComputeTotalPCF
export OptimizeLinearStability, OptimizeRotation
export SimulatePhaseDynamics, SimulateFullDynamics
export SaveData, LoadData
export ComputationScript

# --- 結果型 ---
export ComputationResult, Success, Failure
export is_success, unwrap, unwrap_or, map_result, bind_result

# --- バリデーション ---
export ValidationError, ValidationResult
export validate_fhn_params, validate_network_params
export validate_coupling_tensor, validate_epsilon
export validate_sps, validate_psf, validate_pcf
export validate_optimization_result
export is_valid, assert_valid

# --- ユーティリティ関数 ---
export get_intra_K, get_default_external_input
export create_uniform_coupling_tensor
export get_oscillator_state, set_oscillator_state!
export index_to_phase, phase_to_index
export normalize_phase, normalize_phase_positive

# --- FHN力学 ---
export fhn_single_dynamics, fhn_single_jacobian
export network_intrinsic_dynamics, network_full_dynamics
export network_vector_field
export network_intrinsic_jacobian, network_full_jacobian
export three_network_dynamics
export higher_order_interaction

# --- 結合関数 ---
export diffusive_coupling, compute_intra_coupling
export higher_order_H, higher_order_H_full, higher_order_H_derivatives
export compute_interaction_tensor, compute_inter_coupling
export ThreeNetworkState, compute_three_network_coupling
export compute_three_network_dynamics
export create_full_coupling_matrix, create_ring_coupling_matrix
export create_diagonal_coupling_tensor
export compute_coupling_frobenius_norm, is_antisymmetric

# --- 位相縮約 ---
export compute_individual_pcf, compute_all_individual_pcfs
export compute_total_pcf, compute_total_pcf_direct
export compute_partial_derivative_phi, compute_partial_derivative_psi
export compute_gamma_derivatives
export compute_linear_stability, compute_rotation_characteristic
export compute_stability_metrics
export reduced_phase_dynamics, reduced_phase_dynamics_linear
export compute_linearized_jacobian, compute_eigenvalues
export detect_phase_from_state, extract_phases_from_trajectory

# --- 最適化 ---
export compute_gamma_partial_derivatives, compute_all_gamma_derivatives
export optimize_linear_stability, optimize_linear_stability_constrained
export optimize_rotation, optimize_rotation_constrained
export optimize_combined
export evaluate_coupling_performance, compare_coupling_strategies
export compute_sensitivity_to_coupling

# --- 随伴方程式（分離ヤコビアン） ---
export AdjointState, OscillatorJacobians, CouplingJacobians
export compute_intrinsic_jacobian, compute_all_intrinsic_jacobians
export compute_coupling_jacobians
export adjoint_rhs_intrinsic, adjoint_rhs_coupling_self, adjoint_rhs_coupling_cross
export compute_adjoint_rhs, normalize_adjoint_state
export solve_adjoint_equation, compute_psf_separated
export create_adjoint_solver, compose_psf_pipeline
export solve_adjoint_diffusive_v2, compute_psf_v2
export solve_adjoint_diffusive_v3, compute_psf_v3
export compute_psf, PSFVersion, PSF_V2, PSF_V3

# --- インタープリター ---
export interpret
export run_relaxation_and_extract_sps, run_psf_computation
export compute_all_individual_pcfs_parallel
export run_full_pipeline, run_comparison_experiment
export interpret_script
export get_predefined_coupling
export simulate_one_period, sample_periodic_orbit_from_simulation

# --- IO ---
export generate_filename
export save_sps, save_psf, save_total_pcf
export save_optimization_result, save_phase_time_series
export save_pipeline_results
export load_sps, load_psf, load_total_pcf, load_optimization_result
export list_data_files, get_file_info

# =============================================================================
# 高レベルAPI
# =============================================================================

"""
    run_optimization(N::Int; coupling_type, target_value, save_results)

高レベル最適化実行関数

# Arguments
- `N::Int`: 振動子数 (4 or 10)
- `coupling_type::Symbol`: 結合タイプ (:uniform, :linear, :rotation)
- `target_value::Float64`: 最適化目標値
- `save_results::Bool`: 結果を保存するか

# Returns
- パイプライン実行結果

# Example
```julia
results = run_optimization(4; coupling_type=:linear, target_value=-0.1)
```
"""
function run_optimization(N::Int=N_OSC_DEFAULT;
    coupling_type::Symbol=:uniform,
    target_value::Float64=0.0,
    save_results::Bool=false,
    output_dir::String="output")
    # パラメータ作成
    params = NetworkParams(N)

    # パイプライン実行
    results = run_full_pipeline(params;
        coupling_type=coupling_type,
        target_value=target_value)

    # 結果保存
    if save_results
        save_pipeline_results(output_dir, results)
    end

    return results
end

"""
    compare_strategies(N::Int; target_q, target_mu, save_results)

3つの結合戦略を比較

# Example
```julia
comparison = compare_strategies(4; target_q=-0.1, target_mu=0.5)
```
"""
function compare_strategies(N::Int=N_OSC_DEFAULT;
    target_q::Float64=-0.1,
    target_mu::Float64=0.5,
    save_results::Bool=false,
    output_dir::String="output")
    params = NetworkParams(N)

    results = run_comparison_experiment(params;
        target_q=target_q,
        target_mu=target_mu)

    # TODO: 結果保存の実装

    return results
end

"""
    simulate_phase_dynamics(pcf::TotalPCF, phi1_0, phi2_0, epsilon, t_end)

位相ダイナミクスシミュレーション

# Example
```julia
pts = simulate_phase_dynamics(total_pcf, 0.1, -0.1, 2e-4, 10000.0)
```
"""
function simulate_phase_dynamics(pcf::TotalPCF,
    phi1_0::Float64,
    phi2_0::Float64,
    epsilon::Float64,
    t_end::Float64)
    cmd = SimulatePhaseDynamics(pcf, phi1_0, phi2_0, epsilon, 0.0, t_end)
    result = interpret(cmd)

    if is_success(result)
        return unwrap(result)
    else
        error("Simulation failed: $(result.error)")
    end
end

# =============================================================================
# モジュール初期化
# =============================================================================

function __init__()
    @info "HONetSync v2.0.0 loaded"
    @info "FDD Architecture: Domain → Logic → Interpreters → Runtime"
end

end # module HONetSync
