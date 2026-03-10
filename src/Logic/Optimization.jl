"""
    Logic/Optimization.jl - 結合最適化 [Pure]

高次結合テンソルの解析的最適化を行う純粋関数群

最適化目標:
1. 線形安定性: Λ = -3/2(Γ₁ + Γ₂) を最小化（最も負に）
2. 回転特性: R = √3/2|Γ₁ - Γ₂| を最大化

解析解（論文 eq. 18, 21）:
- 線形安定性最適化: C_ijk = q(γ₁_ijk + γ₂_ijk) / Σ(γ₁_ijk + γ₂_ijk)²
- 回転特性最適化: C_ijk = μ(γ₁_ijk - γ₂_ijk) / Σ(γ₁_ijk - γ₂_ijk)²

Reference: Kuroiwa & Kato, NOLTA2025
"""

using LinearAlgebra

# =============================================================================
# 導関数の計算
# =============================================================================

"""
    compute_gamma_partial_derivatives(pcf::IndividualPCF, Ntheta::Int) -> (γ₁, γ₂)

個別PCFの偏微分を同期点(0,0)で計算（Pure関数）

γ₁_ijk = ∂γ_ijk/∂φ at (0,0)
γ₂_ijk = ∂γ_ijk/∂ψ at (0,0)
"""
function compute_gamma_partial_derivatives(pcf::IndividualPCF, Ntheta::Int)
    gamma = pcf.gamma
    dphi = 2π / (Ntheta - 1)

    # 中心差分法
    # φ方向: (γ(dphi, 0) - γ(-dphi, 0)) / (2*dphi)
    gamma1 = (gamma[2, 1] - gamma[Ntheta, 1]) / (2 * dphi)

    # ψ方向: (γ(0, dpsi) - γ(0, -dpsi)) / (2*dpsi)
    gamma2 = (gamma[1, 2] - gamma[1, Ntheta]) / (2 * dphi)

    return (gamma1, gamma2)
end

"""
    compute_all_gamma_derivatives(individual_pcfs, Ntheta) -> (γ₁_tensor, γ₂_tensor)

全ての個別PCFの導関数をテンソル形式で計算（Pure関数）

# Returns
- `gamma1_tensor::Array{Float64,3}`: γ₁_ijk (N×N×N)
- `gamma2_tensor::Array{Float64,3}`: γ₂_ijk (N×N×N)
"""
function compute_all_gamma_derivatives(individual_pcfs::Vector{IndividualPCF}, N::Int, Ntheta::Int)
    gamma1_tensor = zeros(N, N, N)
    gamma2_tensor = zeros(N, N, N)

    for pcf in individual_pcfs
        gamma1, gamma2 = compute_gamma_partial_derivatives(pcf, Ntheta)
        gamma1_tensor[pcf.i, pcf.j, pcf.k] = gamma1
        gamma2_tensor[pcf.i, pcf.j, pcf.k] = gamma2
    end

    return (gamma1_tensor, gamma2_tensor)
end

# =============================================================================
# 線形安定性最適化
# =============================================================================

"""
    optimize_linear_stability(individual_pcfs, target_q, N, Ntheta) -> OptimizationResult

線形安定性を最適化する結合テンソルを計算（Pure関数）

目標: Γ₁ + Γ₂ = q を達成

解析解（論文 eq. 18）:
C_ijk = q * (γ₁_ijk + γ₂_ijk) / Σ_{i,j,k} (γ₁_ijk + γ₂_ijk)²

# Arguments
- `individual_pcfs::Vector{IndividualPCF}`: 個別位相結合関数
- `target_q::Float64`: 目標値 q（負の値で安定化）
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数

# Returns
- `OptimizationResult`: 最適化された結合テンソルと達成値
"""
function optimize_linear_stability(individual_pcfs::Vector{IndividualPCF},
                                    target_q::Float64, N::Int, Ntheta::Int)
    # 全導関数を計算
    gamma1_tensor, gamma2_tensor = compute_all_gamma_derivatives(individual_pcfs, N, Ntheta)

    # γ₁ + γ₂ の和テンソル
    gamma_sum = gamma1_tensor .+ gamma2_tensor

    # 正規化定数: Σ(γ₁ + γ₂)²
    norm_constant = sum(gamma_sum .^ 2)

    if norm_constant < 1e-10
        # 正規化定数がほぼゼロの場合
        @warn "Normalization constant is near zero, optimization may be ill-conditioned"
        norm_constant = 1e-10
    end

    # 最適結合テンソル
    C_opt = target_q * gamma_sum / norm_constant

    # 達成値を計算
    # Γ₁ + Γ₂ = Σ C_ijk (γ₁_ijk + γ₂_ijk)
    achieved_q = sum(C_opt .* gamma_sum)

    # フロベニウスノルム
    frobenius_norm = sqrt(sum(C_opt .^ 2))

    return OptimizationResult(
        N,
        C_opt,
        :linear,
        target_q,
        achieved_q,
        frobenius_norm
    )
end

"""
    optimize_linear_stability_constrained(individual_pcfs, target_q, norm_constraint, N, Ntheta)
        -> OptimizationResult

ノルム制約付き線形安定性最適化（Pure関数）

||C||_F = norm_constraint という制約のもとで
Γ₁ + Γ₂ を最大化（負方向に）
"""
function optimize_linear_stability_constrained(individual_pcfs::Vector{IndividualPCF},
                                                target_q::Float64,
                                                norm_constraint::Float64,
                                                N::Int, Ntheta::Int)
    gamma1_tensor, gamma2_tensor = compute_all_gamma_derivatives(individual_pcfs, N, Ntheta)
    gamma_sum = gamma1_tensor .+ gamma2_tensor

    # 正規化して指定ノルムにスケール
    current_norm = sqrt(sum(gamma_sum .^ 2))

    if current_norm < 1e-10
        @warn "Gradient tensor norm is near zero"
        C_opt = zeros(N, N, N)
    else
        # 方向を保ちつつノルム制約を満たす
        direction = gamma_sum / current_norm
        C_opt = norm_constraint * direction * sign(target_q)
    end

    achieved_q = sum(C_opt .* gamma_sum)
    frobenius_norm = sqrt(sum(C_opt .^ 2))

    return OptimizationResult(
        N,
        C_opt,
        :linear,
        target_q,
        achieved_q,
        frobenius_norm
    )
end

# =============================================================================
# 回転特性最適化
# =============================================================================

"""
    optimize_rotation(individual_pcfs, target_mu, N, Ntheta) -> OptimizationResult

回転特性を最適化する結合テンソルを計算（Pure関数）

目標: |Γ₁ - Γ₂| = μ を達成

解析解（論文 eq. 21）:
C_ijk = μ * (γ₁_ijk - γ₂_ijk) / Σ_{i,j,k} (γ₁_ijk - γ₂_ijk)²

# Arguments
- `individual_pcfs::Vector{IndividualPCF}`: 個別位相結合関数
- `target_mu::Float64`: 目標値 μ（非負）
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数

# Returns
- `OptimizationResult`: 最適化された結合テンソルと達成値
"""
function optimize_rotation(individual_pcfs::Vector{IndividualPCF},
                            target_mu::Float64, N::Int, Ntheta::Int)
    # 全導関数を計算
    gamma1_tensor, gamma2_tensor = compute_all_gamma_derivatives(individual_pcfs, N, Ntheta)

    # γ₁ - γ₂ の差テンソル
    gamma_diff = gamma1_tensor .- gamma2_tensor

    # 正規化定数: Σ(γ₁ - γ₂)²
    norm_constant = sum(gamma_diff .^ 2)

    if norm_constant < 1e-10
        @warn "Normalization constant is near zero, optimization may be ill-conditioned"
        norm_constant = 1e-10
    end

    # 最適結合テンソル
    C_opt = target_mu * gamma_diff / norm_constant

    # 達成値を計算
    # Γ₁ - Γ₂ = Σ C_ijk (γ₁_ijk - γ₂_ijk)
    achieved_diff = sum(C_opt .* gamma_diff)
    achieved_mu = abs(achieved_diff)

    # フロベニウスノルム
    frobenius_norm = sqrt(sum(C_opt .^ 2))

    return OptimizationResult(
        N,
        C_opt,
        :rotation,
        target_mu,
        achieved_mu,
        frobenius_norm
    )
end

"""
    optimize_rotation_constrained(individual_pcfs, target_mu, norm_constraint, N, Ntheta)
        -> OptimizationResult

ノルム制約付き回転特性最適化（Pure関数）
"""
function optimize_rotation_constrained(individual_pcfs::Vector{IndividualPCF},
                                        target_mu::Float64,
                                        norm_constraint::Float64,
                                        N::Int, Ntheta::Int)
    gamma1_tensor, gamma2_tensor = compute_all_gamma_derivatives(individual_pcfs, N, Ntheta)
    gamma_diff = gamma1_tensor .- gamma2_tensor

    current_norm = sqrt(sum(gamma_diff .^ 2))

    if current_norm < 1e-10
        @warn "Gradient tensor norm is near zero"
        C_opt = zeros(N, N, N)
    else
        direction = gamma_diff / current_norm
        C_opt = norm_constraint * direction
    end

    achieved_diff = sum(C_opt .* gamma_diff)
    achieved_mu = abs(achieved_diff)
    frobenius_norm = sqrt(sum(C_opt .^ 2))

    return OptimizationResult(
        N,
        C_opt,
        :rotation,
        target_mu,
        achieved_mu,
        frobenius_norm
    )
end

# =============================================================================
# 複合最適化（線形安定性 + 回転特性）
# =============================================================================

"""
    optimize_combined(individual_pcfs, target_q, target_mu, alpha, N, Ntheta)
        -> OptimizationResult

線形安定性と回転特性の複合最適化（Pure関数）

目標関数: α * (Γ₁ + Γ₂) + (1-α) * |Γ₁ - Γ₂|

# Arguments
- `alpha::Float64`: 線形安定性の重み (0 ≤ α ≤ 1)
"""
function optimize_combined(individual_pcfs::Vector{IndividualPCF},
                            target_q::Float64, target_mu::Float64,
                            alpha::Float64, N::Int, Ntheta::Int)
    @assert 0 <= alpha <= 1 "alpha must be in [0, 1]"

    gamma1_tensor, gamma2_tensor = compute_all_gamma_derivatives(individual_pcfs, N, Ntheta)

    gamma_sum = gamma1_tensor .+ gamma2_tensor
    gamma_diff = gamma1_tensor .- gamma2_tensor

    # 重み付き勾配
    gradient = alpha * gamma_sum .+ (1 - alpha) * gamma_diff

    norm_constant = sum(gradient .^ 2)
    if norm_constant < 1e-10
        norm_constant = 1e-10
    end

    # スケールファクター（目標値に基づく）
    scale = alpha * target_q + (1 - alpha) * target_mu

    C_opt = scale * gradient / norm_constant

    achieved_q = sum(C_opt .* gamma_sum)
    achieved_mu = abs(sum(C_opt .* gamma_diff))
    frobenius_norm = sqrt(sum(C_opt .^ 2))

    # 複合最適化では線形安定性を主指標として返す
    return OptimizationResult(
        N,
        C_opt,
        :linear,  # 主最適化タイプ
        target_q,
        achieved_q,
        frobenius_norm
    )
end

# =============================================================================
# 最適化ユーティリティ
# =============================================================================

"""
    evaluate_coupling_performance(C, individual_pcfs, N, Ntheta) -> NamedTuple

結合テンソルの性能を評価（Pure関数）

# Returns
- `Gamma1::Float64`: Γ₁
- `Gamma2::Float64`: Γ₂
- `Lambda::Float64`: 線形安定性 Λ = -3/2(Γ₁ + Γ₂)
- `R::Float64`: 回転特性 R = √3/2|Γ₁ - Γ₂|
- `frobenius_norm::Float64`: ||C||_F
"""
function evaluate_coupling_performance(C::Array{Float64,3},
                                        individual_pcfs::Vector{IndividualPCF},
                                        N::Int, Ntheta::Int)
    gamma1_tensor, gamma2_tensor = compute_all_gamma_derivatives(individual_pcfs, N, Ntheta)

    Gamma1 = sum(C .* gamma1_tensor)
    Gamma2 = sum(C .* gamma2_tensor)

    Lambda = -1.5 * (Gamma1 + Gamma2)
    R = (sqrt(3) / 2) * abs(Gamma1 - Gamma2)
    frobenius_norm = sqrt(sum(C .^ 2))

    return (
        Gamma1 = Gamma1,
        Gamma2 = Gamma2,
        Lambda = Lambda,
        R = R,
        frobenius_norm = frobenius_norm,
        is_stable = Lambda < 0
    )
end

"""
    compare_coupling_strategies(individual_pcfs, N, Ntheta;
                                 uniform_strength=0.0112, target_q=-0.1, target_mu=0.5)
        -> Dict{Symbol, NamedTuple}

異なる結合戦略の比較（Pure関数）

# Strategies
- :uniform: 一様結合
- :linear_opt: 線形安定性最適化
- :rotation_opt: 回転特性最適化
"""
function compare_coupling_strategies(individual_pcfs::Vector{IndividualPCF},
                                      N::Int, Ntheta::Int;
                                      uniform_strength::Float64=0.0112,
                                      target_q::Float64=-0.1,
                                      target_mu::Float64=0.5)
    results = Dict{Symbol, NamedTuple}()

    # 一様結合
    C_uniform = create_diagonal_coupling_tensor(N; strength=uniform_strength)
    results[:uniform] = evaluate_coupling_performance(C_uniform, individual_pcfs, N, Ntheta)

    # 線形安定性最適化
    opt_linear = optimize_linear_stability(individual_pcfs, target_q, N, Ntheta)
    results[:linear_opt] = evaluate_coupling_performance(opt_linear.C_opt, individual_pcfs, N, Ntheta)

    # 回転特性最適化
    opt_rotation = optimize_rotation(individual_pcfs, target_mu, N, Ntheta)
    results[:rotation_opt] = evaluate_coupling_performance(opt_rotation.C_opt, individual_pcfs, N, Ntheta)

    return results
end

# =============================================================================
# 感度解析
# =============================================================================

"""
    compute_sensitivity_to_coupling(individual_pcfs, base_C, perturbation_scale, N, Ntheta)
        -> Matrix{Float64}

結合テンソルの各成分に対する感度を計算（Pure関数）

∂Λ/∂C_ijk および ∂R/∂C_ijk を返す
"""
function compute_sensitivity_to_coupling(individual_pcfs::Vector{IndividualPCF},
                                          base_C::Array{Float64,3},
                                          N::Int, Ntheta::Int)
    gamma1_tensor, gamma2_tensor = compute_all_gamma_derivatives(individual_pcfs, N, Ntheta)

    # Λ = -3/2(Γ₁ + Γ₂) = -3/2 Σ C_ijk(γ₁_ijk + γ₂_ijk)
    # ∂Λ/∂C_ijk = -3/2 (γ₁_ijk + γ₂_ijk)
    dLambda_dC = -1.5 * (gamma1_tensor .+ gamma2_tensor)

    # R = √3/2 |Γ₁ - Γ₂|
    # ∂R/∂C_ijk = √3/2 * sign(Γ₁ - Γ₂) * (γ₁_ijk - γ₂_ijk)
    Gamma1 = sum(base_C .* gamma1_tensor)
    Gamma2 = sum(base_C .* gamma2_tensor)
    sign_diff = sign(Gamma1 - Gamma2)
    dR_dC = (sqrt(3) / 2) * sign_diff * (gamma1_tensor .- gamma2_tensor)

    return (dLambda_dC = dLambda_dC, dR_dC = dR_dC)
end

"""
    find_optimal_q_for_stability(individual_pcfs, min_Lambda, N, Ntheta) -> Float64

所望の安定性 Λ_target < 0 を達成するための q 値を計算（Pure関数）
"""
function find_optimal_q_for_stability(individual_pcfs::Vector{IndividualPCF},
                                       target_Lambda::Float64, N::Int, Ntheta::Int)
    # Λ = -3/2(Γ₁ + Γ₂) = -3/2 * q (最適化結合を使用した場合)
    # よって q = -2/3 * Λ_target
    return -2/3 * target_Lambda
end

"""
    find_optimal_mu_for_rotation(target_R) -> Float64

所望の回転特性 R_target を達成するための μ 値を計算（Pure関数）
"""
function find_optimal_mu_for_rotation(target_R::Float64)
    # R = √3/2 |Γ₁ - Γ₂| = √3/2 * μ (最適化結合を使用した場合)
    # よって μ = 2/√3 * R_target
    return 2 / sqrt(3) * target_R
end

