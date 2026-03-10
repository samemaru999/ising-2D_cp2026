"""
    Logic/PhaseReduction.jl - 位相縮約理論 [Pure]

集団振動ネットワークの位相縮約に関する純粋関数群

主要概念:
- 安定周期解 (SPS): X_s(θ)
- 位相感受関数 (PSF): Q(θ)
- 個別位相結合関数: γ_ijk(φ, ψ)
- 全体位相結合関数: Γ(φ, ψ)

Reference: Kuroiwa & Kato, NOLTA2025
"""

using LinearAlgebra

# =============================================================================
# 個別位相結合関数 γ_ijk
# =============================================================================

"""
    compute_individual_pcf(sps::StablePeriodicSolution, psf::PhaseSensitivityFunction,
                            i::Int, j::Int, k::Int) -> IndividualPCF

個別位相結合関数 γ_ijk(φ, ψ) を計算（Pure関数）

γ_ijk(φ, ψ) = (1/T)∫₀^T Q_i(θ) · H_ijk(X_i(θ), X_j(θ-φ), X_k(θ-ψ)) dθ

ここで:
- Q_i(θ): 振動子iの位相感受関数
- X_i(θ): 安定周期解上の振動子iの状態
- H_ijk: 高次相互作用関数

# Arguments
- `sps::StablePeriodicSolution`: 安定周期解
- `psf::PhaseSensitivityFunction`: 位相感受関数
- `i, j, k::Int`: 振動子インデックス

# Returns
- `IndividualPCF`: γ_ijk(φ, ψ) を含む構造体
"""
function compute_individual_pcf(sps::StablePeriodicSolution, psf::PhaseSensitivityFunction,
                                 i::Int, j::Int, k::Int)
    Ntheta = sps.Ntheta
    T = sps.T
    omega = sps.omega
    N = sps.N

    # γ_ijk(φ, ψ) を Ntheta × Ntheta の行列として計算
    gamma = zeros(Ntheta, Ntheta)

    # 各位相差 (φ, ψ) に対して積分を計算
    for phi_idx in 1:Ntheta
        for psi_idx in 1:Ntheta
            # 積分を計算
            integral = 0.0

            for theta_idx in 1:Ntheta
                # θ での値
                theta = index_to_phase(theta_idx, Ntheta)

                # 位相差
                phi = index_to_phase(phi_idx, Ntheta)
                psi = index_to_phase(psi_idx, Ntheta)

                # X(θ) から v成分を取得
                v_i = sps.Xs[2i, theta_idx]

                # X(θ - φ) から v_j を取得（周期境界条件）
                theta_minus_phi_idx = phase_to_index(theta - phi, Ntheta)
                v_j = sps.Xs[2j, theta_minus_phi_idx]

                # X(θ - ψ) から v_k を取得
                theta_minus_psi_idx = phase_to_index(theta - psi, Ntheta)
                v_k = sps.Xs[2k, theta_minus_psi_idx]

                # H_ijk の v成分（u成分は0）
                H_v = higher_order_H(v_i, v_j, v_k)

                # Q_i(θ) の v成分（u成分との内積を取る）
                Q_u = psf.Q[2i-1, theta_idx]
                Q_v = psf.Q[2i, theta_idx]

                # H = (0, H_v) なので内積は Q_v * H_v
                integrand = Q_v * H_v

                integral += integrand
            end

            # 台形積分の補正と周期での平均化
            gamma[phi_idx, psi_idx] = integral / Ntheta
        end
    end

    return IndividualPCF(gamma, i, j, k)
end

"""
    compute_all_individual_pcfs(sps, psf) -> Vector{IndividualPCF}

全ての (i, j, k) の組み合わせに対して γ_ijk を計算（Pure関数）
"""
function compute_all_individual_pcfs(sps::StablePeriodicSolution,
                                      psf::PhaseSensitivityFunction)
    N = sps.N
    pcfs = IndividualPCF[]

    for i in 1:N
        for j in 1:N
            for k in 1:N
                pcf = compute_individual_pcf(sps, psf, i, j, k)
                push!(pcfs, pcf)
            end
        end
    end

    return pcfs
end

# =============================================================================
# 全体位相結合関数 Γ
# =============================================================================

"""
    compute_total_pcf(individual_pcfs::Vector{IndividualPCF}, C::Array{Float64,3},
                       N::Int, Ntheta::Int) -> TotalPCF

全体位相結合関数 Γ(φ, ψ) を計算（Pure関数）

Γ(φ, ψ) = Σ_{i,j,k} C_{ijk} γ_ijk(φ, ψ)

# Arguments
- `individual_pcfs::Vector{IndividualPCF}`: 個別PCFのリスト
- `C::Array{Float64,3}`: 結合テンソル
- `N::Int`: 振動子数
- `Ntheta::Int`: 位相分割数

# Returns
- `TotalPCF`: Γ(φ, ψ) と導関数を含む構造体
"""
function compute_total_pcf(individual_pcfs::Vector{IndividualPCF}, C::Array{Float64,3},
                            N::Int, Ntheta::Int)
    Gamma = zeros(Ntheta, Ntheta)

    # 個別PCFを結合テンソルで重み付け和
    for pcf in individual_pcfs
        i, j, k = pcf.i, pcf.j, pcf.k
        Gamma .+= C[i, j, k] * pcf.gamma
    end

    # 導関数を計算（数値微分）
    dGamma_dphi = compute_partial_derivative_phi(Gamma, Ntheta)
    dGamma_dpsi = compute_partial_derivative_psi(Gamma, Ntheta)

    return TotalPCF(N, Ntheta, Gamma, C, dGamma_dphi, dGamma_dpsi)
end

"""
    compute_total_pcf_direct(sps, psf, C) -> TotalPCF

個別PCFを経由せずに直接Γを計算（メモリ効率版）（Pure関数）
"""
function compute_total_pcf_direct(sps::StablePeriodicSolution, psf::PhaseSensitivityFunction,
                                   C::Array{Float64,3})
    N = sps.N
    Ntheta = sps.Ntheta
    Gamma = zeros(Ntheta, Ntheta)

    for phi_idx in 1:Ntheta
        for psi_idx in 1:Ntheta
            phi = index_to_phase(phi_idx, Ntheta)
            psi = index_to_phase(psi_idx, Ntheta)

            gamma_sum = 0.0

            # 全ての (i, j, k) に対して計算
            for i in 1:N
                for j in 1:N
                    for k in 1:N
                        if C[i, j, k] == 0.0
                            continue  # 最適化: ゼロ結合はスキップ
                        end

                        integral = 0.0
                        for theta_idx in 1:Ntheta
                            theta = index_to_phase(theta_idx, Ntheta)

                            v_i = sps.Xs[2i, theta_idx]
                            theta_minus_phi_idx = phase_to_index(theta - phi, Ntheta)
                            v_j = sps.Xs[2j, theta_minus_phi_idx]
                            theta_minus_psi_idx = phase_to_index(theta - psi, Ntheta)
                            v_k = sps.Xs[2k, theta_minus_psi_idx]

                            H_v = higher_order_H(v_i, v_j, v_k)
                            Q_v = psf.Q[2i, theta_idx]

                            integral += Q_v * H_v
                        end

                        gamma_sum += C[i, j, k] * integral / Ntheta
                    end
                end
            end

            Gamma[phi_idx, psi_idx] = gamma_sum
        end
    end

    dGamma_dphi = compute_partial_derivative_phi(Gamma, Ntheta)
    dGamma_dpsi = compute_partial_derivative_psi(Gamma, Ntheta)

    return TotalPCF(N, Ntheta, Gamma, C, dGamma_dphi, dGamma_dpsi)
end

# =============================================================================
# 導関数の計算（位相結合関数の勾配）
# =============================================================================

"""
    compute_partial_derivative_phi(Gamma, Ntheta) -> Float64

∂Γ/∂φ を (0, 0) で計算（Pure関数）

中心差分法を使用
"""
function compute_partial_derivative_phi(Gamma::Matrix{Float64}, Ntheta::Int)
    # φ方向の刻み幅
    dphi = 2π / (Ntheta - 1)

    # 中心差分: (Γ(h, 0) - Γ(-h, 0)) / (2h)
    # インデックス1が φ=0 に対応
    # インデックス2が φ=dphi に対応
    # インデックス Ntheta が φ=2π-dphi に対応（周期境界で φ=-dphi と同等）

    Gamma_plus = Gamma[2, 1]        # Γ(dphi, 0)
    Gamma_minus = Gamma[Ntheta, 1]  # Γ(-dphi, 0) ≈ Γ(2π-dphi, 0)

    return (Gamma_plus - Gamma_minus) / (2 * dphi)
end

"""
    compute_partial_derivative_psi(Gamma, Ntheta) -> Float64

∂Γ/∂ψ を (0, 0) で計算（Pure関数）
"""
function compute_partial_derivative_psi(Gamma::Matrix{Float64}, Ntheta::Int)
    dpsi = 2π / (Ntheta - 1)

    Gamma_plus = Gamma[1, 2]        # Γ(0, dpsi)
    Gamma_minus = Gamma[1, Ntheta]  # Γ(0, -dpsi)

    return (Gamma_plus - Gamma_minus) / (2 * dpsi)
end

"""
    compute_gamma_derivatives(individual_pcfs, Ntheta) -> Dict{Tuple{Int,Int,Int}, Tuple{Float64, Float64}}

各γ_ijkの導関数を計算（Pure関数）

γ₁_ijk = ∂γ_ijk/∂φ at (0,0)
γ₂_ijk = ∂γ_ijk/∂ψ at (0,0)
"""
function compute_gamma_derivatives(individual_pcfs::Vector{IndividualPCF}, Ntheta::Int)
    derivatives = Dict{Tuple{Int,Int,Int}, Tuple{Float64, Float64}}()

    for pcf in individual_pcfs
        gamma1 = compute_partial_derivative_phi(pcf.gamma, Ntheta)
        gamma2 = compute_partial_derivative_psi(pcf.gamma, Ntheta)
        derivatives[(pcf.i, pcf.j, pcf.k)] = (gamma1, gamma2)
    end

    return derivatives
end

# =============================================================================
# 安定性指標・回転特性の計算
# =============================================================================

"""
    compute_linear_stability(Gamma1::Float64, Gamma2::Float64) -> Float64

線形安定性指標 Λ を計算（Pure関数）

Λ = -3/2 (Γ₁ + Γ₂)

Λ < 0 のとき同期解は線形安定
"""
function compute_linear_stability(Gamma1::Float64, Gamma2::Float64)
    return -1.5 * (Gamma1 + Gamma2)
end

"""
    compute_rotation_characteristic(Gamma1::Float64, Gamma2::Float64) -> Float64

回転特性 R を計算（Pure関数）

R = (√3/2) |Γ₁ - Γ₂|
"""
function compute_rotation_characteristic(Gamma1::Float64, Gamma2::Float64)
    return (sqrt(3) / 2) * abs(Gamma1 - Gamma2)
end

"""
    compute_stability_metrics(pcf::TotalPCF) -> NamedTuple

安定性指標と回転特性を計算（Pure関数）
"""
function compute_stability_metrics(pcf::TotalPCF)
    Lambda = compute_linear_stability(pcf.dGamma_dphi, pcf.dGamma_dpsi)
    R = compute_rotation_characteristic(pcf.dGamma_dphi, pcf.dGamma_dpsi)

    return (
        Lambda = Lambda,
        R = R,
        Gamma1 = pcf.dGamma_dphi,
        Gamma2 = pcf.dGamma_dpsi,
        is_stable = Lambda < 0
    )
end

# =============================================================================
# 縮約位相方程式
# =============================================================================

"""
    reduced_phase_dynamics(phi1, phi2, pcf::TotalPCF, epsilon) -> (dphi1, dphi2)

縮約位相方程式を計算（Pure関数）

dφ₁/dt = εΓ(φ₁, φ₂) - εΓ(φ₂, φ₁)
dφ₂/dt = εΓ(φ₂, -φ₁-φ₂) - εΓ(-φ₁-φ₂, φ₂)

（簡略化版：同期点近傍での線形近似）
dφ₁/dt ≈ ε(Γ₁ - Γ₂)φ₁ + ε(Γ₂)φ₂
"""
function reduced_phase_dynamics(phi1::Float64, phi2::Float64,
                                  pcf::TotalPCF, epsilon::Float64)
    Ntheta = pcf.Ntheta

    # φ₁, φ₂ から対応するインデックスを取得
    phi1_idx = phase_to_index(phi1, Ntheta)
    phi2_idx = phase_to_index(phi2, Ntheta)

    # Γ(φ₁, φ₂)
    Gamma_12 = pcf.Gamma[phi1_idx, phi2_idx]

    # Γ(φ₂, φ₁)
    Gamma_21 = pcf.Gamma[phi2_idx, phi1_idx]

    # 位相差ダイナミクス
    # 第3ネットワークとの位相差も考慮する必要があるが、
    # 簡略化として2自由度系を考える

    # φ₃ = -(φ₁ + φ₂) (位相和の保存)
    phi3 = -(phi1 + phi2)
    phi3_idx = phase_to_index(normalize_phase_positive(phi3), Ntheta)

    # Γ(φ₂, φ₃)
    Gamma_23 = pcf.Gamma[phi2_idx, phi3_idx]

    # Γ(φ₃, φ₂)
    Gamma_32 = pcf.Gamma[phi3_idx, phi2_idx]

    dphi1 = epsilon * (Gamma_12 - Gamma_21)
    dphi2 = epsilon * (Gamma_23 - Gamma_32)

    return (dphi1, dphi2)
end

"""
    reduced_phase_dynamics_linear(phi1, phi2, Gamma1, Gamma2, epsilon) -> (dphi1, dphi2)

線形化された縮約位相方程式（Pure関数）

同期点 (0, 0) 近傍での線形近似（論文 NOLTA2025 式(11)）:
dφ/dt = εJφ

J = | -2Γ₁ - Γ₂    Γ₁ - Γ₂  |
    |  Γ₂ - Γ₁    -2Γ₂ - Γ₁ |

固有値: λ = -3(Γ₁+Γ₂)/2 ± i(√3/2)|Γ₁-Γ₂|
Λ < 0 (Γ₁ + Γ₂ > 0) のとき収束
"""
function reduced_phase_dynamics_linear(phi1::Float64, phi2::Float64,
                                        Gamma1::Float64, Gamma2::Float64,
                                        epsilon::Float64)
    # 論文 NOLTA2025 式(11) のヤコビアン
    j11 = -2 * Gamma1 - Gamma2
    j12 = Gamma1 - Gamma2
    j21 = Gamma2 - Gamma1
    j22 = -2 * Gamma2 - Gamma1

    dphi1 = epsilon * (j11 * phi1 + j12 * phi2)
    dphi2 = epsilon * (j21 * phi1 + j22 * phi2)

    return (dphi1, dphi2)
end

"""
    compute_linearized_jacobian(Gamma1, Gamma2, epsilon) -> Matrix{Float64}

線形化ヤコビアン行列を計算（Pure関数）（論文 NOLTA2025 式(11)）

εJ where J = | -2Γ₁ - Γ₂    Γ₁ - Γ₂  |
             |  Γ₂ - Γ₁    -2Γ₂ - Γ₁ |

固有値: λ = -3(Γ₁+Γ₂)/2 ± i(√3/2)|Γ₁-Γ₂|
"""
function compute_linearized_jacobian(Gamma1::Float64, Gamma2::Float64, epsilon::Float64)
    J = [
        -2*Gamma1 - Gamma2  Gamma1 - Gamma2
        Gamma2 - Gamma1     -2*Gamma2 - Gamma1
    ]
    return epsilon * J
end

"""
    compute_eigenvalues(Gamma1, Gamma2, epsilon) -> Tuple{ComplexF64, ComplexF64}

線形化ヤコビアンの固有値を計算（Pure関数）（論文 NOLTA2025 式(11)）

λ₁,₂ = ε[-3(Γ₁+Γ₂)/2 ± i(√3/2)|Γ₁-Γ₂|]

Re(λ) = Λ = -3ε(Γ₁+Γ₂)/2  （安定性指標）
Im(λ) = ±εR = ±ε(√3/2)|Γ₁-Γ₂|  （回転特性）
"""
function compute_eigenvalues(Gamma1::Float64, Gamma2::Float64, epsilon::Float64)
    J = compute_linearized_jacobian(Gamma1, Gamma2, epsilon)
    eigenvals = eigvals(J)
    return (eigenvals[1], eigenvals[2])
end

# =============================================================================
# 位相検出（フルダイナミクスから位相を抽出）
# =============================================================================

"""
    detect_phase_from_state(X, sps::StablePeriodicSolution) -> Float64

状態ベクトルから位相を検出（Pure関数）

周期解軌道への最近傍点を見つけて対応する位相を返す
"""
function detect_phase_from_state(X::Vector{Float64}, sps::StablePeriodicSolution)
    Ntheta = sps.Ntheta

    # 最小距離を持つ位相インデックスを探索
    min_dist = Inf
    best_idx = 1

    for k in 1:Ntheta
        Xs_k = sps.Xs[:, k]
        dist = norm(X - Xs_k)
        if dist < min_dist
            min_dist = dist
            best_idx = k
        end
    end

    return index_to_phase(best_idx, Ntheta)
end

"""
    extract_phases_from_trajectory(trajectory, sps) -> Vector{Float64}

軌道から位相時系列を抽出（Pure関数）
"""
function extract_phases_from_trajectory(trajectory::Matrix{Float64},
                                         sps::StablePeriodicSolution)
    n_steps = size(trajectory, 2)
    phases = zeros(n_steps)

    for t in 1:n_steps
        X_t = trajectory[:, t]
        phases[t] = detect_phase_from_state(X_t, sps)
    end

    return phases
end

# =============================================================================
# 補助関数
# =============================================================================

"""
    interpolate_on_orbit(sps::StablePeriodicSolution, theta::Float64) -> Vector{Float64}

位相θにおける周期解の状態を線形補間で取得（Pure関数）
"""
function interpolate_on_orbit(sps::StablePeriodicSolution, theta::Float64)
    Ntheta = sps.Ntheta
    theta_normalized = normalize_phase_positive(theta)

    # 位相から連続インデックスを計算
    idx_continuous = (Ntheta - 1) * theta_normalized / (2π) + 1

    # 補間用の2点のインデックス
    idx_low = floor(Int, idx_continuous)
    idx_high = ceil(Int, idx_continuous)

    # 周期境界条件
    idx_low = mod1(idx_low, Ntheta)
    idx_high = mod1(idx_high, Ntheta)

    # 補間係数
    t = idx_continuous - floor(idx_continuous)

    if idx_low == idx_high
        return sps.Xs[:, idx_low]
    else
        return (1 - t) * sps.Xs[:, idx_low] + t * sps.Xs[:, idx_high]
    end
end

