"""
    Logic/AdjointEquation.jl - 随伴方程式の純粋関数実装 [Pure]

式(A1)に基づく随伴方程式の分離実装:
    ω dQ_i/dθ = -J_i(θ)^T Q_i(θ) - Σ_j M_ij(θ)^T Q_i(θ) - Σ_j N_ij(θ) Q_j(θ)

ここで:
    J_i = ∂f_i/∂x_i        (固有ダイナミクスのヤコビアン)
    M_ij = ∂g_ij/∂x_i      (結合の自己項ヤコビアン)
    N_ij = ∂g_ij/∂x_j      (結合の交差項ヤコビアン)

FDD: すべての関数は純粋（副作用なし、参照透明性あり）
"""

using StaticArrays
using LinearAlgebra

# =============================================================================
# 型定義
# =============================================================================

"""
    AdjointState{N}

N振動子系の随伴変数の状態を表す不変データ型
"""
struct AdjointState{N}
    Q_u::SVector{N,Float64}  # u成分の随伴変数
    Q_v::SVector{N,Float64}  # v成分の随伴変数
end

"""
    OscillatorJacobians

各振動子のヤコビアン成分（2x2行列の4成分）
"""
struct OscillatorJacobians
    J11::Float64  # ∂(du)/∂u
    J12::Float64  # ∂(du)/∂v
    J21::Float64  # ∂(dv)/∂u
    J22::Float64  # ∂(dv)/∂v
end

"""
    CouplingJacobians{N}

結合ヤコビアン M_ij と N_ij
"""
struct CouplingJacobians{N}
    M_self::SVector{N,Float64}   # Σ_j M_ij (対角成分、v成分のみ非ゼロ)
    N_cross::SMatrix{N,N,Float64}  # N_ij 行列 (v成分のみ非ゼロ)
end

"""
    PSFVersion{V}

PSF計算のバージョンを型パラメータで表現するファントム型

使用例:
    compute_psf(sps, params, Val(2))  # V2版（旧Juliaコード準拠）
    compute_psf(sps, params, Val(3))  # V3版（Cコード準拠、dt=0.001）
"""
struct PSFVersion{V} end

# 便利なエイリアス
const PSF_V2 = Val{2}()
const PSF_V3 = Val{3}()

# =============================================================================
# 固有ダイナミクスのヤコビアン J_i (Pure)
# =============================================================================

"""
    compute_intrinsic_jacobian(u::Float64, v::Float64, fhn::FHNParams) -> OscillatorJacobians

単一FHN振動子の固有ヤコビアン J_i = ∂f_i/∂x_i を計算

FHN方程式:
    du/dt = δ(a + v - bu)
    dv/dt = v - v³/3 - u + I

ヤコビアン:
    J = | -δb    δ   |
        | -1   1-v²  |
"""
compute_intrinsic_jacobian(u::Float64, v::Float64, fhn::FHNParams) =
    OscillatorJacobians(
        -fhn.delta * fhn.b,  # J11 = ∂(du)/∂u = -δb
        fhn.delta,           # J12 = ∂(du)/∂v = δ
        -1.0,                # J21 = ∂(dv)/∂u = -1
        1.0 - v^2            # J22 = ∂(dv)/∂v = 1 - v²
    )

"""
    compute_all_intrinsic_jacobians(X::Vector{Float64}, fhn::FHNParams, N::Int)
        -> Vector{OscillatorJacobians}

全振動子の固有ヤコビアンを計算（マップ操作）
"""
compute_all_intrinsic_jacobians(X::Vector{Float64}, fhn::FHNParams, N::Int) =
    [compute_intrinsic_jacobian(X[2i-1], X[2i], fhn) for i in 1:N]

# =============================================================================
# 結合ヤコビアン M_ij, N_ij (Pure)
# =============================================================================

"""
    compute_coupling_jacobians(K::Matrix{Float64}) -> CouplingJacobians

拡散結合 g_ij = K_ij(v_j - v_i) のヤコビアンを計算

M_ij = ∂g_ij/∂x_i = [0, -K_ij]^T  (v_i への寄与)
N_ij = ∂g_ij/∂x_j = [0, +K_ij]^T  (v_j からの寄与)

# Returns
- `CouplingJacobians{N}`: M_self (対角項の和) と N_cross (交差項行列)
"""
function compute_coupling_jacobians(K::Matrix{Float64})
    N = size(K, 1)

    # M_self[i] = -Σ_j K_ij (自己項の和、v成分のみ)
    M_self = SVector{N}([-sum(K[i, :]) for i in 1:N])

    # N_cross[i,j] = K_ij (交差項、v成分のみ)
    N_cross = SMatrix{N,N}(K)

    return CouplingJacobians{N}(M_self, N_cross)
end

# =============================================================================
# 随伴方程式の右辺 (Pure)
# =============================================================================

"""
    adjoint_rhs_intrinsic(Q::AdjointState{N}, J::OscillatorJacobians, i::Int)
        -> (dQ_u, dQ_v)

固有ダイナミクス項: -J_i^T Q_i

J^T = | J11  J21 |
      | J12  J22 |

-J^T Q = -| J11*Q_u + J21*Q_v |
          | J12*Q_u + J22*Q_v |
"""
function adjoint_rhs_intrinsic(Q_u_i::Float64, Q_v_i::Float64, J::OscillatorJacobians)
    dQ_u = -(J.J11 * Q_u_i + J.J21 * Q_v_i)
    dQ_v = -(J.J12 * Q_u_i + J.J22 * Q_v_i)
    return (dQ_u, dQ_v)
end

"""
    adjoint_rhs_coupling_self(Q_v_i::Float64, M_self_i::Float64) -> (dQ_u, dQ_v)

結合自己項: -M_ii^T Q_i

M_ij は v 成分のみ非ゼロなので:
-M^T Q の v 成分 = -M_self_i * Q_v_i
(u 成分への寄与は 0)
"""
adjoint_rhs_coupling_self(Q_v_i::Float64, M_self_i::Float64) =
    (0.0, -M_self_i * Q_v_i)

"""
    adjoint_rhs_coupling_cross(Q_v, N_cross_row) -> (dQ_u, dQ_v)

結合交差項: -Σ_j N_ij Q_j

N_ij は v 成分のみ非ゼロなので:
-Σ_j N_ij Q_v_j = -Σ_j K_ij Q_v_j
"""
function adjoint_rhs_coupling_cross(Q_v::AbstractVector{Float64},
    N_cross_row::AbstractVector{Float64})
    dQ_v = -dot(N_cross_row, Q_v)
    return (0.0, dQ_v)
end

"""
    compute_adjoint_rhs(Q::AdjointState{N}, Js::Vector{OscillatorJacobians},
                        coupling::CouplingJacobians{N}, omega::Float64)
        -> AdjointState{N}

随伴方程式 (A1) の完全な右辺を計算

dQ/dθ = (1/ω) * [-J_i^T Q_i - Σ_j M_ij^T Q_i - Σ_j N_ij Q_j]
"""
function compute_adjoint_rhs(Q::AdjointState{N},
    Js::Vector{OscillatorJacobians},
    coupling::CouplingJacobians{N},
    omega::Float64) where N
    dQ_u = zeros(MVector{N,Float64})
    dQ_v = zeros(MVector{N,Float64})

    for i in 1:N
        # Term 1: 固有ダイナミクス -J_i^T Q_i
        dQ_u_1, dQ_v_1 = adjoint_rhs_intrinsic(Q.Q_u[i], Q.Q_v[i], Js[i])

        # Term 2: 結合自己項 -Σ_j M_ij^T Q_i
        dQ_u_2, dQ_v_2 = adjoint_rhs_coupling_self(Q.Q_v[i], coupling.M_self[i])

        # Term 3: 結合交差項 -Σ_j N_ij Q_j
        N_row = SVector{N}([coupling.N_cross[i, j] for j in 1:N])
        dQ_u_3, dQ_v_3 = adjoint_rhs_coupling_cross(Q.Q_v, N_row)

        # 合計 (ωで割る)
        dQ_u[i] = (dQ_u_1 + dQ_u_2 + dQ_u_3) / omega
        dQ_v[i] = (dQ_v_1 + dQ_v_2 + dQ_v_3) / omega
    end

    return AdjointState{N}(SVector{N}(dQ_u), SVector{N}(dQ_v))
end

# =============================================================================
# 正規化条件 (A4) (Pure)
# =============================================================================

"""
    compute_normalization_factor(Q::AdjointState{N}, F_u::Vector{Float64},
                                  F_v::Vector{Float64}) -> Float64

正規化因子を計算: Σ_i (Q_u_i * F_u_i + Q_v_i * F_v_i)

式(A4): Σ_i Q_i · (dx_i/dθ) = 1
"""
compute_normalization_factor(Q::AdjointState{N}, F_u::AbstractVector,
    F_v::AbstractVector) where N =
    dot(Q.Q_u, F_u) + dot(Q.Q_v, F_v)

"""
    normalize_adjoint_state(Q::AdjointState{N}, F_u::Vector{Float64},
                            F_v::Vector{Float64}, omega::Float64) -> AdjointState{N}

随伴変数を正規化: Q_normalized = Q * ω / normalization_factor

式(A4)を満たすように正規化
"""
function normalize_adjoint_state(Q::AdjointState{N}, F_u::AbstractVector,
    F_v::AbstractVector, omega::Float64) where N
    norm_factor = compute_normalization_factor(Q, F_u, F_v)

    if abs(norm_factor) < 1e-10
        return Q  # 数値的に不安定な場合はそのまま返す
    end

    scale = omega / norm_factor
    return AdjointState{N}(Q.Q_u * scale, Q.Q_v * scale)
end

# =============================================================================
# 1ステップの後退Euler積分 (Pure)
# =============================================================================

"""
    backward_euler_step(Q::AdjointState{N}, X::Vector{Float64},
                        fhn::FHNParams, K::Matrix{Float64},
                        dt::Float64, omega::Float64) -> AdjointState{N}

後退Euler法による1ステップ積分

Q(θ - Δθ) = Q(θ) + dt * dQ/dθ
"""
function backward_euler_step(Q::AdjointState{N}, X::Vector{Float64},
    fhn::FHNParams, K::Matrix{Float64},
    dt::Float64, omega::Float64) where N
    # 各振動子のヤコビアンを計算
    Js = compute_all_intrinsic_jacobians(X, fhn, N)

    # 結合ヤコビアンを計算
    coupling = compute_coupling_jacobians(K)

    # 右辺を計算
    dQ = compute_adjoint_rhs(Q, Js, coupling, omega)

    # Euler更新
    new_Q_u = Q.Q_u + dt * dQ.Q_u
    new_Q_v = Q.Q_v + dt * dQ.Q_v

    return AdjointState{N}(new_Q_u, new_Q_v)
end

# =============================================================================
# ベクトル場の計算 (Pure)
# =============================================================================

"""
    compute_vector_field_components(X::Vector{Float64}, params::NetworkParams)
        -> (F_u::Vector{Float64}, F_v::Vector{Float64})

ベクトル場 F(X) の u, v 成分を分離して計算
"""
function compute_vector_field_components(X::Vector{Float64}, params::NetworkParams)
    N = params.N
    F = network_full_dynamics(X, params)

    F_u = [F[2i-1] for i in 1:N]
    F_v = [F[2i] for i in 1:N]

    return (F_u, F_v)
end

# =============================================================================
# 完全な随伴方程式ソルバー (Pure - フォールドとして実装)
# =============================================================================

"""
    solve_adjoint_equation(sps::StablePeriodicSolution, params::NetworkParams,
                           n_iterations::Int) -> Matrix{Float64}

随伴方程式を解いて位相感受関数を計算

foldl を使った関数型実装：初期状態から反復的に解を改善

# Returns
- `Q::Matrix{Float64}`: 位相感受関数 (Ds*N × Ntheta)
"""
function solve_adjoint_equation(sps::StablePeriodicSolution, params::NetworkParams,
    n_iterations::Int)
    N = sps.N
    Ntheta = sps.Ntheta
    T = sps.T
    omega = sps.omega
    Xs = sps.Xs
    K = params.K
    fhn = params.fhn

    dt = T / (Ntheta - 1)

    # 初期状態を作成（全振動子に対して）
    initial_Q = AdjointState{N}(
        SVector{N}(ones(N)),
        SVector{N}(ones(N))
    )

    # ベクトル場を全時点で事前計算（純粋なマップ操作）
    vector_fields = [compute_vector_field_components(Xs[:, k], params) for k in 1:Ntheta]

    # 1回の反復を定義（純粋関数）
    function single_iteration(Q_init)
        # 後退積分（reduceを使用）
        function backward_step(Q, k::Int)
            X_k = Xs[:, k]
            F_u, F_v = vector_fields[k]

            # Euler ステップ
            Q_new = backward_euler_step(Q, X_k, fhn, K, dt, omega)

            # 正規化
            return normalize_adjoint_state(Q_new, F_u, F_v, omega)
        end

        # Ntheta から 2 まで後退（foldr相当）
        return foldr(backward_step, Ntheta:-1:2; init=Q_init)
    end

    # n_iterations 回反復（foldl を使用）
    final_Q = foldl((Q, _) -> single_iteration(Q), 1:n_iterations; init=initial_Q)

    # 最終的な後退積分で全時点の Q を記録
    Q_history = Vector{AdjointState{N}}(undef, Ntheta)
    Q_history[end] = final_Q

    for k in Ntheta:-1:2
        X_k = Xs[:, k]
        F_u, F_v = vector_fields[k]

        Q_new = backward_euler_step(Q_history[k], X_k, fhn, K, dt, omega)
        Q_normalized = normalize_adjoint_state(Q_new, F_u, F_v, omega)
        Q_history[k-1] = Q_normalized
    end

    # Matrix形式に変換
    Q_matrix = zeros(Ds * N, Ntheta)
    for k in 1:Ntheta
        for i in 1:N
            Q_matrix[2i-1, k] = Q_history[k].Q_u[i]
            Q_matrix[2i, k] = Q_history[k].Q_v[i]
        end
    end

    return Q_matrix
end

# =============================================================================
# 高階関数によるパイプライン (Pure)
# =============================================================================

"""
    create_adjoint_solver(params::NetworkParams) -> Function

パラメータを部分適用した随伴方程式ソルバーを生成（カリー化）

# Returns
- `(sps, n_iter) -> Matrix{Float64}`: 部分適用されたソルバー関数
"""
create_adjoint_solver(params::NetworkParams) =
    (sps, n_iterations) -> solve_adjoint_equation(sps, params, n_iterations)

"""
    compose_psf_pipeline(params::NetworkParams, n_iterations::Int) -> Function

PSF計算パイプラインを合成

# Returns
- `sps -> PhaseSensitivityFunction`: SPSからPSFへの変換関数
"""
function compose_psf_pipeline(params::NetworkParams, n_iterations::Int)
    solver = create_adjoint_solver(params)

    return function (sps::StablePeriodicSolution)
        Q = solver(sps, n_iterations)
        return PhaseSensitivityFunction(sps.N, sps.Ntheta, Q, sps.T, sps.omega)
    end
end

# =============================================================================
# エクスポート用ラッパー
# =============================================================================

"""
    compute_psf_separated(sps::StablePeriodicSolution, params::NetworkParams,
                          n_iterations::Int=10) -> PhaseSensitivityFunction

分離ヤコビアン法によるPSF計算のメインエントリーポイント

式(A1)の3項構造を維持した実装
"""
function compute_psf_separated(sps::StablePeriodicSolution, params::NetworkParams,
    n_iterations::Int=10)
    Q = solve_adjoint_equation(sps, params, n_iterations)
    return PhaseSensitivityFunction(sps.N, sps.Ntheta, Q, sps.T, sps.omega)
end

# =============================================================================
# 拡散結合専用の簡略化PSF計算 (PSF.m スタイル)
# =============================================================================
#
# 拡散結合 g_ij = K_ij(v_j - v_i) の場合、M_ij と N_ij が結合行列 L に統合され、
# ヤコビアンが以下の形に簡略化される:
#
#   J = | J_uu + L    J_uv  |
#       | J_vu        J_vv  |
#
# ここで L はラプラシアン行列（結合による寄与）
#
# 随伴方程式: dQ/dθ = -(1/ω) J^T Q
# 後退Euler: Q(θ - Δθ) = Q(θ) + dt * J^T * Q(θ)
#
# =============================================================================

"""
    DiffusiveAdjointWorkspace{N}

拡散結合PSF計算のワークスペース（メモリ効率化）

FDD: 計算に必要な中間データを不変構造として保持
"""
struct DiffusiveAdjointWorkspace{N}
    # 結合行列から計算されるラプラシアン
    L_diag::Vector{Float64}      # -Σ_j K_ij (対角成分)
    L_off::Matrix{Float64}       # K_ij (非対角成分)

    # 事前計算されたベクトル場
    F_u::Matrix{Float64}         # F_u[:, k] = u成分のベクトル場 at θ_k
    F_v::Matrix{Float64}         # F_v[:, k] = v成分のベクトル場 at θ_k
end

"""
    create_diffusive_workspace(sps::StablePeriodicSolution, params::NetworkParams)
        -> DiffusiveAdjointWorkspace

拡散結合用ワークスペースを作成（Pure関数）
"""
function create_diffusive_workspace(sps::StablePeriodicSolution, params::NetworkParams)
    N = sps.N
    Ntheta = sps.Ntheta
    K = params.K

    # ラプラシアン行列の構築（結合行列から）
    # L_diag[i] = -Σ_j K_ij（行和の負）
    L_diag = [-sum(K[i, :]) for i in 1:N]
    L_off = copy(K)

    # ベクトル場の事前計算
    F_u = zeros(N, Ntheta)
    F_v = zeros(N, Ntheta)

    for k in 1:Ntheta
        X_k = sps.Xs[:, k]
        F_k = network_full_dynamics(X_k, params)
        for i in 1:N
            F_u[i, k] = F_k[2i-1]
            F_v[i, k] = F_k[2i]
        end
    end

    return DiffusiveAdjointWorkspace{N}(L_diag, L_off, F_u, F_v)
end

"""
    compute_simplified_jacobian_transpose_action(
        Q_u::Vector{Float64}, Q_v::Vector{Float64},
        u::Vector{Float64}, v::Vector{Float64},
        fhn::FHNParams, ws::DiffusiveAdjointWorkspace{N}
    ) -> (dQ_u, dQ_v)

PSF.mスタイルの簡略化ヤコビアン転置の作用を計算（Pure関数）

FHN + 拡散結合のヤコビアン:
    J = | -δb I        δ I         |
        | -I + L      diag(1-v²)   |

ここで L はラプラシアン（拡散結合の寄与）

J^T の作用:
    J^T Q = | -δb Q_u - Q_v + L^T Q_v |
            | δ Q_u + (1-v²) Q_v      |

注: v成分への結合の場合、L は J_vv ブロックに入る
"""
function compute_simplified_jacobian_transpose_action(
    Q_u::Vector{Float64}, Q_v::Vector{Float64},
    u::Vector{Float64}, v::Vector{Float64},
    fhn::FHNParams, ws::DiffusiveAdjointWorkspace{N}) where N

    delta = fhn.delta
    b = fhn.b

    # J^T Q の計算（PSF.mスタイル）
    # dQ_u = -δb * Q_u - Q_v + L^T * Q_v  (L^T Q_v は v からの結合)
    # dQ_v = δ * Q_u + (1-v²) * Q_v

    dQ_u = zeros(N)
    dQ_v = zeros(N)

    for i in 1:N
        # 固有ダイナミクス項
        dQ_u[i] = -delta * b * Q_u[i] - Q_v[i]
        dQ_v[i] = delta * Q_u[i] + (1.0 - v[i]^2) * Q_v[i]

        # 結合項（v成分への拡散結合の場合、J_vv に L が入る）
        # L^T Q_v の計算: 対角項 + 非対角項
        coupling_term = ws.L_diag[i] * Q_v[i]
        for j in 1:N
            if i != j
                coupling_term += ws.L_off[j, i] * Q_v[j]  # L^T なので [j,i]
            end
        end
        dQ_v[i] += coupling_term
    end

    return (dQ_u, dQ_v)
end

"""
    normalize_psf_diffusive!(
        Q_u::Vector{Float64}, Q_v::Vector{Float64},
        F_u::AbstractVector, F_v::AbstractVector,
        omega::Float64
    )

PSFの正規化条件を適用（in-place）

正規化条件: Q · F = ω
"""
function normalize_psf_diffusive!(
    Q_u::Vector{Float64}, Q_v::Vector{Float64},
    F_u::AbstractVector, F_v::AbstractVector,
    omega::Float64)

    norm_factor = dot(Q_u, F_u) + dot(Q_v, F_v)

    if abs(norm_factor) > 1e-10
        scale = omega / norm_factor
        Q_u .*= scale
        Q_v .*= scale
    end
end

"""
    solve_adjoint_diffusive(
        sps::StablePeriodicSolution,
        params::NetworkParams,
        ::DiffusiveCoupling,
        n_iterations::Int
    ) -> Matrix{Float64}

拡散結合専用の随伴方程式ソルバー（PSF.mスタイル）

多重ディスパッチにより DiffusiveCoupling 型で呼び出された場合に使用

# Algorithm (PSF.m に準拠)
1. ワークスペースを作成（ベクトル場・ラプラシアンの事前計算）
2. n_iterations 回の収束反復
3. 各反復で後退Euler積分 + 正規化
"""
function solve_adjoint_diffusive(
    sps::StablePeriodicSolution,
    params::NetworkParams,
    ::DiffusiveCoupling,
    n_iterations::Int)

    N = sps.N
    Ntheta = sps.Ntheta
    T = sps.T
    omega = sps.omega
    Xs = sps.Xs
    fhn = params.fhn

    dt = T / (Ntheta - 1)

    # ワークスペース作成（純粋関数）
    ws = create_diffusive_workspace(sps, params)

    # 初期条件（任意の非ゼロベクトル）
    Q_u = ones(N)
    Q_v = ones(N)

    # 収束反復（PSF.m の rep=1:40 に対応）
    for rep in 1:n_iterations
        # 後退積分（θ: 2π → 0）
        for k in Ntheta:-1:1
            # 状態を取得
            u_k = [Xs[2i-1, k] for i in 1:N]
            v_k = [Xs[2i, k] for i in 1:N]

            # 簡略化ヤコビアン転置の作用
            dQ_u, dQ_v = compute_simplified_jacobian_transpose_action(
                Q_u, Q_v, u_k, v_k, fhn, ws)

            # 後退Euler更新: Q(θ-Δθ) = Q(θ) + dt * J^T * Q(θ)
            Q_u .+= dt .* dQ_u
            Q_v .+= dt .* dQ_v

            # 正規化条件の適用
            F_u_k = @view ws.F_u[:, k]
            F_v_k = @view ws.F_v[:, k]
            normalize_psf_diffusive!(Q_u, Q_v, F_u_k, F_v_k, omega)
        end
    end

    # 最終結果を記録しながら後退積分
    Q_result = zeros(Ds * N, Ntheta)

    for k in Ntheta:-1:1
        # 結果を保存
        for i in 1:N
            Q_result[2i-1, k] = Q_u[i]
            Q_result[2i, k] = Q_v[i]
        end

        if k > 1
            u_k = [Xs[2i-1, k] for i in 1:N]
            v_k = [Xs[2i, k] for i in 1:N]

            dQ_u, dQ_v = compute_simplified_jacobian_transpose_action(
                Q_u, Q_v, u_k, v_k, fhn, ws)

            Q_u .+= dt .* dQ_u
            Q_v .+= dt .* dQ_v

            F_u_k = @view ws.F_u[:, k]
            F_v_k = @view ws.F_v[:, k]
            normalize_psf_diffusive!(Q_u, Q_v, F_u_k, F_v_k, omega)
        end
    end

    # 周期境界条件: Q(2π) = Q(0)
    Q_result[:, end] = Q_result[:, 1]

    return Q_result
end

"""
    compute_psf(sps::StablePeriodicSolution, params::NetworkParams,
                coupling::CouplingType, n_iterations::Int=40)
        -> PhaseSensitivityFunction

PSF計算のメインエントリーポイント（多重ディスパッチ版）

# Arguments
- `sps`: 安定周期解
- `params`: ネットワークパラメータ
- `coupling`: 結合タイプ（DiffusiveCoupling または GeneralCoupling）
- `n_iterations`: 収束反復回数（デフォルト: 40、PSF.mと同じ）

# Returns
- `PhaseSensitivityFunction`: 位相感受関数

# Example
```julia
sps = compute_sps(params)
psf = compute_psf(sps, params, DiffusiveCoupling())  # 簡略化版
psf = compute_psf(sps, params, GeneralCoupling())    # 一般版
```
"""
function compute_psf(sps::StablePeriodicSolution, params::NetworkParams,
    coupling::DiffusiveCoupling, n_iterations::Int=40)
    Q = solve_adjoint_diffusive(sps, params, coupling, n_iterations)
    return PhaseSensitivityFunction(sps.N, sps.Ntheta, Q, sps.T, sps.omega)
end

# GeneralCoupling の場合は既存の分離実装を使用
function compute_psf(sps::StablePeriodicSolution, params::NetworkParams,
    ::GeneralCoupling, n_iterations::Int=10)
    return compute_psf_separated(sps, params, n_iterations)
end

# 注意: デフォルトのcompute_psf関数はファイル末尾で定義（version引数付き）

# =============================================================================
# 旧Juliaコード (phase_sensitivity.jl) 準拠の実装 - V2
# =============================================================================
#
# 変更履歴:
# - 2024: オリジナル実装 (solve_adjoint_diffusive)
#   - 順序: 更新 → 正規化
#   - 問題: Full結合でPSFが10^4オーダーに発散
#
# - V2実装 (solve_adjoint_diffusive_v2):
#   - 順序: 正規化 → 更新（旧Juliaコード準拠）
#   - 参照: sample-codes/myOld-ProjectDir/src/phase_sensitivity.jl
#
# =============================================================================

"""
    solve_adjoint_diffusive_v2(
        sps::StablePeriodicSolution,
        params::NetworkParams,
        n_iterations::Int
    ) -> Matrix{Float64}

拡散結合専用の随伴方程式ソルバー V2（旧Juliaコード phase_sensitivity.jl 準拠）

# 旧実装との違い
- 順序: 正規化を先に実行 → 固有Jacobian更新 → 結合項追加
- 旧実装は: 更新 → 正規化 の順序だった

# Algorithm (phase_sensitivity.jl の update_Q! に準拠)
1. 正規化を先に実行: Q *= ω / (F' * Q)
2. 固有ヤコビアン転置で更新: Q += J^T * Q * dt
3. 結合項を追加: Q += (M + N項) * dt
"""
function solve_adjoint_diffusive_v2(
    sps::StablePeriodicSolution,
    params::NetworkParams,
    n_iterations::Int)

    N = sps.N
    Ntheta = sps.Ntheta
    T = sps.T
    omega = sps.omega
    Xs = sps.Xs
    fhn = params.fhn
    K = params.K

    dt = T / (Ntheta - 1)

    # ベクトル場を事前計算
    F_all = zeros(Ds * N, Ntheta)
    for k in 1:Ntheta
        F_all[:, k] = network_full_dynamics(Xs[:, k], params)
    end

    # 初期条件（任意の非ゼロベクトル）
    Q = ones(Ds * N)

    # 作業用配列（一括更新のため）
    dQ_intrinsic = zeros(Ds * N)
    dQ_coupling = zeros(Ds * N)

    # 収束反復（旧コード: REP回）
    for rep in 1:n_iterations
        # 後退積分（θ: Nθ → 1）
        for tt in Ntheta:-1:1
            X = Xs[:, tt]
            F = F_all[:, tt]

            # Step 1: 正規化を先に実行（旧コード準拠: Q[:] *= ω/s）
            s = dot(F, Q)
            if abs(s) > 1e-10
                Q .*= omega / s
            end

            # Step 2: 固有ヤコビアン転置を計算（旧コード: Q[:] += Jᵢⱼ * Q * dt）
            # 全要素を一括で計算してから更新（旧コードの行列乗算と同等）
            dQ_intrinsic .= 0.0
            for i in 1:N
                v_i = X[2i]
                Q_u_i = Q[2i-1]
                Q_v_i = Q[2i]

                # J^T * Q の計算
                dQ_intrinsic[2i-1] = (-fhn.delta * fhn.b) * Q_u_i + (-1.0) * Q_v_i
                dQ_intrinsic[2i] = fhn.delta * Q_u_i + (1.0 - v_i^2) * Q_v_i
            end
            Q .+= dQ_intrinsic .* dt

            # Step 3: 結合項を計算（旧コード準拠: 一括計算してから更新）
            # tempM[i] = Σ_j M_ij * Q_i  (M_ij[2,2] = -K[i,j] for i≠j)
            # tempN[i] = Σ_j N_ji * Q_j  (N_ji[2,2] = +K[j,i] for j≠i)
            dQ_coupling .= 0.0
            for i in 1:N
                tempM_v = 0.0
                tempN_v = 0.0

                for j in 1:N
                    if i != j
                        # M_ij の (2,2) 成分 = -K[i,j]
                        tempM_v += (-K[i, j]) * Q[2i]

                        # N_ji の (2,2) 成分 = K[j,i]
                        tempN_v += K[j, i] * Q[2j]
                    end
                end

                # v成分のみ
                dQ_coupling[2i] = tempM_v + tempN_v
            end
            Q .+= dQ_coupling .* dt
        end
    end

    # 最終結果を記録しながら後退積分
    Q_result = zeros(Ds * N, Ntheta)
    Q .= 1.0  # 再初期化

    for rep in 1:(n_iterations + 1)
        for tt in Ntheta:-1:1
            X = Xs[:, tt]
            F = F_all[:, tt]

            # Step 1: 正規化
            s = dot(F, Q)
            if abs(s) > 1e-10
                Q .*= omega / s
            end

            # 最終反復で記録
            if rep == n_iterations + 1
                Q_result[:, tt] = copy(Q)
            end

            # Step 2: 固有ヤコビアン転置を計算（一括更新）
            dQ_intrinsic .= 0.0
            for i in 1:N
                v_i = X[2i]
                Q_u_i = Q[2i-1]
                Q_v_i = Q[2i]

                dQ_intrinsic[2i-1] = (-fhn.delta * fhn.b) * Q_u_i + (-1.0) * Q_v_i
                dQ_intrinsic[2i] = fhn.delta * Q_u_i + (1.0 - v_i^2) * Q_v_i
            end
            Q .+= dQ_intrinsic .* dt

            # Step 3: 結合項を計算（一括更新）
            dQ_coupling .= 0.0
            for i in 1:N
                tempM_v = 0.0
                tempN_v = 0.0

                for j in 1:N
                    if i != j
                        tempM_v += (-K[i, j]) * Q[2i]
                        tempN_v += K[j, i] * Q[2j]
                    end
                end

                dQ_coupling[2i] = tempM_v + tempN_v
            end
            Q .+= dQ_coupling .* dt
        end
    end

    return Q_result
end

"""
    compute_psf_v2(sps, params; n_iterations=40) -> PhaseSensitivityFunction

PSF計算 V2版（旧Juliaコード準拠）

順序が異なる:
- V1 (旧): 更新 → 正規化
- V2 (新): 正規化 → 更新
"""
function compute_psf_v2(sps::StablePeriodicSolution, params::NetworkParams;
    n_iterations::Int=40)
    Q = solve_adjoint_diffusive_v2(sps, params, n_iterations)
    return PhaseSensitivityFunction(sps.N, sps.Ntheta, Q, sps.T, sps.omega)
end

# =============================================================================
# V3: Cコード準拠の細かい時間刻み版
# =============================================================================
#
# V2の問題: dt = T/(Ntheta-1) ≈ 0.757 は大きすぎて不安定
# Cコード: DT = 0.001 で約75,000ステップ
#
# V3解決策: 細かい内部時間刻みを使用し、補間で軌道を取得
#
# =============================================================================

"""
    solve_adjoint_diffusive_v3(
        sps::StablePeriodicSolution,
        params::NetworkParams,
        n_iterations::Int;
        dt_internal::Float64=0.001
    ) -> Matrix{Float64}

拡散結合専用の随伴方程式ソルバー V3（Cコード networkPSF.c 準拠の細かい時間刻み）

# V2との違い
- 細かい内部時間刻み dt_internal = 0.001 を使用（Cコード準拠）
- 周期軌道を線形補間で細かくリサンプル
- 数値的に安定

# Algorithm
1. 周期軌道を dt_internal 間隔でリサンプル（線形補間）
2. n_iterations 回の収束反復
3. 正規化 → 固有Jacobian更新 → 結合項追加（旧Juliaコード順序）
"""
function solve_adjoint_diffusive_v3(
    sps::StablePeriodicSolution,
    params::NetworkParams,
    n_iterations::Int;
    dt_internal::Float64=0.001)

    N = sps.N
    Ntheta = sps.Ntheta
    T = sps.T
    omega = sps.omega
    Xs = sps.Xs
    fhn = params.fhn
    K = params.K

    # 細かい時間刻み（Cコード準拠: DT = 0.001）
    n_fine = Int(ceil(T / dt_internal))
    dt = T / n_fine

    # 周期軌道を細かい時間刻みでリサンプル（線形補間）
    Xs_fine = zeros(Ds * N, n_fine + 1)
    F_fine = zeros(Ds * N, n_fine + 1)

    for tt in 0:n_fine
        t = tt * dt
        # 元のグリッドでの位置を計算
        theta_frac = t / T * (Ntheta - 1)
        k_low = Int(floor(theta_frac)) + 1
        k_high = min(k_low + 1, Ntheta)
        alpha = theta_frac - (k_low - 1)

        # 線形補間
        if k_low >= Ntheta
            Xs_fine[:, tt+1] = Xs[:, Ntheta]
        else
            Xs_fine[:, tt+1] = (1 - alpha) * Xs[:, k_low] + alpha * Xs[:, k_high]
        end

        # ベクトル場を計算
        F_fine[:, tt+1] = network_full_dynamics(Xs_fine[:, tt+1], params)
    end

    # 初期条件
    Q = ones(Ds * N)

    # 全増分を格納する配列
    dQ_total = zeros(Ds * N)

    # 結合項用の配列
    dQ_coupling = zeros(Ds * N)

    # 収束反復
    for rep in 1:n_iterations
        # 後退積分（tt: n_fine → 0）
        for tt in n_fine:-1:0
            X = Xs_fine[:, tt+1]
            F = F_fine[:, tt+1]

            # Step 1: 正規化（旧コード準拠: 条件なしで直接正規化）
            s = dot(F, Q)
            Q .*= omega / s

            # Step 2: 固有ヤコビアン更新（旧コード: Q[:] += Jᵢⱼ * Q * dt）
            dQ_total .= 0.0
            for i in 1:N
                v_i = X[2i]
                Q_u_i = Q[2i-1]
                Q_v_i = Q[2i]

                dQ_total[2i-1] = (-fhn.delta * fhn.b) * Q_u_i + (-1.0) * Q_v_i
                dQ_total[2i] = fhn.delta * Q_u_i + (1.0 - v_i^2) * Q_v_i
            end
            Q .+= dQ_total .* dt

            # Step 3: 結合項更新（旧コード準拠: 更新後のQを使用）
            dQ_coupling .= 0.0
            for i in 1:N
                tempM_v = 0.0
                tempN_v = 0.0

                for j in 1:N
                    if i != j
                        tempM_v += (-K[i, j]) * Q[2i]
                        tempN_v += K[j, i] * Q[2j]
                    end
                end

                dQ_coupling[2i] = tempM_v + tempN_v
            end
            Q .+= dQ_coupling .* dt
        end
    end

    # 最終結果を記録（Ntheta点にダウンサンプル）
    Q_result = zeros(Ds * N, Ntheta)
    Q .= 1.0

    for rep in 1:(n_iterations + 1)
        for tt in n_fine:-1:0
            X = Xs_fine[:, tt+1]
            F = F_fine[:, tt+1]

            # Step 1: 正規化（旧コード準拠: 条件なしで直接正規化）
            s = dot(F, Q)
            Q .*= omega / s

            # 最終反復で記録（等間隔位相でサンプル）
            if rep == n_iterations + 1
                for k in 1:Ntheta
                    target_tt = Int(round((k-1) / (Ntheta-1) * n_fine))
                    if target_tt == tt
                        Q_result[:, k] = copy(Q)
                    end
                end
            end

            # Step 2: 固有ヤコビアン更新
            dQ_total .= 0.0
            for i in 1:N
                v_i = X[2i]
                Q_u_i = Q[2i-1]
                Q_v_i = Q[2i]

                dQ_total[2i-1] = (-fhn.delta * fhn.b) * Q_u_i + (-1.0) * Q_v_i
                dQ_total[2i] = fhn.delta * Q_u_i + (1.0 - v_i^2) * Q_v_i
            end
            Q .+= dQ_total .* dt

            # Step 3: 結合項更新（旧コード準拠: 更新後のQを使用）
            dQ_coupling .= 0.0
            for i in 1:N
                tempM_v = 0.0
                tempN_v = 0.0

                for j in 1:N
                    if i != j
                        tempM_v += (-K[i, j]) * Q[2i]
                        tempN_v += K[j, i] * Q[2j]
                    end
                end

                dQ_coupling[2i] = tempM_v + tempN_v
            end
            Q .+= dQ_coupling .* dt
        end
    end

    return Q_result
end

"""
    compute_psf_v3(sps, params; n_iterations=40, dt_internal=0.001)

PSF計算 V3版（Cコード準拠の細かい時間刻み）
"""
function compute_psf_v3(sps::StablePeriodicSolution, params::NetworkParams;
    n_iterations::Int=40, dt_internal::Float64=0.001)
    Q = solve_adjoint_diffusive_v3(sps, params, n_iterations; dt_internal=dt_internal)
    return PhaseSensitivityFunction(sps.N, sps.Ntheta, Q, sps.T, sps.omega)
end

# =============================================================================
# V4: PSF.m 完全準拠版（Nakao-san paper スタイル）
# =============================================================================
#
# 変数対応関係（重要！）:
#   MATLAB PSF.m (Nakao)  |  Our Project (NOLTA2025)  |  役割
#   ----------------------+---------------------------+-------------
#   u                     |  v                        |  膜電位（速い変数）
#   v                     |  u                        |  回復変数（遅い変数）
#   Zu                    |  Zv                       |  膜電位のPSF
#   Zv                    |  Zu                       |  回復変数のPSF
#
# ダイナミクス対応:
#   MATLAB: du/dt = u - u³/3 - v + L·u + I,  dv/dt = δ(u + a - bv)
#   Ours:   dv/dt = v - v³/3 - u + I + L·v,  du/dt = δ(a + v - bu)
#
# ヤコビアン対応:
#   MATLAB J:               Our J:
#   | 1-u² + L   -1   |     | -δb    δ        |
#   | δ          -δb  |     | -1     1-v² + L |
#
# =============================================================================

"""
    solve_adjoint_psf_m_style(
        sps::StablePeriodicSolution,
        params::NetworkParams,
        n_iterations::Int;
        dt_internal::Float64=0.001
    ) -> Matrix{Float64}

PSF.m のアルゴリズムを完全に準拠した随伴方程式ソルバー

# Algorithm (PSF.m に完全準拠)

PSF.m の処理順序:
```matlab
for t = Tp:-1:1
    % 1. ヤコビアン転置を計算
    J_T = [J11, J12; J21, J22];
    z = [z_u; z_v];
    z_dot = J_T' * z;

    % 2. Euler更新（先に更新）
    z_u = z_u + t_step * z_dot(1:nodes);
    z_v = z_v + t_step * z_dot(nodes+1:end);

    % 3. 正規化（更新後に正規化）
    norm_factor = sum(z_u.*Fu(:,t) + z_v.*Fv(:,t));
    z_u = z_u / norm_factor * omega;
    z_v = z_v / norm_factor * omega;
end
```

# 変数対応（MATLAB→Julia）
- MATLAB u → Julia v（膜電位）
- MATLAB v → Julia u（回復変数）
- MATLAB J11 = diag(1-u²)+L → Julia: v成分への寄与
- MATLAB J22 = -δb → Julia: u成分への寄与

# Returns
- `Q::Matrix{Float64}`: 位相感受関数 (Ds*N × Ntheta)
"""
function solve_adjoint_psf_m_style(
    sps::StablePeriodicSolution,
    params::NetworkParams,
    n_iterations::Int;
    dt_internal::Float64=0.001)

    N = sps.N
    Ntheta = sps.Ntheta
    T = sps.T
    omega = sps.omega
    Xs = sps.Xs
    fhn = params.fhn
    K = params.K

    # 細かい時間刻み（PSF.mでは t_step に対応）
    n_fine = Int(ceil(T / dt_internal))
    dt = T / n_fine

    # -------------------------------------------------------------------
    # Step 1: 周期軌道を細かい時間刻みでリサンプル
    # -------------------------------------------------------------------
    Xs_fine = zeros(Ds * N, n_fine + 1)

    for tt in 0:n_fine
        t = tt * dt
        theta_frac = t / T * (Ntheta - 1)
        k_low = Int(floor(theta_frac)) + 1
        k_high = min(k_low + 1, Ntheta)
        alpha = theta_frac - (k_low - 1)

        if k_low >= Ntheta
            Xs_fine[:, tt+1] = Xs[:, Ntheta]
        else
            Xs_fine[:, tt+1] = (1 - alpha) * Xs[:, k_low] + alpha * Xs[:, k_high]
        end
    end

    # -------------------------------------------------------------------
    # Step 2: ベクトル場 F を全時点で計算（PSF.m の Fu, Fv に対応）
    # -------------------------------------------------------------------
    F_u = zeros(N, n_fine + 1)  # u成分（回復変数）
    F_v = zeros(N, n_fine + 1)  # v成分（膜電位）

    for tt in 0:n_fine
        X = Xs_fine[:, tt+1]
        F = network_full_dynamics(X, params)
        for i in 1:N
            F_u[i, tt+1] = F[2i-1]  # du/dt
            F_v[i, tt+1] = F[2i]    # dv/dt
        end
    end

    # -------------------------------------------------------------------
    # Step 3: ラプラシアン行列 L を構築（結合項）
    #
    # 拡散結合: coupling_i = Σ_j K_ij (v_j - v_i)
    # これは dv_i/dt に加算されるので、ヤコビアンの v-v ブロックに入る
    #
    # L_diag[i] = -Σ_j K_ij（対角成分）
    # L_off[i,j] = K_ij（非対角成分）
    # -------------------------------------------------------------------
    L_diag = [-sum(K[i, :]) for i in 1:N]
    L_off = K

    # -------------------------------------------------------------------
    # Step 4: 初期条件（PSF.m: z_u = ones(nodes,1), z_v = ones(nodes,1)）
    # -------------------------------------------------------------------
    z_u = ones(N)  # u成分（回復変数）のPSF
    z_v = ones(N)  # v成分（膜電位）のPSF

    # -------------------------------------------------------------------
    # Step 5: 収束反復（PSF.m: for rep=1:40）
    # -------------------------------------------------------------------
    for _ in 1:n_iterations
        # 後退積分（PSF.m: for t = Tp:-1:1）
        for tt in n_fine:-1:0
            X = Xs_fine[:, tt+1]

            # 状態を取得（我々の規約: u=回復, v=膜電位）
            # u_t は回復変数だが、ヤコビアンでは v_t のみ使用
            v_t = [X[2i] for i in 1:N]    # 膜電位

            # -------------------------------------------------------
            # ヤコビアン転置 J^T の作用を計算
            #
            # 我々のヤコビアン J:
            #   J = | ∂(du)/∂u  ∂(du)/∂v |   | -δb    δ        |
            #       | ∂(dv)/∂u  ∂(dv)/∂v | = | -1     1-v² + L |
            #
            # J^T = | -δb    -1        |
            #       | δ      1-v² + L  |
            #
            # z_dot = J^T * z:
            #   dz_u = -δb * z_u - z_v
            #   dz_v = δ * z_u + (1-v²) * z_v + L^T * z_v
            # -------------------------------------------------------
            dz_u = zeros(N)
            dz_v = zeros(N)

            for i in 1:N
                # 固有ダイナミクス項（J^T の作用）
                dz_u[i] = -fhn.delta * fhn.b * z_u[i] - z_v[i]
                dz_v[i] = fhn.delta * z_u[i] + (1.0 - v_t[i]^2) * z_v[i]

                # 結合項（L^T * z_v）
                # L^T[i,j] = L[j,i] なので、行和ではなく列和に対応
                coupling_term = L_diag[i] * z_v[i]  # 対角項
                for j in 1:N
                    if i != j
                        coupling_term += L_off[j, i] * z_v[j]  # L^T なので [j,i]
                    end
                end
                dz_v[i] += coupling_term
            end

            # -------------------------------------------------------
            # Euler更新（PSF.m: z_u = z_u + t_step*z_dot）
            # -------------------------------------------------------
            z_u .+= dt .* dz_u
            z_v .+= dt .* dz_v

            # -------------------------------------------------------
            # 正規化（PSF.m: norm_factor = sum(z_u.*Fu + z_v.*Fv)）
            # -------------------------------------------------------
            norm_factor = dot(z_u, @view(F_u[:, tt+1])) + dot(z_v, @view(F_v[:, tt+1]))
            z_u .*= omega / norm_factor
            z_v .*= omega / norm_factor
        end
    end

    # -------------------------------------------------------------------
    # Step 6: 最終結果を記録（Ntheta点にダウンサンプル）
    # -------------------------------------------------------------------
    Q_result = zeros(Ds * N, Ntheta)

    # 再度後退積分して全時点を記録
    z_u .= 1.0
    z_v .= 1.0

    for rep in 1:(n_iterations + 1)
        for tt in n_fine:-1:0
            X = Xs_fine[:, tt+1]
            v_t = [X[2i] for i in 1:N]  # 膜電位のみ使用

            # J^T の作用
            dz_u = zeros(N)
            dz_v = zeros(N)

            for i in 1:N
                dz_u[i] = -fhn.delta * fhn.b * z_u[i] - z_v[i]
                dz_v[i] = fhn.delta * z_u[i] + (1.0 - v_t[i]^2) * z_v[i]

                coupling_term = L_diag[i] * z_v[i]
                for j in 1:N
                    if i != j
                        coupling_term += L_off[j, i] * z_v[j]
                    end
                end
                dz_v[i] += coupling_term
            end

            # 更新
            z_u .+= dt .* dz_u
            z_v .+= dt .* dz_v

            # 正規化
            norm_factor = dot(z_u, @view(F_u[:, tt+1])) + dot(z_v, @view(F_v[:, tt+1]))
            z_u .*= omega / norm_factor
            z_v .*= omega / norm_factor

            # 最終反復で記録
            if rep == n_iterations + 1
                for k in 1:Ntheta
                    target_tt_k = Int(round((k-1) / (Ntheta-1) * n_fine))
                    if target_tt_k == tt
                        for i in 1:N
                            Q_result[2i-1, k] = z_u[i]
                            Q_result[2i, k] = z_v[i]
                        end
                    end
                end
            end
        end
    end

    return Q_result
end

"""
    compute_psf_v4(sps, params; n_iterations=40, dt_internal=0.001)

PSF計算 V4版（PSF.m 完全準拠）

# PSF.m との対応

アルゴリズムの順序（PSF.m に準拠）:
1. ヤコビアン転置を計算
2. Euler更新（先に更新）
3. 正規化（更新後に正規化）

変数対応:
- MATLAB u（膜電位）→ Julia v
- MATLAB v（回復変数）→ Julia u
- MATLAB Zu → Julia Zv
- MATLAB Zv → Julia Zu
"""
function compute_psf_v4(sps::StablePeriodicSolution, params::NetworkParams;
    n_iterations::Int=40, dt_internal::Float64=0.001)
    Q = solve_adjoint_psf_m_style(sps, params, n_iterations; dt_internal=dt_internal)
    return PhaseSensitivityFunction(sps.N, sps.Ntheta, Q, sps.T, sps.omega)
end

# =============================================================================
# 汎用インターフェース: compute_psf (多重ディスパッチ版)
# =============================================================================

"""
    compute_psf(sps, params, ::Val{V}; kwargs...) -> PhaseSensitivityFunction

PSF計算の統一インターフェース。Valによるバージョン指定で実装を切り替え。

# Arguments
- `sps::StablePeriodicSolution`: 安定周期解
- `params::NetworkParams`: ネットワークパラメータ
- `::Val{V}`: バージョン指定 (Val(2), Val(3), または Val(4))

# Keyword Arguments
- `n_iterations::Int`: 反復回数（デフォルト: 40）
- `dt_internal::Float64`: 内部時間刻み（V3, V4のみ、デフォルト: 0.001）

# Versions
- V2: 旧Juliaコード準拠（正規化→更新 の順序）
- V3: Cコード準拠（細かい時間刻み）
- V4: PSF.m 完全準拠（更新→正規化 の順序、推奨）

# Examples
```julia
# V4版（PSF.m 完全準拠、推奨）
psf = compute_psf(sps, params, Val(4); n_iterations=40)

# V2版（旧Juliaコード準拠の順序）
psf = compute_psf(sps, params, Val(2); n_iterations=40)

# V3版（Cコード準拠の細かい時間刻み）
psf = compute_psf(sps, params, Val(3); n_iterations=40)
```
"""
function compute_psf(sps::StablePeriodicSolution, params::NetworkParams, ::Val{2};
                     n_iterations::Int=40)
    return compute_psf_v2(sps, params; n_iterations=n_iterations)
end

function compute_psf(sps::StablePeriodicSolution, params::NetworkParams, ::Val{3};
                     n_iterations::Int=40,
                     dt_internal::Float64=0.001)
    return compute_psf_v3(sps, params; n_iterations=n_iterations, dt_internal=dt_internal)
end

function compute_psf(sps::StablePeriodicSolution, params::NetworkParams, ::Val{4};
                     n_iterations::Int=40,
                     dt_internal::Float64=0.001)
    return compute_psf_v4(sps, params; n_iterations=n_iterations, dt_internal=dt_internal)
end

# デフォルトバージョン（V4を推奨）
"""
    compute_psf(sps, params; version=4, kwargs...) -> PhaseSensitivityFunction

キーワード引数でバージョンを指定するインターフェース。

# Versions
- version=2: 旧Juliaコード準拠
- version=3: Cコード準拠の細かい時間刻み
- version=4: PSF.m 完全準拠（デフォルト、推奨）

# Examples
```julia
psf = compute_psf(sps, params; version=4, n_iterations=40)
```
"""
function compute_psf(sps::StablePeriodicSolution, params::NetworkParams;
                     version::Int=4,
                     n_iterations::Int=40,
                     dt_internal::Float64=0.001)
    if version == 2
        return compute_psf_v2(sps, params; n_iterations=n_iterations)
    elseif version == 3
        return compute_psf_v3(sps, params; n_iterations=n_iterations, dt_internal=dt_internal)
    elseif version == 4
        return compute_psf_v4(sps, params; n_iterations=n_iterations, dt_internal=dt_internal)
    else
        error("Unknown PSF version: $version. Supported versions: 2, 3, 4")
    end
end

