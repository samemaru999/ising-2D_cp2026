"""
    Logic/FHN.jl - FitzHugh-Nagumo振動子力学 [Pure]

FHN振動子の力学系を定義する純粋関数群

力学方程式:
    du/dt = δ(a + v - b*u)
    dv/dt = v - v³/3 - u + I

Reference: FitzHugh (1961), Nagumo et al. (1962)
"""

using StaticArrays
using LinearAlgebra

# =============================================================================
# 単一振動子の力学
# =============================================================================

"""
    fhn_single_dynamics(u, v, I, fhn::FHNParams) -> (du, dv)

単一FHN振動子の力学を計算（Pure関数）

# Arguments
- `u::Float64`: 回復変数
- `v::Float64`: 膜電位変数
- `I::Float64`: 外部入力
- `fhn::FHNParams`: FHNパラメータ

# Returns
- `(du, dv)::Tuple{Float64, Float64}`: 時間微分
"""
function fhn_single_dynamics(u::Float64, v::Float64, I::Float64, fhn::FHNParams)
    du = fhn.delta * (fhn.a + v - fhn.b * u)
    dv = v - v^3 / 3.0 - u + I
    return (du, dv)
end

"""
    fhn_single_dynamics(state::OscillatorState, I, fhn) -> OscillatorState

ベクトル形式のFHN力学（Pure関数）
"""
function fhn_single_dynamics(state::OscillatorState, I::Float64, fhn::FHNParams)
    u, v = state[1], state[2]
    du, dv = fhn_single_dynamics(u, v, I, fhn)
    return OscillatorState(du, dv)
end

# =============================================================================
# 単一振動子のヤコビアン
# =============================================================================

"""
    fhn_single_jacobian(u, v, fhn::FHNParams) -> SMatrix{2,2}

単一FHN振動子のヤコビアン行列（Pure関数）

J = | ∂(du)/∂u  ∂(du)/∂v |   | -δb    δ      |
    | ∂(dv)/∂u  ∂(dv)/∂v | = | -1     1 - v² |

# Returns
- `J::SMatrix{2,2,Float64}`: ヤコビアン行列 ∂f/∂x
"""
function fhn_single_jacobian(u::Float64, v::Float64, fhn::FHNParams)
    return SMatrix{2,2}(
        -fhn.delta * fhn.b,  -1.0,
        fhn.delta,           1.0 - v^2
    )
end

"""
    fhn_single_jacobian(state::OscillatorState, fhn) -> SMatrix{2,2}

状態ベクトル形式のヤコビアン計算（Pure関数）
"""
function fhn_single_jacobian(state::OscillatorState, fhn::FHNParams)
    return fhn_single_jacobian(state[1], state[2], fhn)
end

# =============================================================================
# ネットワーク力学（内部結合を含む）
# =============================================================================

"""
    network_intrinsic_dynamics(X, params::NetworkParams) -> Vector{Float64}

ネットワーク内部力学を計算（内部結合なし）（Pure関数）

# Arguments
- `X::Vector{Float64}`: ネットワーク状態 (長さ = Ds * N)
- `params::NetworkParams`: ネットワークパラメータ

# Returns
- `dX::Vector{Float64}`: 時間微分ベクトル
"""
function network_intrinsic_dynamics(X::Vector{Float64}, params::NetworkParams)
    N = params.N
    dX = zeros(Ds * N)

    for i in 1:N
        u = X[2i-1]
        v = X[2i]
        du, dv = fhn_single_dynamics(u, v, params.I[i], params.fhn)
        dX[2i-1] = du
        dX[2i] = dv
    end

    return dX
end

"""
    network_full_dynamics(X, params::NetworkParams) -> Vector{Float64}

ネットワーク内結合を含む完全な力学（Pure関数）

結合項: K_{ij} * (v_j - v_i) を v_i の方程式に加算

# Arguments
- `X::Vector{Float64}`: ネットワーク状態 (長さ = Ds * N)
- `params::NetworkParams`: ネットワークパラメータ

# Returns
- `dX::Vector{Float64}`: 時間微分ベクトル
"""
function network_full_dynamics(X::Vector{Float64}, params::NetworkParams)
    N = params.N
    K = params.K
    dX = zeros(Ds * N)

    for i in 1:N
        u_i = X[2i-1]
        v_i = X[2i]

        # 内部力学
        du, dv = fhn_single_dynamics(u_i, v_i, params.I[i], params.fhn)

        # ネットワーク内結合
        coupling_sum = 0.0
        for j in 1:N
            if j != i
                v_j = X[2j]
                coupling_sum += K[i, j] * (v_j - v_i)
            end
        end

        dX[2i-1] = du
        dX[2i] = dv + coupling_sum
    end

    return dX
end

"""
    network_vector_field(X, params::NetworkParams) -> Vector{Float64}

ネットワークのベクトル場 F(X) を返す（Pure関数）
内部結合を含む完全な力学
"""
network_vector_field(X::Vector{Float64}, params::NetworkParams) = network_full_dynamics(X, params)

# =============================================================================
# ネットワークヤコビアン
# =============================================================================

"""
    network_intrinsic_jacobian(X, params::NetworkParams) -> Matrix{Float64}

内部力学のみのヤコビアン（ブロック対角）（Pure関数）

# Returns
- `J::Matrix{Float64}`: (Ds*N) × (Ds*N) ヤコビアン行列
"""
function network_intrinsic_jacobian(X::Vector{Float64}, params::NetworkParams)
    N = params.N
    D = Ds * N
    J = zeros(D, D)

    for i in 1:N
        u = X[2i-1]
        v = X[2i]
        Ji = fhn_single_jacobian(u, v, params.fhn)

        # 2×2ブロックを埋める
        idx = 2i-1:2i
        J[idx, idx] = Ji
    end

    return J
end

"""
    network_full_jacobian(X, params::NetworkParams) -> Matrix{Float64}

完全なネットワークヤコビアン（内部結合を含む）（Pure関数）

# Returns
- `J::Matrix{Float64}`: (Ds*N) × (Ds*N) ヤコビアン行列
"""
function network_full_jacobian(X::Vector{Float64}, params::NetworkParams)
    N = params.N
    K = params.K
    D = Ds * N
    J = zeros(D, D)

    for i in 1:N
        u_i = X[2i-1]
        v_i = X[2i]

        # 内部力学のヤコビアン
        Ji = fhn_single_jacobian(u_i, v_i, params.fhn)
        idx_i = 2i-1:2i
        J[idx_i, idx_i] = Ji

        # 結合項のヤコビアン
        # dv_i/dt += Σ K_{ij} (v_j - v_i)
        # ∂(dv_i)/∂v_i += -Σ K_{ij}
        # ∂(dv_i)/∂v_j = K_{ij}

        for j in 1:N
            if j != i
                # v_i に対する微分（対角成分への寄与）
                J[2i, 2i] -= K[i, j]

                # v_j に対する微分
                J[2i, 2j] += K[i, j]
            end
        end
    end

    return J
end

# =============================================================================
# 3ネットワーク結合系の力学
# =============================================================================

"""
    three_network_dynamics(X_A, X_B, X_C, three_params::ThreeNetworkParams)
        -> (dX_A, dX_B, dX_C)

3つの結合ネットワークの力学（Pure関数）

位相縮約で使う3ネットワーク間の高次結合を含む

H_ijk(x_i, x_j, x_k) = (0, sin(v_j + v_k - 2*v_i))

# Arguments
- `X_A, X_B, X_C`: 各ネットワークの状態ベクトル
- `three_params::ThreeNetworkParams`: 3ネットワークパラメータ

# Returns
- `(dX_A, dX_B, dX_C)`: 各ネットワークの時間微分
"""
function three_network_dynamics(X_A::Vector{Float64}, X_B::Vector{Float64},
                                 X_C::Vector{Float64}, three_params::ThreeNetworkParams)
    params = three_params.network
    C = three_params.C
    epsilon = three_params.epsilon
    N = params.N

    # 各ネットワークの内部力学
    dX_A = network_full_dynamics(X_A, params)
    dX_B = network_full_dynamics(X_B, params)
    dX_C = network_full_dynamics(X_C, params)

    # 高次結合: H_ijk(x_i, x_j, x_k) = (0, sin(v_j + v_k - 2*v_i))
    # ネットワークA: ε Σ C_ijk H_ijk(x_A_i, x_B_j, x_C_k)
    for i in 1:N
        v_A_i = X_A[2i]
        coupling_A = 0.0
        coupling_B = 0.0
        coupling_C = 0.0

        for j in 1:N
            v_B_j = X_B[2j]
            v_A_j = X_A[2j]
            v_C_j = X_C[2j]

            for k in 1:N
                v_C_k = X_C[2k]
                v_A_k = X_A[2k]
                v_B_k = X_B[2k]

                # ネットワークAへの結合
                coupling_A += C[i, j, k] * sin(v_B_j + v_C_k - 2 * v_A_i)

                # ネットワークBへの結合（対称性: B受け取り A, C から）
                coupling_B += C[i, j, k] * sin(v_C_j + v_A_k - 2 * X_B[2i])

                # ネットワークCへの結合（対称性: C受け取り A, B から）
                coupling_C += C[i, j, k] * sin(v_A_j + v_B_k - 2 * X_C[2i])
            end
        end

        # v成分（偶数インデックス）に結合項を追加
        dX_A[2i] += epsilon * coupling_A
        dX_B[2i] += epsilon * coupling_B
        dX_C[2i] += epsilon * coupling_C
    end

    return (dX_A, dX_B, dX_C)
end

"""
    higher_order_interaction(v_i, v_j, v_k) -> (h_u, h_v)

高次相互作用関数 H_ijk（Pure関数）

H_ijk(x_i, x_j, x_k) = (0, sin(v_j + v_k - 2*v_i))

# Returns
- `(h_u, h_v)`: 相互作用ベクトル
"""
function higher_order_interaction(v_i::Float64, v_j::Float64, v_k::Float64)
    return (0.0, sin(v_j + v_k - 2 * v_i))
end

"""
    higher_order_interaction_jacobian(v_i, v_j, v_k) -> (∂H/∂v_i, ∂H/∂v_j, ∂H/∂v_k)

高次相互作用のヤコビアン（v成分に関する微分）（Pure関数）

# Returns
- `(dH_dvi, dH_dvj, dH_dvk)`: 各変数に関する微分
"""
function higher_order_interaction_jacobian(v_i::Float64, v_j::Float64, v_k::Float64)
    arg = v_j + v_k - 2 * v_i
    c = cos(arg)

    # ∂H_v/∂v_i = -2 * cos(v_j + v_k - 2*v_i)
    # ∂H_v/∂v_j = cos(v_j + v_k - 2*v_i)
    # ∂H_v/∂v_k = cos(v_j + v_k - 2*v_i)
    return (-2 * c, c, c)
end

# =============================================================================
# ヘルパー関数
# =============================================================================

"""
    extract_v_components(X::Vector{Float64}, N::Int) -> Vector{Float64}

状態ベクトルからv成分（膜電位）のみを抽出（Pure関数）
"""
function extract_v_components(X::Vector{Float64}, N::Int)
    return [X[2i] for i in 1:N]
end

"""
    extract_u_components(X::Vector{Float64}, N::Int) -> Vector{Float64}

状態ベクトルからu成分（回復変数）のみを抽出（Pure関数）
"""
function extract_u_components(X::Vector{Float64}, N::Int)
    return [X[2i-1] for i in 1:N]
end

"""
    create_state_vector(u::Vector{Float64}, v::Vector{Float64}) -> Vector{Float64}

u, vベクトルから状態ベクトルを構成（Pure関数）
"""
function create_state_vector(u::Vector{Float64}, v::Vector{Float64})
    @assert length(u) == length(v)
    N = length(u)
    X = zeros(Ds * N)
    for i in 1:N
        X[2i-1] = u[i]
        X[2i] = v[i]
    end
    return X
end

