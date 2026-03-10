"""
    Logic/Coupling.jl - 結合関数 [Pure]

ネットワーク内結合（ペアワイズ）と
ネットワーク間結合（高次相互作用）を定義する純粋関数群
"""

using LinearAlgebra

# =============================================================================
# ネットワーク内結合（拡散結合）
# =============================================================================

"""
    diffusive_coupling(v_i, v_j, K_ij) -> Float64

拡散結合項を計算（Pure関数）

g_ij(x_i, x_j) = K_ij * (v_j - v_i)

# Arguments
- `v_i::Float64`: 振動子iの膜電位
- `v_j::Float64`: 振動子jの膜電位
- `K_ij::Float64`: 結合強度

# Returns
- `Float64`: v成分への結合寄与
"""
function diffusive_coupling(v_i::Float64, v_j::Float64, K_ij::Float64)
    return K_ij * (v_j - v_i)
end

"""
    compute_intra_coupling(X, K, N) -> Vector{Float64}

ネットワーク内結合項全体を計算（Pure関数）

# Arguments
- `X::Vector{Float64}`: ネットワーク状態 (長さ Ds*N)
- `K::Matrix{Float64}`: N×N結合行列
- `N::Int`: 振動子数

# Returns
- `coupling::Vector{Float64}`: 結合寄与ベクトル (長さ Ds*N)
"""
function compute_intra_coupling(X::Vector{Float64}, K::Matrix{Float64}, N::Int)
    coupling = zeros(Ds * N)

    for i in 1:N
        v_i = X[2i]
        coupling_sum = 0.0
        for j in 1:N
            if i != j
                v_j = X[2j]
                coupling_sum += diffusive_coupling(v_i, v_j, K[i, j])
            end
        end
        coupling[2i] = coupling_sum  # v成分のみに影響
    end

    return coupling
end

"""
    intra_coupling_jacobian(K, N) -> Matrix{Float64}

ネットワーク内結合のヤコビアン寄与を計算（Pure関数）

結合 g_ij = K_ij(v_j - v_i) の寄与:
- ∂g_i/∂v_i = -Σ_j K_ij
- ∂g_i/∂v_j = K_ij

# Returns
- `J::Matrix{Float64}`: (Ds*N) × (Ds*N) ヤコビアン寄与
"""
function intra_coupling_jacobian(K::Matrix{Float64}, N::Int)
    D = Ds * N
    J = zeros(D, D)

    for i in 1:N
        v_i_idx = 2i

        for j in 1:N
            if i != j
                # ∂(dv_i)/∂v_i への寄与
                J[v_i_idx, v_i_idx] -= K[i, j]
                # ∂(dv_i)/∂v_j への寄与
                J[v_i_idx, 2j] += K[i, j]
            end
        end
    end

    return J
end

# =============================================================================
# ネットワーク間高次結合（3体相互作用）
# =============================================================================

"""
    higher_order_H(v_i, v_j, v_k) -> Float64

高次相互作用関数のv成分（Pure関数）

H_ijk(x_i, x_j, x_k) = (0, sin(v_j + v_k - 2*v_i))^T

# Returns
- `Float64`: 相互作用のv成分
"""
function higher_order_H(v_i::Float64, v_j::Float64, v_k::Float64)
    return sin(v_j + v_k - 2 * v_i)
end

"""
    higher_order_H_full(v_i, v_j, v_k) -> Tuple{Float64, Float64}

高次相互作用関数の完全形式（u, v両成分）（Pure関数）

H_ijk(x_i, x_j, x_k) = (0, sin(v_j + v_k - 2*v_i))^T
"""
function higher_order_H_full(v_i::Float64, v_j::Float64, v_k::Float64)
    return (0.0, sin(v_j + v_k - 2 * v_i))
end

"""
    higher_order_H_derivatives(v_i, v_j, v_k) -> (dH_dvi, dH_dvj, dH_dvk)

高次相互作用のv成分に関する偏微分（Pure関数）

∂H_v/∂v_i = -2 * cos(v_j + v_k - 2*v_i)
∂H_v/∂v_j = cos(v_j + v_k - 2*v_i)
∂H_v/∂v_k = cos(v_j + v_k - 2*v_i)
"""
function higher_order_H_derivatives(v_i::Float64, v_j::Float64, v_k::Float64)
    arg = v_j + v_k - 2 * v_i
    c = cos(arg)
    return (-2 * c, c, c)
end

"""
    compute_interaction_tensor(X_self, X_source1, X_source2, N) -> Array{Float64,3}

高次相互作用テンソル H[i,j,k] を計算（Pure関数）

# Arguments
- `X_self::Vector{Float64}`: 受け取りネットワークの状態
- `X_source1::Vector{Float64}`: 第1ソースネットワークの状態
- `X_source2::Vector{Float64}`: 第2ソースネットワークの状態
- `N::Int`: 振動子数

# Returns
- `H::Array{Float64,3}`: N×N×N 相互作用テンソル
"""
function compute_interaction_tensor(X_self::Vector{Float64}, X_source1::Vector{Float64},
                                     X_source2::Vector{Float64}, N::Int)
    H = zeros(N, N, N)

    for i in 1:N
        v_i = X_self[2i]
        for j in 1:N
            v_j = X_source1[2j]
            for k in 1:N
                v_k = X_source2[2k]
                H[i, j, k] = higher_order_H(v_i, v_j, v_k)
            end
        end
    end

    return H
end

"""
    compute_inter_coupling(X_self, X_source1, X_source2, C, epsilon, N) -> Vector{Float64}

ネットワーク間高次結合項を計算（Pure関数）

# Arguments
- `X_self::Vector{Float64}`: 受け取りネットワークの状態
- `X_source1::Vector{Float64}`: 第1ソースネットワークの状態
- `X_source2::Vector{Float64}`: 第2ソースネットワークの状態
- `C::Array{Float64,3}`: N×N×N 結合テンソル
- `epsilon::Float64`: 結合強度
- `N::Int`: 振動子数

# Returns
- `coupling::Vector{Float64}`: 結合寄与ベクトル (長さ Ds*N)
"""
function compute_inter_coupling(X_self::Vector{Float64}, X_source1::Vector{Float64},
                                 X_source2::Vector{Float64}, C::Array{Float64,3},
                                 epsilon::Float64, N::Int)
    coupling = zeros(Ds * N)

    for i in 1:N
        v_i = X_self[2i]
        coupling_sum = 0.0

        for j in 1:N
            v_j = X_source1[2j]
            for k in 1:N
                v_k = X_source2[2k]
                H_ijk = higher_order_H(v_i, v_j, v_k)
                coupling_sum += C[i, j, k] * H_ijk
            end
        end

        coupling[2i] = epsilon * coupling_sum
    end

    return coupling
end

# =============================================================================
# 3ネットワーク結合系の完全力学
# =============================================================================

"""
    ThreeNetworkState

3ネットワークの状態を保持する構造体
"""
struct ThreeNetworkState
    X_A::Vector{Float64}
    X_B::Vector{Float64}
    X_C::Vector{Float64}
end

"""
    compute_three_network_coupling(state::ThreeNetworkState, C, epsilon, N)
        -> (coupling_A, coupling_B, coupling_C)

3ネットワーク間の高次結合を計算（Pure関数）

巡回結合構造:
- ネットワークA: B, Cから受け取る
- ネットワークB: C, Aから受け取る
- ネットワークC: A, Bから受け取る

# Returns
- 各ネットワークへの結合寄与ベクトルのタプル
"""
function compute_three_network_coupling(state::ThreeNetworkState, C::Array{Float64,3},
                                         epsilon::Float64, N::Int)
    coupling_A = compute_inter_coupling(state.X_A, state.X_B, state.X_C, C, epsilon, N)
    coupling_B = compute_inter_coupling(state.X_B, state.X_C, state.X_A, C, epsilon, N)
    coupling_C = compute_inter_coupling(state.X_C, state.X_A, state.X_B, C, epsilon, N)

    return (coupling_A, coupling_B, coupling_C)
end

"""
    compute_three_network_dynamics(state::ThreeNetworkState, params::ThreeNetworkParams)
        -> (dX_A, dX_B, dX_C)

3ネットワーク結合系の完全な力学を計算（Pure関数）

# Returns
- 各ネットワークの時間微分ベクトルのタプル
"""
function compute_three_network_dynamics(state::ThreeNetworkState,
                                         three_params::ThreeNetworkParams)
    params = three_params.network
    C = three_params.C
    epsilon = three_params.epsilon
    N = params.N

    # 各ネットワークの内部力学（ネットワーク内結合を含む）
    dX_A = network_full_dynamics(state.X_A, params)
    dX_B = network_full_dynamics(state.X_B, params)
    dX_C = network_full_dynamics(state.X_C, params)

    # ネットワーク間高次結合を追加
    coupling_A, coupling_B, coupling_C = compute_three_network_coupling(state, C, epsilon, N)

    dX_A .+= coupling_A
    dX_B .+= coupling_B
    dX_C .+= coupling_C

    return (dX_A, dX_B, dX_C)
end

# =============================================================================
# 結合行列・テンソルの生成関数
# =============================================================================

"""
    create_full_coupling_matrix(N; strength=1.0) -> Matrix{Float64}

全結合行列を生成（Pure関数）

対角成分は0（自己結合なし）
"""
function create_full_coupling_matrix(N::Int; strength::Float64=1.0)
    K = fill(strength, N, N)
    for i in 1:N
        K[i, i] = 0.0
    end
    return K
end

"""
    create_ring_coupling_matrix(N; strength=1.0) -> Matrix{Float64}

リング結合行列を生成（Pure関数）

各振動子は隣接する2つの振動子とのみ結合
"""
function create_ring_coupling_matrix(N::Int; strength::Float64=1.0)
    K = zeros(N, N)
    for i in 1:N
        j_prev = mod1(i - 1, N)
        j_next = mod1(i + 1, N)
        K[i, j_prev] = strength
        K[i, j_next] = strength
    end
    return K
end

"""
    create_star_coupling_matrix(N; strength=1.0) -> Matrix{Float64}

スター型結合行列を生成（Pure関数）

振動子1がハブ、他はハブとのみ結合
"""
function create_star_coupling_matrix(N::Int; strength::Float64=1.0)
    K = zeros(N, N)
    for i in 2:N
        K[1, i] = strength  # ハブから周辺へ
        K[i, 1] = strength  # 周辺からハブへ
    end
    return K
end

"""
    create_diagonal_coupling_tensor(N; strength=0.0112) -> Array{Float64,3}

対角型高次結合テンソルを生成（Pure関数）

C[i, i, j] = strength for all i, j
（既存研究の一様結合に対応）
"""
function create_diagonal_coupling_tensor(N::Int; strength::Float64=0.0112)
    C = zeros(N, N, N)
    for i in 1:N
        for j in 1:N
            C[i, i, j] = strength
        end
    end
    return C
end

"""
    create_antisymmetric_coupling_tensor(N, base_tensor) -> Array{Float64,3}

反対称結合テンソルを生成（Pure関数）

C_antisym[i,j,k] = (base[i,j,k] - base[i,k,j]) / 2
"""
function create_antisymmetric_coupling_tensor(N::Int, base_tensor::Array{Float64,3})
    C = zeros(N, N, N)
    for i in 1:N
        for j in 1:N
            for k in 1:N
                C[i, j, k] = (base_tensor[i, j, k] - base_tensor[i, k, j]) / 2
            end
        end
    end
    return C
end

# =============================================================================
# 結合特性の計算
# =============================================================================

"""
    compute_coupling_frobenius_norm(C) -> Float64

結合テンソルのフロベニウスノルムを計算（Pure関数）
"""
function compute_coupling_frobenius_norm(C::Array{Float64,3})
    return sqrt(sum(C .^ 2))
end

"""
    compute_coupling_sum(C) -> Float64

結合テンソルの総和を計算（Pure関数）
"""
function compute_coupling_sum(C::Array{Float64,3})
    return sum(C)
end

"""
    is_antisymmetric(C; tol=1e-10) -> Bool

結合テンソルが反対称か判定（Pure関数）

C[i,j,k] = -C[i,k,j] for all i,j,k
"""
function is_antisymmetric(C::Array{Float64,3}; tol::Float64=1e-10)
    N = size(C, 1)
    for i in 1:N
        for j in 1:N
            for k in 1:N
                if abs(C[i, j, k] + C[i, k, j]) > tol
                    return false
                end
            end
        end
    end
    return true
end

"""
    symmetrize_coupling_tensor(C) -> Array{Float64,3}

結合テンソルを対称化（Pure関数）

C_sym[i,j,k] = (C[i,j,k] + C[i,k,j]) / 2
"""
function symmetrize_coupling_tensor(C::Array{Float64,3})
    N = size(C, 1)
    C_sym = zeros(N, N, N)
    for i in 1:N
        for j in 1:N
            for k in 1:N
                C_sym[i, j, k] = (C[i, j, k] + C[i, k, j]) / 2
            end
        end
    end
    return C_sym
end

