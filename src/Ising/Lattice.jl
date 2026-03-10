"""
2次元イジングモデルの格子状態を生成する関数群。
スピンは +1（上向き）と -1（下向き）の2値をとる。
"""

using Statistics

"""
    random_lattice(L::Int) -> Matrix{Int}

L×L のランダムスピン配置を生成する。
各サイトに +1 または -1 をランダムに割り当てる。
"""
function random_lattice(L::Int)
    return rand([-1, 1], L, L)
end

"""
    uniform_lattice(L::Int; spin::Int=1) -> Matrix{Int}

L×L の一様スピン配置を生成する。
全サイトが同じスピン値（デフォルト +1）を持つ。
"""
function uniform_lattice(L::Int; spin::Int=1)
    @assert spin == 1 || spin == -1 "spin must be +1 or -1"
    return fill(spin, L, L)
end

# =============================================================================
# 物理量
# =============================================================================

"""
    magnetization(lattice::Matrix{Int}) -> Float64

1サイトあたりの磁化 M = (1/N) Σ s_i を計算する。
M ∈ [-1, 1] で、|M| ≈ 1 なら秩序相、M ≈ 0 なら無秩序相。
"""
function magnetization(lattice::Matrix{Int})
    return mean(lattice)
end


"""
    energy(lattice::Matrix{Int}; J::Float64=1.0) -> Float64

H = -J Σ_{⟨i,j⟩} s_i s_j を周期境界条件で計算する。
各ボンドは1回だけカウントする（右隣・下隣のみ）。
"""
function energy(lattice::Matrix{Int}; J::Float64=1.0)
    L = size(lattice, 1)
    E = 0
    for j in 1:L, i in 1:L
        s = lattice[i, j]
        # 右隣と下隣だけ見る → 各ボンドを1回だけカウント
        E += s * lattice[mod1(i + 1, L), j]
        E += s * lattice[i, mod1(j + 1, L)]
    end
    return -J * E
end
