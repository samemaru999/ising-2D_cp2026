"""
2次元イジングモデルの物理量を計算する純粋関数群。
"""

using Statistics

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
        E += s * lattice[mod1(i + 1, L), j]
        E += s * lattice[i, mod1(j + 1, L)]
    end
    return -J * E
end

"""
    delta_energy(lattice::Matrix{Int}, i::Int, j::Int; J::Float64=1.0) -> Float64

サイト (i,j) のスピンを反転した場合のエネルギー変化 ΔE を計算する。
ΔE = 2J × s(i,j) × (4近傍スピンの和)
"""
function delta_energy(lattice::Matrix{Int}, i::Int, j::Int; J::Float64=1.0)
    L = size(lattice, 1)
    s = lattice[i, j]
    neighbors = lattice[mod1(i - 1, L), j] +
                lattice[mod1(i + 1, L), j] +
                lattice[i, mod1(j - 1, L)] +
                lattice[i, mod1(j + 1, L)]
    return 2.0 * J * s * neighbors
end
