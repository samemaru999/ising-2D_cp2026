"""
2次元イジングモデルの格子状態を生成する関数群。
スピンは +1（上向き）と -1（下向き）の2値をとる。
"""

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
