"""
メトロポリス法によるモンテカルロ更新。
乱数を使用するため Impure。
"""

"""
    metropolis_step!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0) -> Bool

ランダムに1サイトを選び、メトロポリス法でスピン反転を試みる。
反転が受容された場合 true を返す。
"""
function metropolis_step!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0)
    L = size(lattice, 1)
    i = rand(1:L)
    j = rand(1:L)

    dE = delta_energy(lattice, i, j; J=J)

    if dE <= 0 || rand() < exp(-beta * dE)
        lattice[i, j] *= -1
        return true
    end
    return false
end

"""
    sweep!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0) -> Int

L×L 回のメトロポリス更新試行を行う（1モンテカルロスイープ）。
受容された回数を返す。
"""
function sweep!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0)
    L = size(lattice, 1)
    accepted = 0
    for _ in 1:L*L
        accepted += metropolis_step!(lattice, beta; J=J)
    end
    return accepted
end
