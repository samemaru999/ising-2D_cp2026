#!/usr/bin/env julia
"""
    benchmark.jl - ホットパスのベンチマーク

使用方法:
    julia --project scripts/benchmark.jl
"""

using Ising2D
using BenchmarkTools

function main()
    L = 32
    N = L * L
    beta = 1.0 / 2.269
    J = 1.0
    lattice = random_lattice(L)

    println("=" ^ 55)
    println("Ising2D Benchmark  (L=$L, N=$N, T=2.269)")
    println("=" ^ 55)

    # --- delta_energy ---
    println("\n[1] delta_energy (1回のΔE計算)")
    b1 = @benchmark delta_energy($lattice, 16, 16; J=$J)
    display(b1)

    # --- metropolis_step! ---
    println("\n\n[2] metropolis_step! (1スピン更新試行)")
    b2 = @benchmark metropolis_step!($lattice, $beta; J=$J)
    display(b2)

    # --- sweep! ---
    println("\n\n[3] sweep! (1スイープ = $(N)回の更新試行)")
    b3 = @benchmark sweep!($lattice, $beta; J=$J)
    display(b3)

    # --- 推定 ---
    t_sweep_ns = minimum(b3).time  # ナノ秒
    t_per_spin = t_sweep_ns / N
    println("\n\n" * "=" ^ 55)
    println("Summary")
    println("=" ^ 55)
    println("  sweep! minimum  : $(round(t_sweep_ns / 1000, digits=1)) μs")
    println("  per spin update : $(round(t_per_spin, digits=1)) ns")
    println("  (参考: 1GHz CPU = 1ns/cycle, 現代CPU ≈ 0.2-0.3ns/cycle)")
    if t_per_spin > 500
        println("  → 500ns超: 改善の余地が大きい")
    elseif t_per_spin > 100
        println("  → 100-500ns: 標準的、最適化で数倍改善可能")
    elseif t_per_spin > 30
        println("  → 30-100ns: 良好")
    else
        println("  → 30ns以下: 非常に高速")
    end
end

main()
