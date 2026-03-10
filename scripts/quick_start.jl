#!/usr/bin/env julia
"""
    quick_start.jl - Ising2D クイックスタート

単一温度でのモンテカルロシミュレーションを実行し、物理量を計測する。

使用方法:
    julia --project scripts/quick_start.jl
"""

using Ising2D

function main()
    # =========================================================================
    # パラメータ設定
    # =========================================================================

    L = 32
    J = 1.0
    T = 2.269      # 臨界温度付近
    beta = 1.0 / T
    N = L * L

    n_equilibration = 1000
    n_sweeps = 5000
    sample_interval = 5

    println("Ising2D Monte Carlo Simulation")
    println("=" ^ 50)
    println("  L = $L, J = $J, T = $T")
    println("  equilibration = $n_equilibration MCS")
    println("  measurement   = $n_sweeps MCS (interval = $sample_interval)")

    # =========================================================================
    # 格子生成 + 熱平衡化
    # =========================================================================

    lattice = random_lattice(L)

    println("\n[1] Equilibrating...")
    for _ in 1:n_equilibration
        sweep!(lattice, beta; J=J)
    end
    println("  Done.")

    # =========================================================================
    # 測定
    # =========================================================================

    println("[2] Measuring...")

    avg = ThermalAverages()
    total_accepted = 0

    for sweep_i in 1:n_sweeps
        accepted = sweep!(lattice, beta; J=J)
        total_accepted += accepted

        if sweep_i % sample_interval == 0
            e = energy(lattice; J=J) / N
            m = magnetization(lattice)

            avg.sum_e += e
            avg.sum_e2 += e^2
            avg.sum_abs_m += abs(m)
            avg.sum_m2 += m^2
            avg.sum_m4 += m^4
            avg.n_samples += 1
        end
    end

    # =========================================================================
    # 熱力学量の計算
    # =========================================================================

    n = avg.n_samples
    mean_e = avg.sum_e / n
    mean_e2 = avg.sum_e2 / n
    mean_abs_m = avg.sum_abs_m / n
    mean_m2 = avg.sum_m2 / n
    mean_m4 = avg.sum_m4 / n

    C = N / T^2 * (mean_e2 - mean_e^2)
    chi = N / T * (mean_m2 - mean_abs_m^2)
    U_L = 1.0 - mean_m4 / (3.0 * mean_m2^2)

    acceptance_rate = total_accepted / (n_sweeps * N)

    # =========================================================================
    # 結果表示
    # =========================================================================

    println("\n" * "=" ^ 50)
    println("Results (n_samples = $n)")
    println("=" ^ 50)
    println("  ⟨e⟩   = $(round(mean_e, digits=6))")
    println("  ⟨|m|⟩ = $(round(mean_abs_m, digits=6))")
    println("  C      = $(round(C, digits=6))")
    println("  χ      = $(round(chi, digits=6))")
    println("  U_L    = $(round(U_L, digits=6))")
    println("  acceptance rate = $(round(acceptance_rate, digits=4))")
end

main()
