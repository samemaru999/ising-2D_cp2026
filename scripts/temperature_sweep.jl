#!/usr/bin/env julia
"""
    temperature_sweep.jl - 温度掃引シミュレーション

T = 1.5 ~ 3.5 を 0.05 刻みで掃引し、相転移を観測する。
Tc ≈ 2.269 付近で比熱・帯磁率のピークが見られる。

使用方法:
    julia --project scripts/temperature_sweep.jl
"""

using Ising2D
using Dates

function main()
    # =========================================================================
    # パラメータ設定
    # =========================================================================

    L = 32
    J = 1.0
    temperatures = collect(1.5:0.05:3.5)

    config = SimulationConfig(5000, 2000, 5)
    params = IsingParams(L, J)

    println("Ising2D Temperature Sweep")
    println("=" ^ 50)
    println("  L = $L, J = $J")
    println("  T = $(first(temperatures)) ~ $(last(temperatures)), ΔT = 0.05 ($(length(temperatures)) points)")
    println("  equilibration = $(config.n_equilibration) MCS")
    println("  measurement   = $(config.n_sweeps) MCS (interval = $(config.sample_interval))")

    # =========================================================================
    # 温度掃引実行
    # =========================================================================

    println("\nRunning temperature sweep...")
    result = run_temperature_sweep(params, config, temperatures)
    println("Done.")

    # =========================================================================
    # 結果サマリ
    # =========================================================================

    println("\n" * "=" ^ 50)
    println("Results Summary")
    println("=" ^ 50)
    println("  T        ⟨e⟩        ⟨|m|⟩      C          χ")
    println("  " * "-" ^ 60)
    for r in result.results
        t = r.thermodynamics
        Printf_T = round(r.temperature, digits=3)
        Printf_e = round(t.mean_energy, digits=4)
        Printf_m = round(t.mean_abs_magnetization, digits=4)
        Printf_C = round(t.specific_heat, digits=4)
        Printf_chi = round(t.susceptibility, digits=4)
        println("  $(lpad(Printf_T, 6))  $(lpad(Printf_e, 10))  $(lpad(Printf_m, 8))  $(lpad(Printf_C, 10))  $(lpad(Printf_chi, 10))")
    end

    # =========================================================================
    # データ保存（JLD2）
    # =========================================================================

    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    outdir = joinpath(@__DIR__, "..", "data", "output", timestamp)
    mkpath(outdir)
    filename = joinpath(outdir, "sweep_L$(L).jld2")
    save_result(filename, result)
    println("\nSaved to $filename")

    # =========================================================================
    # プロット
    # =========================================================================

    println("Generating plots...")
    fig = plot_thermodynamics(result)

    figdir = joinpath(@__DIR__, "..", "data", "output", "figs")
    mkpath(figdir)
    figpath = joinpath(figdir, "thermodynamics_L$(L).png")
    save(figpath, fig)
    println("Saved plot to $figpath")
end

main()
