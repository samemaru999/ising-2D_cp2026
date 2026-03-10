#!/usr/bin/env julia
"""
    run_simulation.jl - 温度掃引シミュレーション (CLI対応)

追加した機能 (run_temperature_sweep, compute_thermodynamics) を
コマンドライン引数で柔軟に実行する。

使用方法:
    # デフォルト (L=32, T=1.5:0.05:3.5)
    julia --project scripts/run_simulation.jl

    # 格子サイズ指定
    julia --project scripts/run_simulation.jl --L 16

    # 温度範囲指定
    julia --project scripts/run_simulation.jl --L 64 --Tmin 2.0 --Tmax 2.6 --dT 0.02

    # MCSパラメータ指定
    julia --project scripts/run_simulation.jl --n_sweeps 10000 --n_eq 5000

    # プロットなし (データ保存のみ)
    julia --project scripts/run_simulation.jl --no-plot
"""

using Ising2D
using Dates
using GLMakie

function parse_args(args)
    defaults = Dict(
        :L => 32,
        :J => 1.0,
        :Tmin => 1.5,
        :Tmax => 3.5,
        :dT => 0.05,
        :n_sweeps => 5000,
        :n_eq => 2000,
        :sample_interval => 5,
        :plot => true,
    )

    i = 1
    while i <= length(args)
        key = args[i]
        if key == "--L"
            defaults[:L] = parse(Int, args[i+1]); i += 2
        elseif key == "--J"
            defaults[:J] = parse(Float64, args[i+1]); i += 2
        elseif key == "--Tmin"
            defaults[:Tmin] = parse(Float64, args[i+1]); i += 2
        elseif key == "--Tmax"
            defaults[:Tmax] = parse(Float64, args[i+1]); i += 2
        elseif key == "--dT"
            defaults[:dT] = parse(Float64, args[i+1]); i += 2
        elseif key == "--n_sweeps"
            defaults[:n_sweeps] = parse(Int, args[i+1]); i += 2
        elseif key == "--n_eq"
            defaults[:n_eq] = parse(Int, args[i+1]); i += 2
        elseif key == "--sample_interval"
            defaults[:sample_interval] = parse(Int, args[i+1]); i += 2
        elseif key == "--no-plot"
            defaults[:plot] = false; i += 1
        elseif key in ["--help", "-h"]
            println(@doc run_simulation)
            exit(0)
        else
            error("Unknown argument: $key\nRun with --help for usage.")
        end
    end

    return defaults
end

function main()
    opts = parse_args(ARGS)

    L = opts[:L]
    J = opts[:J]
    temperatures = collect(opts[:Tmin]:opts[:dT]:opts[:Tmax])

    params = IsingParams(L, J)
    config = SimulationConfig(opts[:n_sweeps], opts[:n_eq], opts[:sample_interval])

    # =========================================================================
    # パラメータ表示
    # =========================================================================

    println("Ising2D Temperature Sweep Simulation")
    println("=" ^ 55)
    println("  L = $L, J = $J, N = $(L^2)")
    println("  T = $(opts[:Tmin]) : $(opts[:dT]) : $(opts[:Tmax])  ($(length(temperatures)) points)")
    println("  equilibration = $(config.n_equilibration) MCS")
    println("  measurement   = $(config.n_sweeps) MCS (interval = $(config.sample_interval))")
    println("  plot          = $(opts[:plot])")

    # =========================================================================
    # 温度掃引実行
    # =========================================================================

    println("\nRunning temperature sweep...")
    t_start = time()
    result = run_temperature_sweep(params, config, temperatures)
    elapsed = round(time() - t_start, digits=1)
    println("Done. ($(elapsed)s)")

    # =========================================================================
    # 結果サマリ
    # =========================================================================

    println("\n" * "=" ^ 55)
    println("Results Summary")
    println("-" ^ 55)
    println("  T       ⟨e⟩        ⟨|m|⟩     C          χ")
    println("  " * "-" ^ 51)
    for r in result.results
        t = r.thermodynamics
        println("  $(lpad(round(r.temperature, digits=3), 6))" *
                "  $(lpad(round(t.mean_energy, digits=4), 10))" *
                "  $(lpad(round(t.mean_abs_magnetization, digits=4), 8))" *
                "  $(lpad(round(t.specific_heat, digits=4), 10))" *
                "  $(lpad(round(t.susceptibility, digits=4), 10))")
    end

    # =========================================================================
    # データ保存
    # =========================================================================

    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    outdir = joinpath(@__DIR__, "..", "data", "output", timestamp)
    mkpath(outdir)
    filename = joinpath(outdir, "sweep_L$(L).jld2")
    save_result(filename, result)
    println("\nData saved to $filename")

    # =========================================================================
    # プロット (--no-plot でスキップ可能)
    # =========================================================================

    if opts[:plot]
        println("Generating plots...")
        fig = plot_thermodynamics(result)

        figdir = joinpath(@__DIR__, "..", "data", "output", "figs")
        mkpath(figdir)
        figpath = joinpath(figdir, "thermodynamics_L$(L)_$(timestamp).png")
        save(figpath, fig)
        println("Plot saved to $figpath")
    end

    println("\nAll done.")
end

main()
