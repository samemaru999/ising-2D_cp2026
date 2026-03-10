#!/usr/bin/env julia
"""
    plot_results.jl - 保存済みシミュレーション結果の可視化

data/output/ 内のJLD2ファイルを読み込み、物理量をプロットする。

使用方法:
    julia --project scripts/plot_results.jl <jld2ファイルパス>
    julia --project scripts/plot_results.jl data/output/20260310_153042/L32_T2.269.jld2
"""

using Ising2D
using JLD2
using GLMakie
using Dates

function main()
    if length(ARGS) < 1
        println("Usage: julia --project scripts/plot_results.jl <path-to-jld2>")
        return
    end

    filepath = ARGS[1]
    isfile(filepath) || error("File not found: $filepath")

    println("Loading: $filepath")
    data = load_result(filepath)

    # SingleTemperatureResult の場合（temperature キーで判別）
    if haskey(data, "temperature")
        plot_single(data)
    # TemperatureSweepResult の場合（temperatures キーで判別）
    elseif haskey(data, "temperatures")
        plot_sweep(data)
    else
        error("Unknown data format. Keys: $(keys(data))")
    end
end

function plot_single(data)
    T = data["temperature"]
    thermo = data["thermodynamics"]
    params = data["params"]
    rate = data["acceptance_rate"]

    println("\n  L = $(params.L), J = $(params.J), T = $T")
    println("  ⟨e⟩   = $(round(thermo.mean_energy, digits=6))")
    println("  ⟨|m|⟩ = $(round(thermo.mean_abs_magnetization, digits=6))")
    println("  C      = $(round(thermo.specific_heat, digits=6))")
    println("  χ      = $(round(thermo.susceptibility, digits=6))")
    println("  U_L    = $(round(thermo.binder_cumulant, digits=6))")
    println("  acceptance rate = $(round(rate, digits=4))")

    fig = Figure(size=(600, 400))
    ax = Axis(fig[1, 1]; title="L=$(params.L), T=$T",
              xlabel="Quantity", ylabel="Value")
    labels = ["⟨e⟩", "⟨|m|⟩", "C", "χ", "U_L"]
    values = [thermo.mean_energy, thermo.mean_abs_magnetization,
              thermo.specific_heat, thermo.susceptibility, thermo.binder_cumulant]
    barplot!(ax, 1:5, values; bar_labels=labels)
    ax.xticks = (1:5, labels)

    save_fig(fig, "single_L$(params.L)_T$(round(T, digits=3))")
    display(fig)
    println("\nPress Enter to close...")
    readline()
end

function plot_sweep(data)
    params = data["params"]
    temperatures = data["temperatures"]
    results = data["results"]

    sweep_result = TemperatureSweepResult(params, data["config"], temperatures, results)
    fig = plot_thermodynamics(sweep_result)

    save_fig(fig, "sweep_L$(params.L)")
    display(fig)
    println("\nPress Enter to close...")
    readline()
end

function save_fig(fig::Figure, name::String)
    figdir = joinpath(@__DIR__, "..", "data", "output", "figs")
    mkpath(figdir)
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    filepath = joinpath(figdir, "$(name)_$(timestamp).png")
    save(filepath, fig)
    println("  Figure saved to $filepath")
end

main()
