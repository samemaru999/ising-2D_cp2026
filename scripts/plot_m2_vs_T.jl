#!/usr/bin/env julia
"""
    plot_m2_vs_T.jl - ⟨M²⟩ vs T プロット

quick_mc を用いて温度スイープを行い、⟨M²⟩ の温度依存性をプロットする。

使用方法:
    julia --project scripts/plot_m2_vs_T.jl
    julia --project scripts/plot_m2_vs_T.jl 16    # L=16 を指定
"""

using Ising2D
using GLMakie

function main()
    L = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 32

    temperatures = collect(1.5:0.05:3.5)
    m2_values = Float64[]

    println("Running ⟨M²⟩ vs T  (L=$L, $(length(temperatures)) points)")
    for (i, T) in enumerate(temperatures)
        result = quick_mc(T, L)
        push!(m2_values, result["M2"])
        println("  [$i/$(length(temperatures))] T = $(round(T, digits=2))  ⟨M²⟩ = $(round(result["M2"], digits=6))")
    end

    # Plot
    fig = Figure(size=(700, 500))
    ax = Axis(fig[1, 1];
              title="⟨M²⟩ vs T  (L=$L)",
              xlabel="Temperature T",
              ylabel="⟨M²⟩")
    Tc = 2 / log(1 + sqrt(2))  # exact Tc ≈ 2.269
    vlines!(ax, [Tc]; color=:gray, linestyle=:dash, label="Tc ≈ $(round(Tc, digits=3))")
    scatterlines!(ax, temperatures, m2_values; markersize=6, label="L=$L")
    axislegend(ax; position=:rt)

    # Save
    figdir = joinpath(@__DIR__, "..", "data", "output", "figs")
    mkpath(figdir)
    filepath = joinpath(figdir, "m2_vs_T.png")
    save(filepath, fig)
    println("\nFigure saved to $filepath")
end

main()
