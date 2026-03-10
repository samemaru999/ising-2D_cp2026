"""
シミュレーション結果の可視化（Impure）。
GLMakieを使用する。
"""

using GLMakie

"""
    plot_spin_lattice(lattice::Matrix{Int}) -> Figure

スピン配置のヒートマップを表示する。(F-VIS-01)
"""
function plot_spin_lattice(lattice::Matrix{Int})
    fig = Figure()
    ax = Axis(fig[1, 1]; title="Spin Configuration",
              xlabel="j", ylabel="i", aspect=DataAspect())
    heatmap!(ax, lattice; colormap=[:blue, :red])
    return fig
end

"""
    plot_thermodynamics(result::TemperatureSweepResult) -> Figure

物理量の温度依存性をプロットする。(F-VIS-02)
e(T), |m|(T), C(T), χ(T) の4パネル。
"""
function plot_thermodynamics(result::TemperatureSweepResult)
    T = result.temperatures
    e = [r.thermodynamics.mean_energy for r in result.results]
    m = [r.thermodynamics.mean_abs_magnetization for r in result.results]
    C = [r.thermodynamics.specific_heat for r in result.results]
    chi = [r.thermodynamics.susceptibility for r in result.results]

    fig = Figure(size=(800, 600))

    ax1 = Axis(fig[1, 1]; xlabel="T", ylabel="⟨e⟩", title="Energy")
    scatter!(ax1, T, e; markersize=6)

    ax2 = Axis(fig[1, 2]; xlabel="T", ylabel="⟨|m|⟩", title="Magnetization")
    scatter!(ax2, T, m; markersize=6)

    ax3 = Axis(fig[2, 1]; xlabel="T", ylabel="C", title="Specific Heat")
    scatter!(ax3, T, C; markersize=6)

    ax4 = Axis(fig[2, 2]; xlabel="T", ylabel="χ", title="Susceptibility")
    scatter!(ax4, T, chi; markersize=6)

    return fig
end

"""
    plot_binder_cumulant(results::Vector{TemperatureSweepResult}) -> Figure

複数格子サイズのビンダーキュムラントを重ね描きする。(F-VIS-03)
"""
function plot_binder_cumulant(results::Vector{TemperatureSweepResult})
    fig = Figure()
    ax = Axis(fig[1, 1]; xlabel="T", ylabel="U_L", title="Binder Cumulant")

    for result in results
        T = result.temperatures
        U = [r.thermodynamics.binder_cumulant for r in result.results]
        L = result.params.L
        scatter!(ax, T, U; label="L=$L", markersize=6)
    end
    axislegend(ax)

    return fig
end

"""
    plot_timeseries(energies::Vector{Float64}, magnetizations::Vector{Float64}) -> Figure

エネルギー・磁化のモンテカルロ時系列をプロットする。(F-VIS-05)
"""
function plot_timeseries(energies::Vector{Float64}, magnetizations::Vector{Float64})
    fig = Figure(size=(800, 400))
    steps = 1:length(energies)

    ax1 = Axis(fig[1, 1]; xlabel="MCS", ylabel="e", title="Energy")
    lines!(ax1, steps, energies)

    ax2 = Axis(fig[2, 1]; xlabel="MCS", ylabel="m", title="Magnetization")
    lines!(ax2, steps, magnetizations)

    return fig
end
