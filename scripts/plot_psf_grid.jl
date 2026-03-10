#!/usr/bin/env julia
"""
    plot_psf_grid.jl - 位相感受関数の2x2グリッドプロット

各振動子のv成分 Q_v(θ) を2x2レイアウトで表示

使用方法:
    julia --project=. src/App/plot_psf_grid.jl
"""

using JLD2
using GLMakie

# PSFデータを読み込み
psf_file = "data/psf/psf_N4_Nth101_20251217_211050.jld2"
data = JLD2.load(psf_file)
Q = data["Q"]

Ntheta = size(Q, 2)
N = 4

# θの範囲
theta = range(0, 2π, length=Ntheta)

# 白背景テーマ
set_theme!(theme_light())

fig = Figure(size=(1200, 800), backgroundcolor=:white)

# 2x2グリッドで各振動子のv成分をプロット
for i in 1:N
    row = (i - 1) ÷ 2 + 1
    col = (i - 1) % 2 + 1

    ax = Axis(fig[row, col],
        title = "v_$i",
        xlabel = "θ",
        ylabel = "Qᵥ(θ)",
        xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3/2π", "2π"]),
        backgroundcolor = :white,
        xgridvisible = true,
        ygridvisible = true,
        xgridcolor = (:gray, 0.3),
        ygridcolor = (:gray, 0.3)
    )

    # v成分のインデックス
    v_idx = 2 * i
    Q_v = Q[v_idx, :]

    lines!(ax, theta, Q_v, linewidth=2, color=:steelblue)
end

display(fig)

println("PSF plot displayed!")
println("Press Enter to close...")
readline()
