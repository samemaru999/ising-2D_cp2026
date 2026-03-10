#!/usr/bin/env julia
"""
    plot_pcf_3d.jl - 位相結合関数の3Dサーフェスプロット

Γ(φ, ψ) を3D表面として表示

使用方法:
    julia --project=. scripts/plot_pcf_3d.jl
"""

using JLD2
using GLMakie

# PCFデータを読み込み
pcf_file = "data/pcf/total_pcf_linear_N4_Nth101_20251217_211050.jld2"
data = JLD2.load(pcf_file)
Gamma = data["Gamma"]

Ntheta = size(Gamma, 1)

# φ, ψの範囲を -π から π に設定
phi_range = range(-π, π, length=Ntheta)
psi_range = range(-π, π, length=Ntheta)

# Γ行列をシフト（0~2π → -π~π）
half = Ntheta ÷ 2
Gamma_shifted = circshift(Gamma, (half, half))

# Figure作成
fig = Figure(size=(900, 800), backgroundcolor=:white)

# 3Dサーフェスプロット
ax = Axis3(fig[1, 1],
    title = "Linear stability",
    xlabel = "φ",
    ylabel = "φ",
    zlabel = "Γ(φ, φ)",
    xticks = ([-π, -π/2, 0, π/2, π], ["-π", "-π/2", "0", "π/2", "π"]),
    yticks = ([-π, -π/2, 0, π/2, π], ["-π", "-π/2", "0", "π/2", "π"]),
    azimuth = -0.4π,
    elevation = 0.15π,
    perspectiveness = 0.3
)

# サーフェスプロット
surf = surface!(ax, phi_range, psi_range, Gamma_shifted,
    colormap = :plasma
)

# カラーバー
Colorbar(fig[1, 2], surf, height = Relative(0.7))

# サブタイトル
Label(fig[2, :], "with optimizing linear stability (r = 1)", fontsize=16)

display(fig)

println("3D PCF surface plot displayed!")
println("Press Enter to close...")
readline()
