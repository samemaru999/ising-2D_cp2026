#!/usr/bin/env julia
"""
    plot_sps_psf_per_oscillator.jl - 振動子ごとのSPS/PSFプロット

安定周期解と位相感受関数を計算し、振動子ごとに図を作成

使用方法:
    julia --project=. scripts/plot_sps_psf_per_oscillator.jl
"""

using HONetSync
using GLMakie
using Printf

println("="^60)
println("SPS/PSF Computation and Per-Oscillator Plotting")
println("="^60)

# =============================================================================
# パラメータ設定
# =============================================================================
const N = 10          # 振動子数
const Ntheta = 101    # 位相分割数
const T_RELAX = 500.0 # 緩和時間
const N_ITER_PSF = 50 # PSF計算の反復回数

# =============================================================================
# ネットワークパラメータ作成
# =============================================================================
println("\n[1/3] Setting up network parameters...")

params = NetworkParams(N)
println("  N = $N oscillators")
println("  δ = $(params.fhn.delta), a = $(params.fhn.a), b = $(params.fhn.b)")
println("  I = $(params.I)")

# =============================================================================
# SPS計算
# =============================================================================
println("\n[2/3] Computing Stable Periodic Solution...")

cmd_sps = ComputeSPS(params, Ntheta, T_RELAX)
result_sps = interpret(cmd_sps)

if !is_success(result_sps)
    error("SPS computation failed: $(result_sps.error)")
end

sps = unwrap(result_sps)
println("  -> Period T = $(@sprintf("%.4f", sps.T)) [s]")
println("  -> Angular frequency ω = $(@sprintf("%.4f", sps.omega)) [rad/s]")

# =============================================================================
# PSF計算
# =============================================================================
println("\n[3/3] Computing Phase Sensitivity Function...")

psf = compute_psf(sps, params; version=3, n_iterations=N_ITER_PSF)
println("  -> PSF computed: Q matrix size = $(size(psf.Q))")

# =============================================================================
# プロット作成
# =============================================================================
println("\n" * "="^60)
println("Generating Plots...")
println("="^60)

# 出力ディレクトリ
output_dir = "output/per_oscillator"
mkpath(output_dir)

# 位相/時間軸
theta = range(0, 2π, length=Ntheta)
t_axis = range(0, sps.T, length=Ntheta)

# カラー
colors = [:blue, :orange, :green, :purple, :red, :cyan, :magenta, :brown, :pink, :gray]

# レイアウト設定
n_cols = 5
n_rows = ceil(Int, N / n_cols)

# --- Figure 1: 安定周期解 u(t), v(t) 振動子ごと ---
println("Creating SPS plot (u(t), v(t) per oscillator)...")

fig_sps = Figure(size=(1600, 400 * n_rows), backgroundcolor=:white)

Label(fig_sps[0, 1:n_cols], "Stable Periodic Solution: u(t), v(t) per Oscillator (N=$N, T=$(@sprintf("%.2f", sps.T)))",
      fontsize=20)

for i in 1:N
    row_u = 2 * (i - 1) + 1
    row_v = 2 * i

    grid_row = div(i - 1, n_cols) + 1
    grid_col = mod(i - 1, n_cols) + 1

    ax = Axis(fig_sps[grid_row, grid_col],
        title = "Oscillator $i",
        xlabel = grid_row == n_rows ? "t" : "",
        ylabel = "u, v",
        backgroundcolor = :white
    )

    lines!(ax, t_axis, sps.Xs[row_u, :], linewidth=2, color=:blue, label="u")
    lines!(ax, t_axis, sps.Xs[row_v, :], linewidth=2, color=:red, linestyle=:dash, label="v")

    if i == N
        axislegend(ax, position=:rt, framevisible=false)
    end
end

save(joinpath(output_dir, "sps_per_oscillator_N$(N).png"), fig_sps)
println("  Saved: sps_per_oscillator_N$(N).png")

# --- Figure 2: 位相感受関数 Q_u(θ), Q_v(θ) 振動子ごと ---
println("Creating PSF plot (Q_u, Q_v per oscillator)...")

fig_psf = Figure(size=(1600, 400 * n_rows), backgroundcolor=:white)

Label(fig_psf[0, 1:n_cols], "Phase Sensitivity Function: Q_u(θ), Q_v(θ) per Oscillator (N=$N)",
      fontsize=20)

for i in 1:N
    row_u = 2 * (i - 1) + 1
    row_v = 2 * i

    grid_row = div(i - 1, n_cols) + 1
    grid_col = mod(i - 1, n_cols) + 1

    ax = Axis(fig_psf[grid_row, grid_col],
        title = "Oscillator $i",
        xlabel = grid_row == n_rows ? "θ" : "",
        ylabel = "Q",
        xticks = ([0, π, 2π], ["0", "π", "2π"]),
        backgroundcolor = :white
    )

    # 水平線 y=0
    hlines!(ax, [0], color=:gray, linestyle=:dot, linewidth=1)

    lines!(ax, theta, psf.Q[row_u, :], linewidth=2, color=:blue, label="Q_u")
    lines!(ax, theta, psf.Q[row_v, :], linewidth=2, color=:red, linestyle=:dash, label="Q_v")

    if i == N
        axislegend(ax, position=:rt, framevisible=false)
    end
end

save(joinpath(output_dir, "psf_per_oscillator_N$(N).png"), fig_psf)
println("  Saved: psf_per_oscillator_N$(N).png")

# --- Figure 3: リミットサイクル u-v平面 振動子ごと ---
println("Creating limit cycle plot (u-v plane per oscillator)...")

fig_lc = Figure(size=(1600, 400 * n_rows), backgroundcolor=:white)

Label(fig_lc[0, 1:n_cols], "Limit Cycles in u-v Plane (N=$N)",
      fontsize=20)

for i in 1:N
    row_u = 2 * (i - 1) + 1
    row_v = 2 * i

    grid_row = div(i - 1, n_cols) + 1
    grid_col = mod(i - 1, n_cols) + 1

    ax = Axis(fig_lc[grid_row, grid_col],
        title = "Oscillator $i",
        xlabel = "u",
        ylabel = "v",
        aspect = DataAspect(),
        backgroundcolor = :white
    )

    lines!(ax, sps.Xs[row_u, :], sps.Xs[row_v, :], linewidth=2, color=colors[i])
end

save(joinpath(output_dir, "limit_cycles_N$(N).png"), fig_lc)
println("  Saved: limit_cycles_N$(N).png")

# =============================================================================
# 結果サマリー
# =============================================================================
println("\n" * "="^60)
println("COMPLETED")
println("="^60)
println("""
Summary:
  N = $N oscillators
  Period T = $(@sprintf("%.4f", sps.T)) [s]
  Omega = $(@sprintf("%.4f", sps.omega)) [rad/s]

Output files saved to: $output_dir
  - sps_per_oscillator_N$(N).png
  - psf_per_oscillator_N$(N).png
  - limit_cycles_N$(N).png
""")

println("All plots saved to: $output_dir")
