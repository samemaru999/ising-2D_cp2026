#!/usr/bin/env julia
"""
    test_periodic_orbit.jl - 周期軌道の検証スクリプト

周期時間分の実際のシミュレーションを行い、軌道を確認する。
Nθによるサンプリングではなく、dt=0.001の細かい時間刻みで積分。
"""

using HONetSync
using GLMakie
using Printf
using LinearAlgebra
using Statistics

println("="^60)
println("Periodic Orbit Verification Test")
println("="^60)

# =============================================================================
# パラメータ設定
# =============================================================================
const N = 10
const T_RELAX = 500.0
const DT = 0.001  # 時間刻み

# =============================================================================
# ネットワークパラメータ作成
# =============================================================================
println("\n[1/4] Setting up network parameters...")

params = NetworkParams(N)
println("  N = $N oscillators")
println("  δ = $(params.fhn.delta), a = $(params.fhn.a), b = $(params.fhn.b)")
println("  I = $(params.I)")

# =============================================================================
# 緩和シミュレーション
# =============================================================================
println("\n[2/4] Running relaxation simulation...")

# 初期条件
X0 = zeros(Ds * N)
for i in 1:N
    X0[2i-1] = 0.5 + 0.1 * randn()  # u
    X0[2i] = 0.0 + 0.1 * randn()    # v
end

# 緩和シミュレーション
using DifferentialEquations
function dynamics!(dX, X, p, t)
    dX .= network_full_dynamics(X, params)
end

prob_relax = ODEProblem(dynamics!, X0, (0.0, T_RELAX))
sol_relax = solve(prob_relax, Tsit5(); abstol=1e-8, reltol=1e-8, saveat=0.1)

X_relaxed = sol_relax.u[end]
println("  Relaxation completed. Final state norm = $(norm(X_relaxed))")

# =============================================================================
# 周期推定（ゼロクロス法）
# =============================================================================
println("\n[3/4] Estimating period via zero-crossing...")

# 追加シミュレーションで周期を計測
t_search = 200.0
prob_search = ODEProblem(dynamics!, X_relaxed, (0.0, t_search))
sol_search = solve(prob_search, Tsit5(); abstol=1e-8, reltol=1e-8, saveat=DT)

v1_series = [u[2] for u in sol_search.u]
t_series = sol_search.t

# ゼロクロス検出
origin_threshold = 0.0
crossings = Int[]
for i in 2:length(v1_series)
    if v1_series[i-1] < origin_threshold && v1_series[i] > origin_threshold
        push!(crossings, i)
    end
end

if length(crossings) >= 2
    periods = [t_series[crossings[i]] - t_series[crossings[i-1]] for i in 2:length(crossings)]
    T_estimated = mean(periods)
    println("  Detected $(length(crossings)) crossings")
    println("  Period T = $(@sprintf("%.4f", T_estimated)) [s]")
    println("  Period range: $(@sprintf("%.4f", minimum(periods))) - $(@sprintf("%.4f", maximum(periods)))")
else
    error("Could not detect enough zero-crossings")
end

# 位相原点からスタート（最初のゼロクロス点）
X_origin = sol_search.u[crossings[1]]

# =============================================================================
# 1周期分のシミュレーション（細かい時間刻み）
# =============================================================================
println("\n[4/4] Simulating one period with fine time steps...")

t_orbit, X_orbit = simulate_one_period(X_origin, params, T_estimated; dt=DT)
println("  Simulated $(length(t_orbit)) time points over T = $(@sprintf("%.4f", T_estimated)) [s]")

# 周期性の確認
endpoint_error = norm(X_orbit[1] - X_orbit[end])
println("  Endpoint error = $(@sprintf("%.6e", endpoint_error))")

# =============================================================================
# プロット
# =============================================================================
println("\n" * "="^60)
println("Generating Plots...")
println("="^60)

output_dir = "output/periodic_orbit_test"
mkpath(output_dir)

# --- Figure 1: 全振動子のv(t)を重ねてプロット ---
fig1 = Figure(size=(1200, 600), backgroundcolor=:white)

ax1 = Axis(fig1[1, 1],
    title = "Periodic Orbit: v_i(t) for all oscillators (N=$N)",
    xlabel = "t [s]",
    ylabel = "v",
    backgroundcolor = :white
)

colors = [:blue, :orange, :green, :red, :purple, :cyan, :magenta, :brown, :pink, :gray]
for i in 1:N
    v_i = [X[2i] for X in X_orbit]
    lines!(ax1, t_orbit, v_i, linewidth=1.5, color=colors[i], label="v_$i")
end

axislegend(ax1, position=:rt, framevisible=false, labelsize=10)

save(joinpath(output_dir, "periodic_orbit_v_all.png"), fig1)
println("Saved: periodic_orbit_v_all.png")

# --- Figure 2: 振動子ごとのu(t), v(t) (5x2グリッド) ---
fig2 = Figure(size=(1400, 1200), backgroundcolor=:white)

Label(fig2[0, 1:2], "Periodic Orbit: u(t), v(t) per Oscillator (T=$(@sprintf("%.2f", T_estimated)) s)",
      fontsize=18)

for i in 1:N
    row = (i - 1) ÷ 2 + 1
    col = (i - 1) % 2 + 1

    ax = Axis(fig2[row, col],
        title = "Oscillator $i",
        xlabel = row == 5 ? "t [s]" : "",
        ylabel = "u, v",
        backgroundcolor = :white
    )

    u_i = [X[2i-1] for X in X_orbit]
    v_i = [X[2i] for X in X_orbit]

    lines!(ax, t_orbit, u_i, linewidth=1.5, color=:blue, label="u")
    lines!(ax, t_orbit, v_i, linewidth=1.5, color=:red, linestyle=:dash, label="v")

    if i == N
        axislegend(ax, position=:rt, framevisible=false)
    end
end

save(joinpath(output_dir, "periodic_orbit_uv_per_osc.png"), fig2)
println("Saved: periodic_orbit_uv_per_osc.png")

# --- Figure 3: u-v位相空間 (5x2グリッド) ---
fig3 = Figure(size=(1400, 1200), backgroundcolor=:white)

Label(fig3[0, 1:2], "Limit Cycles in u-v Plane (N=$N)",
      fontsize=18)

for i in 1:N
    row = (i - 1) ÷ 2 + 1
    col = (i - 1) % 2 + 1

    ax = Axis(fig3[row, col],
        title = "Oscillator $i",
        xlabel = "u",
        ylabel = "v",
        aspect = DataAspect(),
        backgroundcolor = :white
    )

    u_i = [X[2i-1] for X in X_orbit]
    v_i = [X[2i] for X in X_orbit]

    lines!(ax, u_i, v_i, linewidth=2, color=colors[i])

    # 開始点をマーク
    scatter!(ax, [u_i[1]], [v_i[1]], markersize=10, color=:black)
end

save(joinpath(output_dir, "limit_cycles_uv.png"), fig3)
println("Saved: limit_cycles_uv.png")

# --- Figure 4: v_1の詳細（位相原点の確認） ---
fig4 = Figure(size=(800, 400), backgroundcolor=:white)

ax4 = Axis(fig4[1, 1],
    title = "v_1(t) with zero-crossing (phase origin)",
    xlabel = "t [s]",
    ylabel = "v_1",
    backgroundcolor = :white
)

v1 = [X[2] for X in X_orbit]
lines!(ax4, t_orbit, v1, linewidth=2, color=:blue)
hlines!(ax4, [0], color=:gray, linestyle=:dash, linewidth=1)
scatter!(ax4, [t_orbit[1]], [v1[1]], markersize=12, color=:red, marker=:circle)
scatter!(ax4, [t_orbit[end]], [v1[end]], markersize=12, color=:green, marker=:star5)

save(joinpath(output_dir, "v1_phase_origin.png"), fig4)
println("Saved: v1_phase_origin.png")

# =============================================================================
# 結果サマリー
# =============================================================================
println("\n" * "="^60)
println("COMPLETED")
println("="^60)
println("""
Summary:
  N = $N oscillators
  Period T = $(@sprintf("%.4f", T_estimated)) [s]
  Time step dt = $DT [s]
  Number of time points = $(length(t_orbit))
  Endpoint error = $(@sprintf("%.6e", endpoint_error))

Output files saved to: $output_dir
  - periodic_orbit_v_all.png
  - periodic_orbit_uv_per_osc.png
  - limit_cycles_uv.png
  - v1_phase_origin.png
""")
