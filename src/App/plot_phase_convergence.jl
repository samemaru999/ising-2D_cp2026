#!/usr/bin/env julia
"""
    plot_phase_convergence.jl - 位相差の収束プロット

二つの位相差 φ と ψ が安定固定点に収束する様子を二次元平面にプロット

使用方法:
    julia --project=. src/App/plot_phase_convergence.jl
"""

using JLD2
using GLMakie

set_theme!(theme_dark())

# =============================================================================
# データ読み込み
# =============================================================================
println("Loading phase dynamics data...")

# Linear と Rotation の両方のデータを読み込む
data_linear = JLD2.load("data/sim/phase_dynamics_linear_N4_Nth101_20251217_211050.jld2")
data_rotation = JLD2.load("data/sim/phase_dynamics_rotation_N4_Nth101_20251217_211050.jld2")

# 位相差を抽出
phi_linear = data_linear["phase_diffs"][1, :]    # φ (位相差1)
psi_linear = data_linear["phase_diffs"][2, :]    # ψ (位相差2)
t_linear = data_linear["time"]

phi_rotation = data_rotation["phase_diffs"][1, :]
psi_rotation = data_rotation["phase_diffs"][2, :]
t_rotation = data_rotation["time"]

println("  Linear data: $(length(t_linear)) points")
println("  Rotation data: $(length(t_rotation)) points")
println("  Time range: $(t_linear[1]) to $(t_linear[end])")

# =============================================================================
# Figure作成
# =============================================================================
println("\nCreating convergence plot...")

fig = Figure(size=(1400, 700))

# =============================================================================
# Plot 1: Linear結合の位相空間軌道
# =============================================================================
ax1 = Axis(fig[1, 1],
    title="Linear Coupling - Phase Space Trajectory",
    xlabel="φ (phase difference 1)",
    ylabel="ψ (phase difference 2)",
    aspect=DataAspect()
)

# 時間に応じた色付けで軌道をプロット
# 時間を正規化して色に変換
n_points = length(t_linear)
time_normalized = (t_linear .- t_linear[1]) ./ (t_linear[end] - t_linear[1])

# 軌道を色付きの線でプロット（時間経過で色変化）
for i in 1:(n_points-1)
    # 時間に応じた色（青→黄色→赤）
    color = RGBf(time_normalized[i], 0.7 * (1 - time_normalized[i]), 1 - time_normalized[i])
    lines!(ax1, phi_linear[i:i+1], psi_linear[i:i+1],
           linewidth=1.5, color=color)
end

# 始点（緑）
scatter!(ax1, [phi_linear[1]], [psi_linear[1]],
    markersize=20, color=:lime, marker=:circle,
    strokewidth=2, strokecolor=:white)

# 終点（赤い星）
scatter!(ax1, [phi_linear[end]], [psi_linear[end]],
    markersize=20, color=:red, marker=:star5,
    strokewidth=2, strokecolor=:white)

# 安定固定点（原点）
scatter!(ax1, [0], [0],
    markersize=25, color=:white, marker=:xcross,
    strokewidth=3)

# 軌道上のいくつかの点を時刻ラベル付きでマーク
time_markers = [0.1, 0.25, 0.5, 0.75]  # 正規化時間でのマーカー位置
for tm in time_markers
    idx = round(Int, tm * (n_points - 1)) + 1
    scatter!(ax1, [phi_linear[idx]], [psi_linear[idx]],
        markersize=10, color=:yellow, marker=:diamond)
    text!(ax1, phi_linear[idx] + 0.02, psi_linear[idx] + 0.02,
        text="t=$(round(Int, t_linear[idx]))", fontsize=10, color=:yellow)
end

# =============================================================================
# Plot 2: Rotation結合の位相空間軌道
# =============================================================================
ax2 = Axis(fig[1, 2],
    title="Rotation Coupling - Phase Space Trajectory",
    xlabel="φ (phase difference 1)",
    ylabel="ψ (phase difference 2)",
    aspect=DataAspect()
)

n_points_rot = length(t_rotation)
time_normalized_rot = (t_rotation .- t_rotation[1]) ./ (t_rotation[end] - t_rotation[1])

# 軌道を色付きの線でプロット
for i in 1:(n_points_rot-1)
    color = RGBf(time_normalized_rot[i], 0.7 * (1 - time_normalized_rot[i]), 1 - time_normalized_rot[i])
    lines!(ax2, phi_rotation[i:i+1], psi_rotation[i:i+1],
           linewidth=1.5, color=color)
end

# 始点
scatter!(ax2, [phi_rotation[1]], [psi_rotation[1]],
    markersize=20, color=:lime, marker=:circle,
    strokewidth=2, strokecolor=:white)

# 終点
scatter!(ax2, [phi_rotation[end]], [psi_rotation[end]],
    markersize=20, color=:red, marker=:star5,
    strokewidth=2, strokecolor=:white)

# 安定固定点
scatter!(ax2, [0], [0],
    markersize=25, color=:white, marker=:xcross,
    strokewidth=3)

# 時刻マーカー
for tm in time_markers
    idx = round(Int, tm * (n_points_rot - 1)) + 1
    scatter!(ax2, [phi_rotation[idx]], [psi_rotation[idx]],
        markersize=10, color=:yellow, marker=:diamond)
    text!(ax2, phi_rotation[idx] + 0.02, psi_rotation[idx] + 0.02,
        text="t=$(round(Int, t_rotation[idx]))", fontsize=10, color=:yellow)
end

# =============================================================================
# カラーバー（時間スケール）
# =============================================================================
Colorbar(fig[1, 3],
    limits=(0, t_linear[end]),
    colormap=cgrad([:blue, :yellow, :red]),
    label="Time",
    height=Relative(0.8))

# =============================================================================
# タイトル
# =============================================================================
Label(fig[0, :],
    "Phase Difference Convergence to Stable Fixed Point (0, 0)",
    fontsize=22,
    font=:bold
)

# 凡例情報
legend_text = """
  ● Start point (green)  ★ End point (red)  ✕ Stable fixed point (white)
  ◆ Time markers (yellow)  Color: blue (t=0) → red (t=end)
"""
Label(fig[2, :], legend_text, fontsize=12, halign=:center)

# 収束情報
final_dist_linear = sqrt(phi_linear[end]^2 + psi_linear[end]^2)
final_dist_rotation = sqrt(phi_rotation[end]^2 + psi_rotation[end]^2)

info_text = """
Final distance from fixed point:
  Linear: $(round(final_dist_linear, digits=6))
  Rotation: $(round(final_dist_rotation, digits=6))
"""
Label(fig[3, :], info_text, fontsize=11, halign=:left)

# =============================================================================
# 表示
# =============================================================================
println("Displaying figure...")
display(fig)

println("\n" * "="^60)
println("Phase convergence plot completed!")
println("="^60)
println("""
The plot shows how the two phase differences (φ, ψ) converge
to the stable fixed point (0, 0) over time.

  - Linear coupling: Optimized for faster convergence
  - Rotation coupling: Optimized for rotating wave solutions

Press Enter to close...
""")

readline()
