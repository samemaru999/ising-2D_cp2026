#!/usr/bin/env julia
"""
    compute_n10_sps_psf.jl - N=10 Tree vs Full Coupling Comparison

N=10振動子ネットワークにおけるTree型結合と全結合での
安定周期解（SPS）と位相感受関数（PSF）の計算・比較

使用方法:
    julia --project=. scripts/compute_n10_sps_psf.jl
"""

using HONetSync
using GLMakie
using Printf
using JLD2

println("="^60)
println("N=10 Network: Tree vs Full Coupling Comparison")
println("="^60)

# =============================================================================
# 共通パラメータ
# =============================================================================
const N = 10
const Ntheta = 101  # 位相分割数
const T_RELAX = 500.0  # 緩和時間
const N_ITER_PSF = 100  # PSF計算の反復回数（収束のため増加）

# =============================================================================
# Tree型結合でのSPS・PSF計算
# =============================================================================
println("\n[1/4] Computing SPS with TREE coupling...")

# Tree型結合行列を使用してNetworkParamsを作成
fhn = FHNParams()

# Tree結合: 全ノードがペースメーカー（I=0.8）
I_ext_tree = fill(I_PACEMAKER, N)
println("  Tree I settings: all nodes I = $(I_PACEMAKER)")

K_tree = get_intra_K(N; topology=:tree)

params_tree = NetworkParams(N, fhn, I_ext_tree, K_tree, TreeCoupling([0, 1, 1, 1, 2, 2, 3, 3, 4, 4]))

# SPS計算
cmd_sps_tree = ComputeSPS(params_tree, Ntheta, T_RELAX)
result_sps_tree = interpret(cmd_sps_tree)

if !is_success(result_sps_tree)
    error("Tree SPS computation failed: $(result_sps_tree.error)")
end

sps_tree = unwrap(result_sps_tree)
println("  -> Period T_tree = $(@sprintf("%.4f", sps_tree.T)) [s]")
println("  -> Angular frequency omega = $(@sprintf("%.4f", sps_tree.omega)) [rad/s]")

# PSF計算（バージョン3: 細かい時間刻み dt=0.001）
println("\n[2/4] Computing PSF with TREE coupling (version=3)...")
psf_tree = compute_psf(sps_tree, params_tree; version=3, n_iterations=N_ITER_PSF)
println("  -> PSF computed: Q matrix size = $(size(psf_tree.Q))")

# =============================================================================
# 全結合でのSPS・PSF計算
# =============================================================================
println("\n[3/4] Computing SPS with FULL coupling...")

# Full結合: Nakao論文設定（I[8,9,10]=0.8, 他は0.2）
I_ext_full = get_default_external_input(N)
println("  Full I settings: I = $I_ext_full")

K_full = get_intra_K(N; topology=:full)
params_full = NetworkParams(N, fhn, I_ext_full, K_full, FullCoupling())

# SPS計算
cmd_sps_full = ComputeSPS(params_full, Ntheta, T_RELAX)
result_sps_full = interpret(cmd_sps_full)

if !is_success(result_sps_full)
    error("Full SPS computation failed: $(result_sps_full.error)")
end

sps_full = unwrap(result_sps_full)
println("  -> Period T_full = $(@sprintf("%.4f", sps_full.T)) [s]")
println("  -> Angular frequency omega = $(@sprintf("%.4f", sps_full.omega)) [rad/s]")

# PSF計算（バージョン3: 細かい時間刻み dt=0.001）
println("\n[4/4] Computing PSF with FULL coupling (version=3)...")
psf_full = compute_psf(sps_full, params_full; version=3, n_iterations=N_ITER_PSF)
println("  -> PSF computed: Q matrix size = $(size(psf_full.Q))")

# =============================================================================
# 結果の保存
# =============================================================================
output_dir = "output/n10_comparison"
mkpath(output_dir)

# JLD2で保存
save_file = joinpath(output_dir, "n10_tree_full_comparison.jld2")
JLD2.jldsave(save_file;
    sps_tree = sps_tree,
    psf_tree = psf_tree,
    sps_full = sps_full,
    psf_full = psf_full,
    K_tree = K_tree,
    K_full = K_full
)
println("\nData saved to: $save_file")

# =============================================================================
# プロット
# =============================================================================
println("\n" * "="^60)
println("Generating Plots...")
println("="^60)

# 位相軸
theta = range(0, 2π, length=Ntheta)

# カラーマップ
colors_tree = [:steelblue, :royalblue, :dodgerblue, :deepskyblue, :lightskyblue,
               :cornflowerblue, :mediumslateblue, :slateblue, :darkslateblue, :midnightblue]
colors_full = [:crimson, :firebrick, :indianred, :salmon, :lightsalmon,
               :coral, :tomato, :orangered, :red, :darkred]

# --- Figure 1: 安定周期解 (SPS) の比較 ---
fig1 = Figure(size=(1600, 900), backgroundcolor=:white)

# Tree coupling SPS
ax1 = Axis(fig1[1, 1],
    title = "Stable Periodic Solution (Tree Coupling)",
    xlabel = "Phase θ",
    ylabel = "v_i(θ)",
    xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]),
    backgroundcolor = :white
)

for i in 1:N
    v_idx = 2 * i
    lines!(ax1, theta, sps_tree.Xs[v_idx, :], linewidth=2,
           color=colors_tree[i], label="v_$i")
end
axislegend(ax1, position=:rt, framevisible=false, labelsize=10)

# Full coupling SPS
ax2 = Axis(fig1[1, 2],
    title = "Stable Periodic Solution (Full Coupling)",
    xlabel = "Phase θ",
    ylabel = "v_i(θ)",
    xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]),
    backgroundcolor = :white
)

for i in 1:N
    v_idx = 2 * i
    lines!(ax2, theta, sps_full.Xs[v_idx, :], linewidth=2,
           color=colors_full[i], label="v_$i")
end
axislegend(ax2, position=:rt, framevisible=false, labelsize=10)

# Tree/Full 比較 (v_1のみ)
ax3 = Axis(fig1[2, 1:2],
    title = "SPS Comparison: v_1 (Tree vs Full)",
    xlabel = "Phase θ",
    ylabel = "v_1(θ)",
    xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]),
    backgroundcolor = :white
)

lines!(ax3, theta, sps_tree.Xs[2, :], linewidth=3, color=:steelblue,
       label="Tree (T=$(@sprintf("%.2f", sps_tree.T)))")
lines!(ax3, theta, sps_full.Xs[2, :], linewidth=3, color=:crimson, linestyle=:dash,
       label="Full (T=$(@sprintf("%.2f", sps_full.T)))")
axislegend(ax3, position=:rt, framevisible=false)

save(joinpath(output_dir, "sps_comparison_n10.png"), fig1, px_per_unit=2)
println("Saved: sps_comparison_n10.png")

# --- Figure 2: 位相感受関数 (PSF) の比較 ---
fig2 = Figure(size=(1600, 900), backgroundcolor=:white)

# Tree coupling PSF
ax4 = Axis(fig2[1, 1],
    title = "Phase Sensitivity Function (Tree Coupling)",
    xlabel = "Phase θ",
    ylabel = "Q_v_i(θ)",
    xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]),
    backgroundcolor = :white
)

for i in 1:N
    v_idx = 2 * i
    lines!(ax4, theta, psf_tree.Q[v_idx, :], linewidth=2,
           color=colors_tree[i], label="Q_v$i")
end
axislegend(ax4, position=:rt, framevisible=false, labelsize=10)

# Full coupling PSF
ax5 = Axis(fig2[1, 2],
    title = "Phase Sensitivity Function (Full Coupling)",
    xlabel = "Phase θ",
    ylabel = "Q_v_i(θ)",
    xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]),
    backgroundcolor = :white
)

for i in 1:N
    v_idx = 2 * i
    lines!(ax5, theta, psf_full.Q[v_idx, :], linewidth=2,
           color=colors_full[i], label="Q_v$i")
end
axislegend(ax5, position=:rt, framevisible=false, labelsize=10)

# Tree/Full 比較 (Q_v_1のみ)
ax6 = Axis(fig2[2, 1:2],
    title = "PSF Comparison: Q_v_1 (Tree vs Full)",
    xlabel = "Phase θ",
    ylabel = "Q_v_1(θ)",
    xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]),
    backgroundcolor = :white
)

lines!(ax6, theta, psf_tree.Q[2, :], linewidth=3, color=:steelblue, label="Tree")
lines!(ax6, theta, psf_full.Q[2, :], linewidth=3, color=:crimson, linestyle=:dash, label="Full")
axislegend(ax6, position=:rt, framevisible=false)

save(joinpath(output_dir, "psf_comparison_n10.png"), fig2, px_per_unit=2)
println("Saved: psf_comparison_n10.png")

# --- Figure 3: PSF Tree Coupling (5x2グリッド) ---
fig3_tree = Figure(size=(1400, 1600), backgroundcolor=:white)

Label(fig3_tree[0, 1:2], "Phase Sensitivity Functions (Tree Coupling)",
      fontsize=20, font=:bold)

for i in 1:N
    row = (i - 1) ÷ 2 + 1
    col = (i - 1) % 2 + 1

    ax = Axis(fig3_tree[row, col],
        title = "Oscillator $i",
        xlabel = row == 5 ? "Phase θ" : "",
        ylabel = "Q_v_$i(θ)",
        xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]),
        backgroundcolor = :white
    )

    v_idx = 2 * i
    lines!(ax, theta, psf_tree.Q[v_idx, :], linewidth=2, color=colors_tree[i])
end

save(joinpath(output_dir, "psf_tree_n10.png"), fig3_tree, px_per_unit=2)
println("Saved: psf_tree_n10.png")

# --- Figure: PSF Tree Coupling (1つの図に全振動子を重ねる) ---
fig_tree_combined = Figure(size=(800, 600), backgroundcolor=:white)

ax_tree_combined = Axis(fig_tree_combined[1, 1],
    title = "Phase Sensitivity Functions (Tree Coupling)",
    xlabel = "θ",
    ylabel = "Qᵥⁱ",
    xticks = ([0, π, 2π], ["0", "π", "2π"]),
    backgroundcolor = :white
)

# 水平線 (y=0)
hlines!(ax_tree_combined, [0], color=:black, linestyle=:dot, linewidth=1)

# #1 (ルートノード) - 赤
lines!(ax_tree_combined, theta, psf_tree.Q[2, :], linewidth=2, color=:red, label="#1 (root)")

# #2-#4 (レベル1ノード) - 青
for i in 2:4
    v_idx = 2 * i
    lines!(ax_tree_combined, theta, psf_tree.Q[v_idx, :], linewidth=2, color=:blue,
           label = i == 2 ? "#2-#4 (level 1)" : nothing)
end

# #5-#10 (リーフノード) - 緑
for i in 5:10
    v_idx = 2 * i
    lines!(ax_tree_combined, theta, psf_tree.Q[v_idx, :], linewidth=2, color=:green,
           label = i == 5 ? "#5-#10 (leaves)" : nothing)
end

axislegend(ax_tree_combined, position=:rt, framevisible=false)

save(joinpath(output_dir, "psf_tree_combined_n10.png"), fig_tree_combined, px_per_unit=2)
println("Saved: psf_tree_combined_n10.png")

# --- Figure 4: PSF Full Coupling (5x2グリッド) ---
fig3_full = Figure(size=(1400, 1600), backgroundcolor=:white)

Label(fig3_full[0, 1:2], "Phase Sensitivity Functions (Full Coupling)",
      fontsize=20, font=:bold)

for i in 1:N
    row = (i - 1) ÷ 2 + 1
    col = (i - 1) % 2 + 1

    ax = Axis(fig3_full[row, col],
        title = "Oscillator $i",
        xlabel = row == 5 ? "Phase θ" : "",
        ylabel = "Q_v_$i(θ)",
        xticks = ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]),
        backgroundcolor = :white
    )

    v_idx = 2 * i
    lines!(ax, theta, psf_full.Q[v_idx, :], linewidth=2, color=colors_full[i])
end

save(joinpath(output_dir, "psf_full_n10.png"), fig3_full, px_per_unit=2)
println("Saved: psf_full_n10.png")

# --- Figure 4: 結合行列のヒートマップ ---
fig4 = Figure(size=(1200, 500), backgroundcolor=:white)

ax7 = Axis(fig4[1, 1],
    title = "Tree Coupling Matrix",
    xlabel = "j",
    ylabel = "i",
    aspect = 1
)
hm1 = heatmap!(ax7, K_tree, colormap=:RdBu)
Colorbar(fig4[1, 2], hm1, label="K_ij")

ax8 = Axis(fig4[1, 3],
    title = "Full Coupling Matrix",
    xlabel = "j",
    ylabel = "i",
    aspect = 1
)
hm2 = heatmap!(ax8, K_full, colormap=:RdBu)
Colorbar(fig4[1, 4], hm2, label="K_ij")

save(joinpath(output_dir, "coupling_matrices_n10.png"), fig4, px_per_unit=2)
println("Saved: coupling_matrices_n10.png")

# =============================================================================
# 結果サマリー
# =============================================================================
println("\n" * "="^60)
println("COMPUTATION COMPLETED")
println("="^60)
println("""
Summary:
  Tree Coupling:
    - Period T = $(@sprintf("%.4f", sps_tree.T)) [s]
    - Omega = $(@sprintf("%.4f", sps_tree.omega)) [rad/s]

  Full Coupling:
    - Period T = $(@sprintf("%.4f", sps_full.T)) [s]
    - Omega = $(@sprintf("%.4f", sps_full.omega)) [rad/s]

Output files saved to: $output_dir
  - n10_tree_full_comparison.jld2
  - sps_comparison_n10.png
  - psf_comparison_n10.png
  - psf_tree_n10.png
  - psf_tree_combined_n10.png
  - psf_full_n10.png
  - coupling_matrices_n10.png
""")

# プロットを表示
display(fig1)
display(fig2)
display(fig3_tree)
display(fig_tree_combined)
display(fig3_full)
display(fig4)

println("Press Enter to close plots...")
readline()
