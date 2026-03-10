#!/usr/bin/env julia
"""
    visualize.jl - HONetSync 可視化スクリプト

GLMakieを使用してシミュレーション結果を可視化

使用方法:
    julia --project=. scripts/visualize.jl
"""

using HONetSync
using GLMakie
using Printf

# GLMakieのテーマ設定
set_theme!(theme_dark())

println("HONetSync Visualization")
println("="^50)

# =============================================================================
# Step 1: データ計算
# =============================================================================
println("\n[Step 1] Computing simulation data...")

N = 4
Ntheta = 101  # 可視化用に粗い解像度
params = NetworkParams(N)

# SPS計算
println("  Computing SPS...")
cmd_sps = ComputeSPS(params, Ntheta, 500.0)
result_sps = interpret(cmd_sps)
sps = unwrap(result_sps)
println("  Period T = $(round(sps.T, digits=4))")

# PSF計算
println("  Computing PSF...")
cmd_psf = ComputePSF(sps, 40)
result_psf = interpret(cmd_psf)
psf = unwrap(result_psf)

# Individual PCF計算（一部）
println("  Computing Individual PCFs...")
individual_pcfs = IndividualPCF[]
for i in 1:N, j in 1:N, k in 1:N
    cmd = ComputeIndividualPCF(sps, psf, i, j, k)
    result = interpret(cmd)
    if is_success(result)
        push!(individual_pcfs, unwrap(result))
    end
end
println("  Computed $(length(individual_pcfs)) PCFs")

# Total PCF計算（事前定義結合使用）
println("  Computing Total PCF...")
C_linear = get_predefined_coupling(:linear, N)
cmd_total = ComputeTotalPCF(individual_pcfs, C_linear, N, Ntheta)
result_total = interpret(cmd_total)
total_pcf = unwrap(result_total)

# 位相ダイナミクスシミュレーション
println("  Simulating phase dynamics...")
epsilon = 2e-4
t_end = 5000.0
cmd_phase = SimulatePhaseDynamics(total_pcf, 0.3, -0.2, epsilon, 0.0, t_end)
result_phase = interpret(cmd_phase)
pts = unwrap(result_phase)

println("  Data computation completed!")

# =============================================================================
# Step 2: Figure作成
# =============================================================================
println("\n[Step 2] Creating visualizations...")

fig = Figure(size=(1600, 1200))

# =============================================================================
# Plot 1: リミットサイクル (u-v平面)
# =============================================================================
ax1 = Axis(fig[1, 1],
    title="Limit Cycle (Oscillator 1)",
    xlabel="u",
    ylabel="v",
    aspect=DataAspect()
)

# 振動子1のリミットサイクルを抽出
Ds = 2  # 状態次元
u_cycle = sps.Xs[1, :]  # 振動子1のu成分
v_cycle = sps.Xs[2, :]  # 振動子1のv成分

lines!(ax1, u_cycle, v_cycle, linewidth=2, color=:cyan)
scatter!(ax1, [u_cycle[1]], [v_cycle[1]], markersize=15, color=:red, marker=:star5)

# 他の振動子も薄く表示
for osc in 2:N
    u_osc = sps.Xs[Ds*(osc-1)+1, :]
    v_osc = sps.Xs[Ds*(osc-1)+2, :]
    lines!(ax1, u_osc, v_osc, linewidth=1, color=(:white, 0.3))
end

# =============================================================================
# Plot 2: 位相感受関数 Q(θ)
# =============================================================================
ax2 = Axis(fig[1, 2],
    title="Phase Sensitivity Function Q(θ)",
    xlabel="θ / 2π",
    ylabel="Q"
)

theta = range(0, 1, length=Ntheta)

# 振動子1のQ成分
Q_u = psf.Q[1, :]
Q_v = psf.Q[2, :]

lines!(ax2, theta, Q_u, linewidth=2, color=:orange, label="Qᵤ")
lines!(ax2, theta, Q_v, linewidth=2, color=:green, label="Qᵥ")
axislegend(ax2, position=:rt)
hlines!(ax2, [0], color=:white, linestyle=:dash, linewidth=0.5)

# =============================================================================
# Plot 3: Total Phase Coupling Function Γ(φ, ψ)
# =============================================================================
ax3 = Axis(fig[1, 3],
    title="Total PCF Γ(φ, ψ)",
    xlabel="φ / 2π",
    ylabel="ψ / 2π",
    aspect=DataAspect()
)

phi_range = range(0, 1, length=Ntheta)
psi_range = range(0, 1, length=Ntheta)

hm = heatmap!(ax3, phi_range, psi_range, total_pcf.Gamma,
    colormap=:viridis)
Colorbar(fig[1, 4], hm, label="Γ(φ, ψ)")

# 同期点(0,0)をマーク
scatter!(ax3, [0], [0], markersize=15, color=:red, marker=:xcross)

# =============================================================================
# Plot 4: 位相ダイナミクス時系列
# =============================================================================
ax4 = Axis(fig[2, 1:2],
    title="Phase Dynamics (ε = $epsilon)",
    xlabel="Time",
    ylabel="Phase difference"
)

lines!(ax4, pts.t, pts.phi1, linewidth=1.5, color=:cyan, label="φ₁")
lines!(ax4, pts.t, pts.phi2, linewidth=1.5, color=:magenta, label="φ₂")
hlines!(ax4, [0], color=:white, linestyle=:dash, linewidth=0.5)
axislegend(ax4, position=:rt)

# =============================================================================
# Plot 5: 位相空間軌道 (φ₁-φ₂平面)
# =============================================================================
ax5 = Axis(fig[2, 3:4],
    title="Phase Space Trajectory",
    xlabel="φ₁",
    ylabel="φ₂",
    aspect=DataAspect()
)

# 軌道をプロット
lines!(ax5, pts.phi1, pts.phi2, linewidth=1, color=:yellow, alpha=0.7)

# 始点と終点
scatter!(ax5, [pts.phi1[1]], [pts.phi2[1]],
    markersize=15, color=:green, marker=:circle, label="Start")
scatter!(ax5, [pts.phi1[end]], [pts.phi2[end]],
    markersize=15, color=:red, marker=:star5, label="End")

# 同期点
scatter!(ax5, [0], [0], markersize=20, color=:white, marker=:xcross, label="Sync")
axislegend(ax5, position=:lt)

# =============================================================================
# タイトルと情報表示
# =============================================================================
Label(fig[0, :],
    "HONetSync - Higher-Order Network Synchronization (N=$N, Linear Optimized)",
    fontsize=24,
    font=:bold
)

# 安定性情報
Gamma1 = total_pcf.dGamma_dphi
Gamma2 = total_pcf.dGamma_dpsi
Lambda = compute_linear_stability(Gamma1, Gamma2)
R = compute_rotation_characteristic(Gamma1, Gamma2)

info_text = @sprintf("Γ₁ = %.4f, Γ₂ = %.4f\nΛ = %.4f, R = %.4f",
    Gamma1, Gamma2, Lambda, R)
Label(fig[3, :], info_text, fontsize=14)

# =============================================================================
# 表示
# =============================================================================
println("  Displaying figure...")
display(fig)

println("\n" * "="^50)
println("Visualization completed!")
println("="^50)
println("""
Press Enter to close the window and exit...
(Or close the window manually)
""")

# ウィンドウを開いたまま待機
readline()
