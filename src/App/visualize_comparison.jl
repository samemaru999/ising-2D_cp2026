#!/usr/bin/env julia
"""
    visualize_comparison.jl - 結合戦略比較の可視化

3つの結合戦略（Uniform, Linear, Rotation）の比較可視化

使用方法:
    julia --project=. scripts/visualize_comparison.jl
"""

using HONetSync
using GLMakie
using Printf

set_theme!(theme_dark())

println("HONetSync Strategy Comparison Visualization")
println("="^60)

# =============================================================================
# データ計算
# =============================================================================
println("\n[Step 1] Computing data for comparison...")

N = 4
Ntheta = 51
params = NetworkParams(N)

# SPS & PSF計算
println("  Computing SPS...")
sps = unwrap(interpret(ComputeSPS(params, Ntheta, 100.0)))

println("  Computing PSF...")
psf = unwrap(interpret(ComputePSF(sps, 5)))

# Individual PCF計算
println("  Computing Individual PCFs...")
individual_pcfs = IndividualPCF[]
for i in 1:N, j in 1:N, k in 1:N
    result = interpret(ComputeIndividualPCF(sps, psf, i, j, k))
    is_success(result) && push!(individual_pcfs, unwrap(result))
end

# 3つの戦略で計算
strategies = [:uniform, :linear, :rotation]
strategy_names = ["Uniform", "Linear Optimized", "Rotation Optimized"]
colors = [:gray, :cyan, :magenta]

results = Dict{Symbol, NamedTuple}()

println("  Computing Total PCFs for each strategy...")
for (strategy, name) in zip(strategies, strategy_names)
    C = get_predefined_coupling(strategy, N)
    total_pcf = unwrap(interpret(ComputeTotalPCF(individual_pcfs, C, N, Ntheta)))

    Gamma1 = total_pcf.dGamma_dphi
    Gamma2 = total_pcf.dGamma_dpsi
    Lambda = compute_linear_stability(Gamma1, Gamma2)
    R = compute_rotation_characteristic(Gamma1, Gamma2)

    results[strategy] = (
        C = C,
        total_pcf = total_pcf,
        Gamma1 = Gamma1,
        Gamma2 = Gamma2,
        Lambda = Lambda,
        R = R,
        name = name
    )

    @printf("  %s: Λ = %.4f, R = %.4f\n", name, Lambda, R)
end

# 位相ダイナミクスシミュレーション
println("  Simulating phase dynamics...")
epsilon = 2e-4
t_end = 8000.0
phi1_0, phi2_0 = 0.3, -0.2

phase_results = Dict{Symbol, PhaseTimeSeries}()
for strategy in strategies
    total_pcf = results[strategy].total_pcf
    cmd = SimulatePhaseDynamics(total_pcf, phi1_0, phi2_0, epsilon, 0.0, t_end)
    phase_results[strategy] = unwrap(interpret(cmd))
end

# =============================================================================
# 可視化
# =============================================================================
println("\n[Step 2] Creating comparison figures...")

fig = Figure(size=(1800, 1200))

# --- Row 1: Phase Coupling Functions ---
for (idx, strategy) in enumerate(strategies)
    ax = Axis(fig[1, idx],
        title = "Γ(φ,ψ) - $(results[strategy].name)",
        xlabel = "φ / 2π",
        ylabel = "ψ / 2π",
        aspect = DataAspect()
    )

    phi_range = range(0, 1, length=Ntheta)
    psi_range = range(0, 1, length=Ntheta)

    hm = heatmap!(ax, phi_range, psi_range, results[strategy].total_pcf.Gamma,
        colormap = :viridis)

    scatter!(ax, [0], [0], markersize=12, color=:red, marker=:xcross)

    if idx == 3
        Colorbar(fig[1, 4], hm, label="Γ")
    end
end

# --- Row 2: Phase Dynamics Comparison ---
ax_phi1 = Axis(fig[2, 1:2],
    title = "Phase Difference φ₁(t)",
    xlabel = "Time",
    ylabel = "φ₁"
)

ax_phi2 = Axis(fig[2, 3:4],
    title = "Phase Difference φ₂(t)",
    xlabel = "Time",
    ylabel = "φ₂"
)

for (strategy, color) in zip(strategies, colors)
    pts = phase_results[strategy]
    name = results[strategy].name
    lines!(ax_phi1, pts.t, pts.phi1, linewidth=1.5, color=color, label=name)
    lines!(ax_phi2, pts.t, pts.phi2, linewidth=1.5, color=color, label=name)
end

hlines!(ax_phi1, [0], color=:white, linestyle=:dash, linewidth=0.5)
hlines!(ax_phi2, [0], color=:white, linestyle=:dash, linewidth=0.5)
axislegend(ax_phi1, position=:rt)
axislegend(ax_phi2, position=:rt)

# --- Row 3: Phase Space Trajectories ---
for (idx, (strategy, color)) in enumerate(zip(strategies, colors))
    ax = Axis(fig[3, idx],
        title = "Phase Space - $(results[strategy].name)",
        xlabel = "φ₁",
        ylabel = "φ₂",
        aspect = DataAspect()
    )

    pts = phase_results[strategy]

    lines!(ax, pts.phi1, pts.phi2, linewidth=1, color=color, alpha=0.8)
    scatter!(ax, [pts.phi1[1]], [pts.phi2[1]],
        markersize=12, color=:green, marker=:circle)
    scatter!(ax, [pts.phi1[end]], [pts.phi2[end]],
        markersize=12, color=:red, marker=:star5)
    scatter!(ax, [0], [0], markersize=15, color=:white, marker=:xcross)
end

# --- Row 3, Col 4: Stability Comparison Bar Chart ---
ax_bar = Axis(fig[3, 4],
    title = "Stability Comparison",
    xlabel = "Strategy",
    ylabel = "Value",
    xticks = (1:3, strategy_names)
)

Lambda_values = [results[s].Lambda for s in strategies]
R_values = [results[s].R for s in strategies]

barplot!(ax_bar, [0.8, 1.8, 2.8], Lambda_values,
    color = :cyan, width=0.35, label="Λ (stability)")
barplot!(ax_bar, [1.2, 2.2, 3.2], R_values,
    color = :magenta, width=0.35, label="R (rotation)")

hlines!(ax_bar, [0], color=:white, linestyle=:dash, linewidth=0.5)
axislegend(ax_bar, position=:rt)

# --- Title ---
Label(fig[0, :],
    "HONetSync - Coupling Strategy Comparison (N=$N, ε=$epsilon)",
    fontsize=24,
    font=:bold
)

# --- Info Table ---
info_text = ""
for strategy in strategies
    r = results[strategy]
    info_text *= @sprintf("%s: Γ₁=%.3f, Γ₂=%.3f, Λ=%.4f, R=%.4f\n",
        r.name, r.Gamma1, r.Gamma2, r.Lambda, r.R)
end
Label(fig[4, :], info_text, fontsize=12, halign=:left)

# =============================================================================
# 表示
# =============================================================================
println("  Displaying figure...")
display(fig)

println("\n" * "="^60)
println("Comparison visualization completed!")
println("="^60)
println("""
Key observations:
  - Linear optimized: Maximizes |Λ| for faster convergence to sync
  - Rotation optimized: Maximizes R for rotating wave solutions
  - Uniform: Baseline with equal coupling strengths

Press Enter to close...
""")

readline()
