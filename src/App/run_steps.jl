#!/usr/bin/env julia
"""
    run_steps.jl - Step-by-step execution of HONetSync pipeline
    
    Ref: README.md "REPL での対話的実行"
"""

using HONetSync
using Printf
using LinearAlgebra
using Random

# Set random seed for reproducibility
Random.seed!(1234)

# Output directory
OUTPUT_DIR = "output_steps"
if !isdir(OUTPUT_DIR)
    mkdir(OUTPUT_DIR)
end

println("HONetSync Step-by-Step Execution")
println("="^60)

# =============================================================================
# Step 1: Network Parameters
# =============================================================================
println("\n[Step 1] Setting up Network Parameters")
N = 4
params = NetworkParams(N)
println("  N = $N")
println("  FHN Parameters: δ=$(params.fhn.delta), a=$(params.fhn.a), b=$(params.fhn.b)")

# =============================================================================
# Step 2: Stable Periodic Solution (SPS)
# =============================================================================
println("\n[Step 2] Computing Stable Periodic Solution (SPS)")
cmd_sps = ComputeSPS(params, 101, 2000.0)
@time result_sps = interpret(cmd_sps)

if !is_success(result_sps)
    error("SPS failed: $(result_sps.error)")
end
sps = unwrap(result_sps)
println("  Period T = $(sps.T)")
println("  Omega ω = $(sps.omega)")
save_sps(joinpath(OUTPUT_DIR, "sps.jld2"), sps)

# =============================================================================
# Step 3: Phase Sensitivity Function (PSF)
# =============================================================================
println("\n[Step 3] Computing Phase Sensitivity Function (PSF)")
# Using DiffusiveCoupling as in README example
cmd_psf = ComputePSF(sps, 40, DiffusiveCoupling())
@time result_psf = interpret(cmd_psf)

if !is_success(result_psf)
    error("PSF failed: $(result_psf.error)")
end
psf = unwrap(result_psf)
println("  PSF computed. Size: $(size(psf.Q))")
save_psf(joinpath(OUTPUT_DIR, "psf.jld2"), psf)

# =============================================================================
# Step 4: Individual Phase Coupling Functions (Individual PCF)
# =============================================================================
println("\n[Step 4] Computing Individual PCFs")
println("  This may take a while for $N^3 combinations...")
individual_pcfs = IndividualPCF[]

# Optimization: Parallel computation if available, but keep simple loop for script
@time begin
    for i in 1:N, j in 1:N, k in 1:N
        cmd = ComputeIndividualPCF(sps, psf, i, j, k)
        res = interpret(cmd)
        if is_success(res)
            push!(individual_pcfs, unwrap(res))
        end
    end
end
println("  Computed $(length(individual_pcfs)) PCFs")

# =============================================================================
# Step 5: Optimization
# =============================================================================
println("\n[Step 5] Optimization")
Ntheta = sps.Ntheta

# 5a. Linear Stability
println("  [5a] Linear Stability Optimization (target q = -0.1)")
target_q = -0.1
cmd_linear = OptimizeLinearStability(individual_pcfs, target_q, N, Ntheta)
@time result_linear = interpret(cmd_linear)
opt_linear = unwrap(result_linear)
println("    Achieved Λ = $(opt_linear.achieved)")
save_optimization_result(joinpath(OUTPUT_DIR, "opt_linear.jld2"), opt_linear)

# 5b. Rotation
println("  [5b] Rotation Optimization (target μ = 0.5)")
target_mu = 0.5
cmd_rotation = OptimizeRotation(individual_pcfs, target_mu, N, Ntheta)
@time result_rotation = interpret(cmd_rotation)
opt_rotation = unwrap(result_rotation)
println("    Achieved R = $(opt_rotation.achieved)")
save_optimization_result(joinpath(OUTPUT_DIR, "opt_rotation.jld2"), opt_rotation)

# =============================================================================
# Step 6: Total PCF (using Linear Optimization result)
# =============================================================================
println("\n[Step 6] Computing Total PCF (using Linear Opt result)")
cmd_total = ComputeTotalPCF(individual_pcfs, opt_linear.C_opt, N, Ntheta)
@time result_total = interpret(cmd_total)
total_pcf = unwrap(result_total)

metrics = compute_stability_metrics(total_pcf)
println("  Stability Metrics:")
println("    Λ = $(metrics.Lambda)")
println("    R = $(metrics.R)")
save_total_pcf(joinpath(OUTPUT_DIR, "total_pcf_linear.jld2"), total_pcf)

# =============================================================================
# Step 7: Phase Dynamics Simulation
# =============================================================================
println("\n[Step 7] Phase Dynamics Simulation")
epsilon = 2e-4
t_end = 10000.0
phi1_0, phi2_0 = 0.3, -0.2

cmd_phase = SimulatePhaseDynamics(total_pcf, phi1_0, phi2_0, epsilon, 0.0, t_end)
@time result_phase = interpret(cmd_phase)
pts = unwrap(result_phase)

println("  Simulation done.")
println("  Final Phase Diff: φ₁=$(pts.phi1[end]), φ₂=$(pts.phi2[end])")
save_phase_time_series(joinpath(OUTPUT_DIR, "phase_dynamics.jld2"), pts)

println("\n" * "="^60)
println("All steps completed successfully.")
println("Results saved to $OUTPUT_DIR/")
println("="^60)
