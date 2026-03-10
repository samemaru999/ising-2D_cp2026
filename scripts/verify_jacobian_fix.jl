#!/usr/bin/env julia
"""
    verify_jacobian_fix.jl - Verify that Jacobian fix matches paper equation (11)

Tests:
1. Eigenvalue formula verification
2. Convergence simulation with Γ₁ + Γ₂ > 0
3. Divergence simulation with Γ₁ + Γ₂ < 0
"""

using HONetSync
using LinearAlgebra
using DifferentialEquations

println("Phase Dynamics Verification: Corrected Jacobian (NOLTA2025 Eq.11)")
println("="^60)

# Test cases
test_cases = [
    ("Strong convergence", 0.1, 0.1),      # Γ₁+Γ₂ = 0.2 > 0, Λ = -0.3 < 0
    ("Weak convergence", 0.05, 0.05),      # Γ₁+Γ₂ = 0.1 > 0, Λ = -0.15 < 0
    ("Asymmetric convergence", 0.15, 0.05),# Γ₁+Γ₂ = 0.2 > 0, Λ = -0.3 < 0
    ("Divergence", -0.1, -0.1),            # Γ₁+Γ₂ = -0.2 < 0, Λ = 0.3 > 0
]

println("\n1. Eigenvalue Formula Verification")
println("-"^60)

for (name, G1, G2) in test_cases
    println("\nCase: $name (Γ₁=$G1, Γ₂=$G2)")

    # Theoretical values
    Lambda_theory = -1.5 * (G1 + G2)
    R_theory = (sqrt(3) / 2) * abs(G1 - G2)

    # Numerical eigenvalues (ε=1 for comparison)
    J = compute_linearized_jacobian(G1, G2, 1.0)
    eigenvals = eigvals(J)

    Re_lambda = real(eigenvals[1])
    Im_lambda = abs(imag(eigenvals[1]))

    println("  Theory: Λ = $Lambda_theory, R = $R_theory")
    println("  Numerical: Re(λ) = $Re_lambda, |Im(λ)| = $Im_lambda")
    println("  Match: Re=$(isapprox(Re_lambda, Lambda_theory, atol=1e-10)), Im=$(isapprox(Im_lambda, R_theory, atol=1e-10))")
    println("  Stable: $(Lambda_theory < 0)")
end

println("\n\n2. Convergence Simulation")
println("-"^60)

# Simulation parameters
epsilon = 1.0
t_end = 50.0
phi0 = [0.3, -0.2]  # Initial phase differences

function simulate_phase_dynamics(G1, G2, epsilon, phi0, t_end)
    function dynamics!(dphi, phi, p, t)
        dphi1, dphi2 = reduced_phase_dynamics_linear(phi[1], phi[2], G1, G2, epsilon)
        dphi[1] = dphi1
        dphi[2] = dphi2
    end

    prob = ODEProblem(dynamics!, phi0, (0.0, t_end))
    sol = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-10, saveat=0.5)
    return sol
end

for (name, G1, G2) in test_cases
    println("\nCase: $name (Γ₁=$G1, Γ₂=$G2)")

    sol = simulate_phase_dynamics(G1, G2, epsilon, phi0, t_end)

    initial_dist = norm(phi0)
    final_phi = [sol.u[end][1], sol.u[end][2]]
    final_dist = norm(final_phi)
    ratio = final_dist / initial_dist * 100

    println("  Initial: φ = ($(phi0[1]), $(phi0[2])), dist = $initial_dist")
    println("  Final:   φ = ($(round(final_phi[1], digits=6)), $(round(final_phi[2], digits=6))), dist = $(round(final_dist, digits=6))")
    println("  Ratio: $(round(ratio, digits=2))%")

    Lambda = -1.5 * (G1 + G2)
    if Lambda < 0
        println("  Expected: Convergence (Λ = $Lambda < 0)")
        println("  Result: $(final_dist < initial_dist ? "✓ CONVERGED" : "✗ DID NOT CONVERGE")")
    else
        println("  Expected: Divergence (Λ = $Lambda > 0)")
        println("  Result: $(final_dist > initial_dist ? "✓ DIVERGED" : "✗ DID NOT DIVERGE")")
    end
end

println("\n\n3. Jacobian Matrix Check")
println("-"^60)

G1, G2 = 0.1, 0.1
J = compute_linearized_jacobian(G1, G2, 1.0)
println("For Γ₁=$G1, Γ₂=$G2 (ε=1):")
println("J = ")
println("  [$(J[1,1])  $(J[1,2])]")
println("  [$(J[2,1])  $(J[2,2])]")
println("\nExpected (from paper Eq.11):")
println("J = [-2Γ₁-Γ₂    Γ₁-Γ₂ ] = [-$(2*G1+G2)   $(G1-G2)]")
println("    [ Γ₂-Γ₁   -2Γ₂-Γ₁ ]   [ $(G2-G1)  -$(2*G2+G1)]")

println("\n\nVerification complete!")
println("Press Enter to close...")
readline()
