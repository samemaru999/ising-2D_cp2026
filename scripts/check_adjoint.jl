
using HONetSync
using LinearAlgebra
using Printf

function check_psf_residual()
    println("Checking PSF Residual...")
    flush(stdout)

    # Setup
    N = 4
    Ntheta = 101
    params = NetworkParams(N)

    # 1. Compute SPS
    println("Computing SPS (t_relax=500)...")
    cmd_sps = ComputeSPS(params, Ntheta, 500.0)
    res_sps = interpret(cmd_sps)
    if !is_success(res_sps)
        println("SPS failed: $(res_sps.error)")
        return
    end
    sps = unwrap(res_sps)
    println("SPS computed. T=$(sps.T)")

    # 2. Compute PSF
    println("Computing PSF (20 iterations)...")
    cmd_psf = ComputePSF(sps, 20)
    res_psf = interpret(cmd_psf)
    if !is_success(res_psf)
        println("PSF failed: $(res_psf.error)")
        return
    end
    psf = unwrap(res_psf)
    println("PSF computed.")

    # 3. Check Residual
    # Equation: dQ/dt = -J^T Q
    # dQ/dtheta * omega = -J^T Q
    # (Q_{k+1} - Q_k) / (dt * omega) approx -J^T Q
    # We use Forward Diff in time (Backward in index) to match Euler?
    # Code used: Q_{k-1} = Q_k + dt * J_k^T Q_k
    # So (Q_{k-1} - Q_k)/dt = J_k^T Q_k  => -(Q_k - Q_{k-1})/dt = J^T Q
    # Residual = (Q_{k-1} - Q_k)/dt - J_k^T Q_k

    omega = sps.omega
    T = sps.T
    dt = T / (Ntheta - 1)

    residuals = Float64[]

    for k in Ntheta:-1:2
        # Data at k
        Q_k = psf.Q[:, k]
        X_k = sps.Xs[:, k]
        # SimulationInterpreter calls network_full_jacobian

        # Calculate J_k manually or use internal function if exposed?
        # Since I cannot call internal function easily from script without `using HONetSync.Logic` (if exported)
        # Assuming exported or available. FHN.jl exports `network_full_jacobian`.

        J_k_mat = HONetSync.network_full_jacobian(X_k, params)

        # Data at k-1
        Q_km1 = psf.Q[:, k-1]

        lhs = (Q_km1 - Q_k) / dt
        rhs = (J_k_mat' * Q_k)

        res = norm(lhs - rhs)
        push!(residuals, res)
    end

    max_res = maximum(residuals)
    mean_res = sum(residuals) / length(residuals)

    println("Max Residual: $max_res")
    println("Mean Residual: $mean_res")

    # Check magnitude of derivative
    avg_mag = sum(norm.(eachcol(psf.Q))) / Ntheta
    println("Avg |Q|: $avg_mag")

    rel_res = mean_res / avg_mag
    println("Relative Residual: $(rel_res * 100)%")

    if rel_res > 0.01
        println("CONCLUSION: Residual is significant (>1%). Accuracy might be low.")
    else
        println("CONCLUSION: Residual is small. Accuracy matches Euler scheme.")
    end
end

check_psf_residual()
