"""
    Logic Layer ユニットテスト

Pure関数（FHN力学、結合、位相縮約、最適化）のテスト
"""

using Test
using HONetSync
using LinearAlgebra
using StaticArrays

@testset "FHN Dynamics" begin
    fhn = FHNParams()

    @testset "Single Oscillator" begin
        # 静止点付近でのダイナミクス
        u, v, I = 0.0, 0.0, 0.0
        du, dv = fhn_single_dynamics(u, v, I, fhn)
        @test du ≈ fhn.delta * fhn.a  # δ(a + v - bu) = δa when u=v=0
        @test dv ≈ 0.0  # v - v³/3 - u + I = 0 when all zero

        # SVector形式
        state = OscillatorState(u, v)
        dstate = fhn_single_dynamics(state, I, fhn)
        @test dstate[1] ≈ du
        @test dstate[2] ≈ dv
    end

    @testset "Single Oscillator Jacobian" begin
        u, v = 0.5, -0.5
        J = fhn_single_jacobian(u, v, fhn)

        @test size(J) == (2, 2)
        @test J[1, 1] ≈ -fhn.delta * fhn.b  # ∂(du)/∂u
        @test J[1, 2] ≈ fhn.delta           # ∂(du)/∂v
        @test J[2, 1] ≈ -1.0                # ∂(dv)/∂u
        @test J[2, 2] ≈ 1.0 - v^2           # ∂(dv)/∂v
    end

    @testset "Network Dynamics" begin
        N = 4
        params = NetworkParams(N)
        X = zeros(Ds * N)
        X[2] = 0.5  # v_1 = 0.5

        # 内部力学のみ
        dX_intrinsic = network_intrinsic_dynamics(X, params)
        @test length(dX_intrinsic) == Ds * N

        # 完全力学（結合込み）
        dX_full = network_full_dynamics(X, params)
        @test length(dX_full) == Ds * N

        # ベクトル場API
        dX_vf = network_vector_field(X, params)
        @test dX_vf ≈ dX_full
    end

    @testset "Network Jacobian" begin
        N = 4
        params = NetworkParams(N)
        X = randn(Ds * N)

        J_intrinsic = network_intrinsic_jacobian(X, params)
        @test size(J_intrinsic) == (Ds * N, Ds * N)

        # ブロック対角構造の確認
        for i in 1:N
            for j in 1:N
                if i != j
                    idx_i = 2i-1:2i
                    idx_j = 2j-1:2j
                    @test all(J_intrinsic[idx_i, idx_j] .== 0)
                end
            end
        end

        J_full = network_full_jacobian(X, params)
        @test size(J_full) == (Ds * N, Ds * N)
    end
end

@testset "Coupling Functions" begin
    @testset "Diffusive Coupling" begin
        v_i, v_j, K_ij = 0.5, 1.0, 0.2
        coupling = diffusive_coupling(v_i, v_j, K_ij)
        @test coupling ≈ K_ij * (v_j - v_i)
    end

    @testset "Higher-Order Interaction" begin
        v_i, v_j, v_k = 0.0, π/4, π/4
        H = higher_order_H(v_i, v_j, v_k)
        @test H ≈ sin(v_j + v_k - 2*v_i)

        H_full = higher_order_H_full(v_i, v_j, v_k)
        @test H_full[1] ≈ 0.0
        @test H_full[2] ≈ H
    end

    @testset "Coupling Matrix Generation" begin
        N = 4

        # 全結合
        K_full = create_full_coupling_matrix(N)
        @test size(K_full) == (N, N)
        @test all(diag(K_full) .== 0)
        @test all(K_full[1, 2:end] .== 1.0)

        # リング結合
        K_ring = create_ring_coupling_matrix(N)
        @test K_ring[1, N] == 1.0  # 周期境界
        @test K_ring[1, 2] == 1.0

        # 対角型テンソル
        C = create_diagonal_coupling_tensor(N)
        @test size(C) == (N, N, N)
        for i in 1:N
            @test C[i, i, 1] > 0
        end
    end

    @testset "Coupling Tensor Properties" begin
        N = 4
        C = create_diagonal_coupling_tensor(N)

        # フロベニウスノルム
        norm_C = compute_coupling_frobenius_norm(C)
        @test norm_C > 0

        # 反対称性チェック
        @test !is_antisymmetric(C)

        # 反対称化
        C_antisym = create_antisymmetric_coupling_tensor(N, C)
        @test is_antisymmetric(C_antisym)
    end
end

@testset "Phase Reduction" begin
    @testset "Stability Metrics" begin
        # 線形安定性
        Gamma1, Gamma2 = -0.1, -0.05
        Lambda = compute_linear_stability(Gamma1, Gamma2)
        @test Lambda ≈ -1.5 * (Gamma1 + Gamma2)
        @test Lambda > 0  # この場合不安定

        Gamma1_stable, Gamma2_stable = 0.1, 0.05
        Lambda_stable = compute_linear_stability(Gamma1_stable, Gamma2_stable)
        @test Lambda_stable < 0  # 安定

        # 回転特性
        R = compute_rotation_characteristic(Gamma1, Gamma2)
        @test R ≈ (sqrt(3)/2) * abs(Gamma1 - Gamma2)
        @test R >= 0
    end

    @testset "Linearized Jacobian" begin
        Gamma1, Gamma2 = 0.1, 0.05
        epsilon = 2e-4

        J = compute_linearized_jacobian(Gamma1, Gamma2, epsilon)
        @test size(J) == (2, 2)

        # 固有値
        lambda1, lambda2 = compute_eigenvalues(Gamma1, Gamma2, epsilon)
        @test real(lambda1) + real(lambda2) ≈ 2 * epsilon * (Gamma1 - Gamma2)
    end

    @testset "Reduced Phase Dynamics (Linear)" begin
        Gamma1, Gamma2 = 0.1, 0.05
        epsilon = 2e-4
        phi1, phi2 = 0.1, -0.1

        dphi1, dphi2 = reduced_phase_dynamics_linear(phi1, phi2, Gamma1, Gamma2, epsilon)

        # 線形近似なのでスケーリングが正しいことを確認
        @test abs(dphi1) < 1e-3  # 小さな変化率
        @test abs(dphi2) < 1e-3
    end
end

@testset "Optimization" begin
    @testset "Analytical Solution Properties" begin
        # モック個別PCFを作成
        N = 2
        Ntheta = 11

        # 単純な個別PCFを手動作成
        pcfs = IndividualPCF[]
        for i in 1:N
            for j in 1:N
                for k in 1:N
                    gamma = randn(Ntheta, Ntheta) * 0.01
                    push!(pcfs, IndividualPCF(gamma, i, j, k))
                end
            end
        end

        # 導関数計算
        gamma1_tensor, gamma2_tensor = compute_all_gamma_derivatives(pcfs, N, Ntheta)
        @test size(gamma1_tensor) == (N, N, N)
        @test size(gamma2_tensor) == (N, N, N)
    end

    @testset "Target Value Calculations" begin
        # 所望の安定性を達成するためのq値
        target_Lambda = -1.0
        q = find_optimal_q_for_stability(IndividualPCF[], target_Lambda, 4, 101)
        @test q ≈ -2/3 * target_Lambda

        # 所望の回転特性を達成するためのμ値
        target_R = 0.5
        mu = find_optimal_mu_for_rotation(target_R)
        @test mu ≈ 2/sqrt(3) * target_R
    end
end

@testset "Three Network Dynamics" begin
    N = 4
    params = NetworkParams(N)
    C = create_uniform_coupling_tensor(N)
    epsilon = 2e-4
    three_params = ThreeNetworkParams(params, C, epsilon)

    X_A = randn(Ds * N) * 0.1
    X_B = randn(Ds * N) * 0.1
    X_C = randn(Ds * N) * 0.1

    # 3ネットワーク力学
    dX_A, dX_B, dX_C = three_network_dynamics(X_A, X_B, X_C, three_params)

    @test length(dX_A) == Ds * N
    @test length(dX_B) == Ds * N
    @test length(dX_C) == Ds * N

    # ThreeNetworkState API
    state = ThreeNetworkState(X_A, X_B, X_C)
    dX_A2, dX_B2, dX_C2 = compute_three_network_dynamics(state, three_params)

    @test dX_A2 ≈ dX_A
    @test dX_B2 ≈ dX_B
    @test dX_C2 ≈ dX_C
end

@testset "Helper Functions" begin
    N = 4
    X = collect(1.0:8.0)  # [1,2,3,4,5,6,7,8]

    # v成分抽出
    v = extract_v_components(X, N)
    @test v == [2.0, 4.0, 6.0, 8.0]

    # u成分抽出
    u = extract_u_components(X, N)
    @test u == [1.0, 3.0, 5.0, 7.0]

    # 状態ベクトル構成
    X_reconstructed = create_state_vector(u, v)
    @test X_reconstructed == X
end

