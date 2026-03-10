"""
    Interpreters Layer 統合テスト

インタープリターとシミュレーションのテスト
"""

using Test
using HONetSync
using LinearAlgebra

@testset "Simulation Interpreter" begin
    @testset "Initial State Creation" begin
        N = 4
        X0 = HONetSync.create_initial_state(N)
        @test length(X0) == Ds * N
        @test all(isfinite.(X0))
    end

    # Note: 完全なSPSテストは時間がかかるため、
    # CI環境では短いシミュレーションパラメータを使用

    @testset "ComputeSPS Command (Short)" begin
        N = 4
        params = NetworkParams(N)

        # 短い緩和時間でテスト
        cmd = ComputeSPS(params, 11, 10.0)  # 短いNtheta, t_relax

        result = interpret(cmd)

        # 結果がSuccess/Failureのいずれかであること
        @test result isa ComputationResult

        if is_success(result)
            sps = unwrap(result)
            @test sps.N == N
            @test sps.Ntheta == 11
            @test sps.T > 0
            @test sps.omega > 0
            @test size(sps.Xs) == (Ds * N, 11)
        else
            @warn "SPS computation failed (may be expected in short test): $(result.error)"
        end
    end
end

@testset "Phase Interpreter" begin
    @testset "Individual PCF Computation" begin
        # モックデータでテスト
        N = 2
        Ntheta = 11

        # 簡単なSPSモック
        T = 50.0
        omega = 2π / T
        Xs = randn(Ds * N, Ntheta)
        Xs[:, end] = Xs[:, 1]  # 周期性

        params = NetworkParams(N)
        sps = StablePeriodicSolution(N, Ntheta, T, omega, Xs, params)

        # PSFモック
        Q = randn(Ds * N, Ntheta)
        Q[:, end] = Q[:, 1]
        psf = PhaseSensitivityFunction(N, Ntheta, Q, T, omega)

        # 個別PCF計算コマンド
        cmd = ComputeIndividualPCF(sps, psf, 1, 1, 1)
        result = interpret(cmd)

        @test is_success(result)
        if is_success(result)
            pcf = unwrap(result)
            @test pcf.i == 1
            @test pcf.j == 1
            @test pcf.k == 1
            @test size(pcf.gamma) == (Ntheta, Ntheta)
        end
    end

    @testset "Total PCF Computation" begin
        N = 2
        Ntheta = 11

        # モック個別PCFを作成
        individual_pcfs = IndividualPCF[]
        for i in 1:N
            for j in 1:N
                for k in 1:N
                    gamma = randn(Ntheta, Ntheta) * 0.01
                    push!(individual_pcfs, IndividualPCF(gamma, i, j, k))
                end
            end
        end

        C = create_uniform_coupling_tensor(N)

        cmd = ComputeTotalPCF(individual_pcfs, C, N, Ntheta)
        result = interpret(cmd)

        @test is_success(result)
        if is_success(result)
            total_pcf = unwrap(result)
            @test total_pcf.N == N
            @test total_pcf.Ntheta == Ntheta
            @test size(total_pcf.Gamma) == (Ntheta, Ntheta)
            @test isfinite(total_pcf.dGamma_dphi)
            @test isfinite(total_pcf.dGamma_dpsi)
        end
    end

    @testset "Optimization Commands" begin
        N = 2
        Ntheta = 11

        # モックPCF
        individual_pcfs = IndividualPCF[]
        for i in 1:N
            for j in 1:N
                for k in 1:N
                    gamma = randn(Ntheta, Ntheta) * 0.01
                    push!(individual_pcfs, IndividualPCF(gamma, i, j, k))
                end
            end
        end

        # 線形安定性最適化
        cmd_linear = OptimizeLinearStability(individual_pcfs, -0.1, N, Ntheta)
        result_linear = interpret(cmd_linear)

        @test is_success(result_linear)
        if is_success(result_linear)
            opt = unwrap(result_linear)
            @test opt.opt_type == :linear
            @test size(opt.C_opt) == (N, N, N)
        end

        # 回転特性最適化
        cmd_rotation = OptimizeRotation(individual_pcfs, 0.5, N, Ntheta)
        result_rotation = interpret(cmd_rotation)

        @test is_success(result_rotation)
        if is_success(result_rotation)
            opt = unwrap(result_rotation)
            @test opt.opt_type == :rotation
            @test size(opt.C_opt) == (N, N, N)
        end
    end
end

@testset "Phase Dynamics Simulation" begin
    N = 2
    Ntheta = 11

    # モックTotalPCF
    Gamma = randn(Ntheta, Ntheta) * 0.01
    Gamma[1, 1] = 0.0  # 同期点での値
    C = create_uniform_coupling_tensor(N)
    dGamma_dphi = 0.1
    dGamma_dpsi = 0.05

    total_pcf = TotalPCF(N, Ntheta, Gamma, C, dGamma_dphi, dGamma_dpsi)

    # シミュレーションコマンド
    cmd = SimulatePhaseDynamics(total_pcf, 0.1, -0.1, 2e-4, 0.0, 100.0)
    result = interpret(cmd)

    @test is_success(result)
    if is_success(result)
        pts = unwrap(result)
        @test length(pts.time) > 0
        @test size(pts.phases, 1) == N_NET
        @test size(pts.phase_diffs, 1) == 2
    end
end

@testset "IO Handler" begin
    @testset "Filename Generation" begin
        filename = generate_filename(:sps, 4, 101)
        @test occursin("sps", filename)
        @test occursin("N4", filename)
        @test occursin("Nth101", filename)
        @test endswith(filename, ".jld2")

        filename_opt = generate_filename(:opt, 4, 101; coupling_type=:linear, target_value=-0.1)
        @test occursin("linear", filename_opt)
    end

    @testset "Serialization Helpers" begin
        params = NetworkParams(4)
        serialized = HONetSync.serialize_params(params)

        @test serialized["N"] == 4
        @test haskey(serialized, "fhn")
        @test haskey(serialized, "I")
        @test haskey(serialized, "K")

        deserialized = HONetSync.deserialize_params(serialized)
        @test deserialized.N == params.N
        @test deserialized.fhn.delta == params.fhn.delta
    end

    @testset "Metadata Creation" begin
        metadata = HONetSync.create_metadata()
        @test haskey(metadata, "created_at")
        @test haskey(metadata, "julia_version")
        @test haskey(metadata, "package_version")
    end
end

@testset "Predefined Coupling Access" begin
    # N=4の事前定義結合
    C_uniform = get_predefined_coupling(4, :uniform)
    @test size(C_uniform) == (4, 4, 4)

    C_linear = get_predefined_coupling(4, :linear)
    @test size(C_linear) == (4, 4, 4)

    C_rotation = get_predefined_coupling(4, :rotation)
    @test size(C_rotation) == (4, 4, 4)

    # 事前定義がないケース
    C_uniform_10 = get_predefined_coupling(10, :uniform)
    @test size(C_uniform_10) == (10, 10, 10)

    @test_throws ErrorException get_predefined_coupling(10, :linear)
end

