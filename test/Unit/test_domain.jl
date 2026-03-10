"""
    Domain Layer ユニットテスト

Pure関数とドメイン型のテスト
"""

using Test
using HONetSync
using LinearAlgebra

@testset "Constants" begin
    @test N_OSC_DEFAULT == 4
    @test N_NET == 3
    @test Ds == 2
    @test Nθ_DEFAULT == 101

    @test δ_DEFAULT > 0
    @test 0 < b_DEFAULT < 1

    @test size(INTRA_K_4) == (4, 4)
    @test size(INTRA_K_10) == (10, 10)

    # 対角成分がゼロであること
    @test all(diag(INTRA_K_4) .== 0)
    @test all(diag(INTRA_K_10) .== 0)
end

@testset "FHNParams" begin
    # 正常な作成
    fhn = FHNParams(0.08, 0.7, 0.8)
    @test fhn.delta == 0.08
    @test fhn.a == 0.7
    @test fhn.b == 0.8

    # デフォルトコンストラクタ
    fhn_default = FHNParams()
    @test fhn_default.delta == δ_DEFAULT
    @test fhn_default.a == a_DEFAULT
    @test fhn_default.b == b_DEFAULT

    # バリデーション: delta > 0
    @test_throws DomainError FHNParams(-0.1, 0.7, 0.8)
    @test_throws DomainError FHNParams(0.0, 0.7, 0.8)

    # バリデーション: 0 < b < 1
    @test_throws DomainError FHNParams(0.08, 0.7, 0.0)
    @test_throws DomainError FHNParams(0.08, 0.7, 1.0)
    @test_throws DomainError FHNParams(0.08, 0.7, 1.5)
end

@testset "NetworkParams" begin
    N = 4
    fhn = FHNParams()
    I = get_default_external_input(N)
    K = get_intra_K(N)

    # 正常な作成
    params = NetworkParams(N, fhn, I, K, FullCoupling())
    @test params.N == N
    @test params.fhn == fhn
    @test length(params.I) == N
    @test size(params.K) == (N, N)

    # デフォルトコンストラクタ
    params_default = NetworkParams(N)
    @test params_default.N == N

    # バリデーション: N > 0
    @test_throws DomainError NetworkParams(0, fhn, I, K, FullCoupling())

    # バリデーション: 次元不整合
    @test_throws DimensionMismatch NetworkParams(N, fhn, zeros(3), K, FullCoupling())
    @test_throws DimensionMismatch NetworkParams(N, fhn, I, zeros(3, 3), FullCoupling())
end

@testset "ThreeNetworkParams" begin
    N = 4
    params = NetworkParams(N)
    C = create_uniform_coupling_tensor(N)
    epsilon = 2e-4

    three_params = ThreeNetworkParams(params, C, epsilon)
    @test three_params.network.N == N
    @test size(three_params.C) == (N, N, N)
    @test three_params.epsilon == epsilon

    # バリデーション: 次元不整合
    C_wrong = zeros(3, 3, 3)
    @test_throws DimensionMismatch ThreeNetworkParams(params, C_wrong, epsilon)
end

@testset "TopologyType" begin
    # FullCoupling
    @test FullCoupling() isa TopologyType

    # TreeCoupling
    parent_indices = [0, 1, 1, 2]  # ルート, 1の子, 1の子, 2の子
    tree = TreeCoupling(parent_indices)
    @test tree isa TopologyType
    @test tree.parent_indices == parent_indices

    # 不正な親インデックス
    @test_throws DomainError TreeCoupling([0, 5, 1, 2])  # 範囲外
    @test_throws DomainError TreeCoupling([1, 2, 3, 4])  # 自己参照
end

@testset "eDSL Commands" begin
    N = 4
    params = NetworkParams(N)

    # ComputeSPS
    cmd_sps = ComputeSPS(params)
    @test cmd_sps.Ntheta == Nθ_DEFAULT
    @test cmd_sps.t_relax == TREL

    cmd_sps2 = ComputeSPS(params, 51, 100.0)
    @test cmd_sps2.Ntheta == 51
    @test cmd_sps2.t_relax == 100.0

    # バリデーション
    @test_throws DomainError ComputeSPS(params, 0, 100.0)
    @test_throws DomainError ComputeSPS(params, 101, -1.0)
end

@testset "ComputationResult" begin
    # Success
    result_success = Success(42)
    @test is_success(result_success)
    @test unwrap(result_success) == 42
    @test unwrap_or(result_success, 0) == 42

    # Failure
    result_failure = Failure{Int}("error message")
    @test !is_success(result_failure)
    @test_throws ErrorException unwrap(result_failure)
    @test unwrap_or(result_failure, 0) == 0

    # map_result
    doubled = map_result(x -> x * 2, result_success)
    @test is_success(doubled)
    @test unwrap(doubled) == 84

    failed_map = map_result(x -> x * 2, result_failure)
    @test !is_success(failed_map)
end

@testset "Utility Functions" begin
    # get_intra_K
    K4 = get_intra_K(4)
    @test size(K4) == (4, 4)
    K10 = get_intra_K(10)
    @test size(K10) == (10, 10)
    @test_throws ErrorException get_intra_K(5)

    # get_default_external_input
    I = get_default_external_input(4)
    @test length(I) == 4
    @test I[1] == I_PACEMAKER
    @test all(I[2:end] .== I_DEFAULT)

    # create_uniform_coupling_tensor
    C = create_uniform_coupling_tensor(4)
    @test size(C) == (4, 4, 4)

    # index_to_phase / phase_to_index
    Ntheta = 101
    @test index_to_phase(1, Ntheta) ≈ 0.0
    @test index_to_phase(Ntheta, Ntheta) ≈ 2π

    @test phase_to_index(0.0, Ntheta) == 1
    @test phase_to_index(2π, Ntheta) == Ntheta

    # normalize_phase
    @test normalize_phase(0.0) ≈ 0.0
    @test normalize_phase(π) ≈ π
    @test normalize_phase(3π) ≈ π
    @test normalize_phase(-π/2) ≈ -π/2

    # normalize_phase_positive
    @test normalize_phase_positive(0.0) ≈ 0.0
    @test normalize_phase_positive(2π) ≈ 0.0 atol=1e-10
    @test normalize_phase_positive(-π/2) ≈ 3π/2
end

@testset "Validation Functions" begin
    # validate_fhn_params
    @test validate_fhn_params(0.08, 0.7, 0.8) === nothing
    @test validate_fhn_params(-0.1, 0.7, 0.8) !== nothing
    @test validate_fhn_params(0.08, 0.7, 1.5) !== nothing

    # validate_network_params
    N = 4
    I = get_default_external_input(N)
    K = get_intra_K(N)
    @test validate_network_params(N, I, K) === nothing
    @test validate_network_params(0, I, K) !== nothing
    @test validate_network_params(N, zeros(3), K) !== nothing

    # validate_coupling_tensor
    C = create_uniform_coupling_tensor(4)
    @test validate_coupling_tensor(C, 4) === nothing
    @test validate_coupling_tensor(C, 5) !== nothing

    C_nan = copy(C)
    C_nan[1, 1, 1] = NaN
    @test validate_coupling_tensor(C_nan, 4) !== nothing

    # validate_epsilon
    @test validate_epsilon(2e-4) === nothing
    @test validate_epsilon(-1e-4) !== nothing
    @test validate_epsilon(0.5) !== nothing
end

