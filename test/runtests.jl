"""
    HONetSync テストスイート

FDD層別のテスト構成:
- Unit/: 単体テスト（Pure関数のテスト）
- Integration/: 統合テスト（インタープリターのテスト）
- E2E/: End-to-Endテスト（パイプライン全体のテスト）
"""

using Test
using HONetSync

@testset "HONetSync.jl" begin

    @testset "Domain Layer" begin
        include("Unit/test_domain.jl")
    end

    @testset "Logic Layer" begin
        include("Unit/test_logic.jl")
    end

    @testset "Integration" begin
        include("Integration/test_interpreters.jl")
    end

end

