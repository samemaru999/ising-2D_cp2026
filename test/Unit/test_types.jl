@testset "IsingParams" begin
    @testset "valid construction" begin
        p = IsingParams(32, 1.0)
        @test p.L == 32
        @test p.J == 1.0
    end

    @testset "rejects L < 2" begin
        @test_throws ArgumentError IsingParams(1, 1.0)
        @test_throws ArgumentError IsingParams(0, 1.0)
        @test_throws ArgumentError IsingParams(-1, 1.0)
    end

    @testset "rejects J = 0" begin
        @test_throws ArgumentError IsingParams(32, 0.0)
    end

    @testset "accepts negative J" begin
        p = IsingParams(4, -1.0)
        @test p.J == -1.0
    end
end

@testset "SimulationConfig" begin
    @testset "valid construction" begin
        c = SimulationConfig(1000, 500, 5)
        @test c.n_sweeps == 1000
        @test c.n_equilibration == 500
        @test c.sample_interval == 5
        @test c.seed === nothing
    end

    @testset "with seed" begin
        c = SimulationConfig(100, 0, 1, 42)
        @test c.seed == 42
    end

    @testset "rejects invalid values" begin
        @test_throws ArgumentError SimulationConfig(0, 0, 1)
        @test_throws ArgumentError SimulationConfig(100, -1, 1)
        @test_throws ArgumentError SimulationConfig(100, 0, 0)
    end
end
