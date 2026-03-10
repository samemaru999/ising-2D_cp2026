@testset "metropolis_step!" begin
    @testset "T = ∞ (beta = 0) → always accept" begin
        Random.seed!(42)
        L = 4
        lat = uniform_lattice(L)
        n_accept = sum(metropolis_step!(copy(lat), 0.0) for _ in 1:100)
        @test n_accept == 100
    end

    @testset "preserves spin values ±1" begin
        Random.seed!(42)
        lat = random_lattice(8)
        for _ in 1:100
            metropolis_step!(lat, 1.0)
        end
        @test all(s -> s == 1 || s == -1, lat)
    end
end

@testset "sweep!" begin
    @testset "returns accepted count in [0, N]" begin
        Random.seed!(42)
        L = 8
        N = L * L
        lat = random_lattice(L)
        accepted = sweep!(lat, 0.5)
        @test 0 <= accepted <= N
    end

    @testset "T = ∞ → all accepted" begin
        Random.seed!(42)
        L = 4
        lat = random_lattice(L)
        accepted = sweep!(lat, 0.0)
        @test accepted == L * L
    end

    @testset "low T, all-up → stays ordered" begin
        Random.seed!(42)
        L = 8
        lat = uniform_lattice(L)
        beta = 10.0  # T = 0.1
        for _ in 1:100
            sweep!(lat, beta)
        end
        @test abs(magnetization(lat)) > 0.9
    end
end

@testset "quick_mc" begin
    @testset "returns Dict with expected keys" begin
        Random.seed!(42)
        result = quick_mc(2.0, 4)
        @test result isa Dict{String, Float64}
        @test haskey(result, "E")
        @test haskey(result, "M")
        @test haskey(result, "M2")
    end

    @testset "physical value ranges" begin
        Random.seed!(42)
        result = quick_mc(2.0, 4)
        @test -2.0 <= result["E"] <= 0.0
        @test 0.0 <= result["M"] <= 1.0
        @test 0.0 <= result["M2"] <= 1.0
    end

    @testset "low T → high M2" begin
        Random.seed!(42)
        result = quick_mc(0.5, 8; n_sweeps=1000, n_equilibration=500)
        @test result["M2"] > 0.8
    end

    @testset "high T → low M2" begin
        Random.seed!(42)
        result = quick_mc(10.0, 8; n_sweeps=1000, n_equilibration=500)
        @test result["M2"] < 0.2
    end

    @testset "invalid arguments" begin
        @test_throws ArgumentError quick_mc(-1.0, 4)
        @test_throws ArgumentError quick_mc(0.0, 4)
        @test_throws ArgumentError quick_mc(2.0, 1)
    end
end
