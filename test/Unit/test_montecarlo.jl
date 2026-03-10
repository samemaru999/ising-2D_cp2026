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
