@testset "magnetization" begin
    @testset "all up → m = 1" begin
        lat = uniform_lattice(4)
        @test magnetization(lat) ≈ 1.0
    end

    @testset "all down → m = -1" begin
        lat = uniform_lattice(4; spin=-1)
        @test magnetization(lat) ≈ -1.0
    end

    @testset "checkerboard → m = 0" begin
        L = 4
        lat = [iseven(i + j) ? 1 : -1 for i in 1:L, j in 1:L]
        @test magnetization(lat) ≈ 0.0
    end
end

@testset "energy" begin
    @testset "all up → E = -2JN" begin
        L = 4
        N = L * L
        J = 1.0
        lat = uniform_lattice(L)
        @test energy(lat; J=J) ≈ -2.0 * J * N
    end

    @testset "all down → E = -2JN (same as all up)" begin
        L = 4
        N = L * L
        lat = uniform_lattice(L; spin=-1)
        @test energy(lat; J=1.0) ≈ -2.0 * N
    end

    @testset "checkerboard → E = +2JN" begin
        L = 4
        N = L * L
        lat = [iseven(i + j) ? 1 : -1 for i in 1:L, j in 1:L]
        @test energy(lat; J=1.0) ≈ 2.0 * N
    end

    @testset "scales with J" begin
        L = 4
        lat = uniform_lattice(L)
        @test energy(lat; J=2.0) ≈ 2.0 * energy(lat; J=1.0)
    end
end

@testset "delta_energy" begin
    @testset "consistency with energy" begin
        L = 8
        lat = random_lattice(L)
        J = 1.0
        i, j = 3, 5

        E_before = energy(lat; J=J)
        dE = delta_energy(lat, i, j; J=J)
        lat[i, j] *= -1
        E_after = energy(lat; J=J)

        @test dE ≈ E_after - E_before
    end

    @testset "all up, flip center → ΔE = 8J" begin
        L = 4
        lat = uniform_lattice(L)
        @test delta_energy(lat, 2, 2; J=1.0) ≈ 8.0
    end

    @testset "all up, flip corner (periodic BC) → ΔE = 8J" begin
        L = 4
        lat = uniform_lattice(L)
        @test delta_energy(lat, 1, 1; J=1.0) ≈ 8.0
    end
end
