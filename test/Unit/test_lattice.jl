@testset "random_lattice" begin
    L = 8
    lat = random_lattice(L)

    @test size(lat) == (L, L)
    @test all(s -> s == 1 || s == -1, lat)
end

@testset "uniform_lattice" begin
    L = 4

    @testset "all up" begin
        lat = uniform_lattice(L)
        @test all(s -> s == 1, lat)
        @test size(lat) == (L, L)
    end

    @testset "all down" begin
        lat = uniform_lattice(L; spin=-1)
        @test all(s -> s == -1, lat)
    end

    @testset "rejects invalid spin" begin
        @test_throws AssertionError uniform_lattice(L; spin=0)
        @test_throws AssertionError uniform_lattice(L; spin=2)
    end
end
