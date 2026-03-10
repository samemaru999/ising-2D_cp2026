using Test

include("../../src/Ising/Lattice.jl")

@testset "random_lattice" begin
    lattice = random_lattice(8)
    @test size(lattice) == (8, 8)
    @test eltype(lattice) == Int
    @test all(s -> s == 1 || s == -1, lattice)
end

@testset "uniform_lattice" begin
    up = uniform_lattice(4)
    @test all(s -> s == 1, up)
    @test size(up) == (4, 4)

    down = uniform_lattice(4; spin=-1)
    @test all(s -> s == -1, down)

    @test_throws AssertionError uniform_lattice(4; spin=0)
end
