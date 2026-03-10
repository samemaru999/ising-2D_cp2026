using Test
using Random
using Ising2D

@testset "Ising2D" begin

    # Unit tests (80%)
    @testset "Types" begin
        include("Unit/test_types.jl")
    end

    @testset "Observables & Averages" begin
        include("Unit/test_observables.jl")
    end

    @testset "Results" begin
        include("Unit/test_results.jl")
    end

    @testset "Lattice" begin
        include("Unit/test_lattice.jl")
    end

    @testset "Physics" begin
        include("Unit/test_physics.jl")
    end

    @testset "MonteCarlo" begin
        include("Unit/test_montecarlo.jl")
    end

    # Integration / E2E (20%)
    @testset "Integration" begin
        include("Ising/test_integration.jl")
    end

    @testset "E2E" begin
        include("E2E/test_simulation.jl")
    end

end
