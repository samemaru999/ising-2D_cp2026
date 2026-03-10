@testset "Observables" begin
    obs = Observables(-1.5, 0.8)
    @test obs.energy_per_site == -1.5
    @test obs.magnetization_per_site == 0.8
end

@testset "ThermalAverages" begin
    @testset "zero constructor" begin
        avg = ThermalAverages()
        @test avg.sum_e == 0.0
        @test avg.sum_e2 == 0.0
        @test avg.sum_abs_m == 0.0
        @test avg.sum_m2 == 0.0
        @test avg.sum_m4 == 0.0
        @test avg.n_samples == 0
    end

    @testset "mutable accumulation" begin
        avg = ThermalAverages()
        avg.sum_e += 1.5
        avg.n_samples += 1
        @test avg.sum_e == 1.5
        @test avg.n_samples == 1
    end
end

@testset "ThermodynamicQuantities" begin
    tq = ThermodynamicQuantities(-1.5, 0.8, 2.0, 10.0, 0.6)
    @test tq.mean_energy == -1.5
    @test tq.mean_abs_magnetization == 0.8
    @test tq.specific_heat == 2.0
    @test tq.susceptibility == 10.0
    @test tq.binder_cumulant == 0.6
end
