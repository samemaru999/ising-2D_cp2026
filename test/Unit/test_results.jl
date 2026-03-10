@testset "BoltzmannTable" begin
    bt = BoltzmannTable(Dict(4 => 0.5, 8 => 0.1), 0.5, 1.0)
    @test bt.factors[4] ≈ 0.5
    @test bt.beta == 0.5
    @test bt.J == 1.0
end

@testset "SingleTemperatureResult" begin
    params = IsingParams(4, 1.0)
    avg = ThermalAverages(10.0, 2.0, 5.0, 1.0, 0.5, 100)
    tq = ThermodynamicQuantities(-1.5, 0.8, 2.0, 10.0, 0.6)
    result = SingleTemperatureResult(params, 2.269, avg, tq, 0.5)

    @test result.temperature == 2.269
    @test result.acceptance_rate == 0.5
    @test result.params.L == 4
end

@testset "TemperatureSweepResult" begin
    params = IsingParams(4, 1.0)
    config = SimulationConfig(100, 50, 1)
    temps = [1.0, 2.0, 3.0]
    results = SingleTemperatureResult[]

    sweep_result = TemperatureSweepResult(params, config, temps, results)
    @test length(sweep_result.temperatures) == 3
    @test isempty(sweep_result.results)
end
