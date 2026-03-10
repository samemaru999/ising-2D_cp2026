@testset "compute_thermodynamics" begin
    # 既知の値で手計算と一致することを確認
    avg = ThermalAverages()
    avg.sum_e = -4.0
    avg.sum_e2 = 8.0
    avg.sum_abs_m = 1.8
    avg.sum_m2 = 1.7
    avg.sum_m4 = 1.5
    avg.n_samples = 2

    N = 16  # L=4
    T = 2.0

    thermo = compute_thermodynamics(avg, N, T)

    # mean_e = -4.0/2 = -2.0
    @test thermo.mean_energy ≈ -2.0
    # mean_abs_m = 1.8/2 = 0.9
    @test thermo.mean_abs_magnetization ≈ 0.9
    # mean_e2 = 8.0/2 = 4.0
    # C = N/T^2 * (mean_e2 - mean_e^2) = 16/4 * (4.0 - 4.0) = 0.0
    @test thermo.specific_heat ≈ 0.0
    # mean_m2 = 1.7/2 = 0.85
    # chi = N/T * (mean_m2 - mean_abs_m^2) = 16/2 * (0.85 - 0.81) = 8 * 0.04 = 0.32
    @test thermo.susceptibility ≈ 0.32
    # mean_m4 = 1.5/2 = 0.75
    # U_L = 1 - 0.75/(3*0.85^2) = 1 - 0.75/2.1675 ≈ 0.65397...
    @test thermo.binder_cumulant ≈ 1.0 - 0.75 / (3.0 * 0.85^2)

    # エラーケース
    @test_throws ArgumentError compute_thermodynamics(ThermalAverages(), 16, 2.0)  # n_samples=0
    @test_throws ArgumentError compute_thermodynamics(avg, 0, 2.0)   # N=0
    @test_throws ArgumentError compute_thermodynamics(avg, 16, -1.0) # T<0
end
