@testset "MC simulation workflow" begin
    Random.seed!(12345)
    L = 8
    N = L * L
    J = 1.0
    T = 2.0
    beta = 1.0 / T

    lat = random_lattice(L)

    # 熱平衡化
    for _ in 1:500
        sweep!(lat, beta; J=J)
    end

    # 測定
    avg = ThermalAverages()
    for _ in 1:200
        sweep!(lat, beta; J=J)
        e = energy(lat; J=J) / N
        m = magnetization(lat)
        avg.sum_e += e
        avg.sum_e2 += e^2
        avg.sum_abs_m += abs(m)
        avg.sum_m2 += m^2
        avg.sum_m4 += m^4
        avg.n_samples += 1
    end

    mean_e = avg.sum_e / avg.n_samples
    mean_abs_m = avg.sum_abs_m / avg.n_samples

    # T < Tc なので秩序相: |m| > 0, e < 0
    @test mean_e < 0
    @test mean_abs_m > 0.3
    @test avg.n_samples == 200
end
