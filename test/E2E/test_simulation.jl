@testset "high T vs low T" begin
    Random.seed!(99)
    L = 8
    N = L * L
    J = 1.0

    function measure(T, n_eq, n_meas)
        beta = 1.0 / T
        lat = random_lattice(L)
        for _ in 1:n_eq
            sweep!(lat, beta; J=J)
        end
        m_sum = 0.0
        for _ in 1:n_meas
            sweep!(lat, beta; J=J)
            m_sum += abs(magnetization(lat))
        end
        return m_sum / n_meas
    end

    m_low = measure(1.0, 500, 200)   # T << Tc → ordered
    m_high = measure(4.0, 500, 200)  # T >> Tc → disordered

    @test m_low > 0.8
    @test m_high < 0.3
end
