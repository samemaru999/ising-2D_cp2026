"""
メトロポリス法によるモンテカルロ更新。
乱数を使用するため Impure。
"""

"""
    metropolis_step!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0) -> Bool

ランダムに1サイトを選び、メトロポリス法でスピン反転を試みる。
反転が受容された場合 true を返す。
"""
function metropolis_step!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0)
    L = size(lattice, 1)
    i = rand(1:L)
    j = rand(1:L)

    dE = delta_energy(lattice, i, j; J=J)

    if dE <= 0 || rand() < exp(-beta * dE)
        lattice[i, j] *= -1
        return true
    end
    return false
end

"""
    sweep!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0) -> Int

L×L 回のメトロポリス更新試行を行う（1モンテカルロスイープ）。
受容された回数を返す。
"""
function sweep!(lattice::Matrix{Int}, beta::Float64; J::Float64=1.0)
    L = size(lattice, 1)
    accepted = 0
    for _ in 1:L*L
        accepted += metropolis_step!(lattice, beta; J=J)
    end
    return accepted
end

"""
    run_single_temperature(lattice::Matrix{Int}, beta::Float64, N::Int, config::SimulationConfig;
                           J::Float64=1.0) -> (ThermalAverages, Int)

平衡化 + 測定ループを実行し、(ThermalAverages, total_accepted) を返す。
lattice は in-place で更新される。
"""
function run_single_temperature(lattice::Matrix{Int}, beta::Float64, N::Int, config::SimulationConfig;
                                J::Float64=1.0)
    # 平衡化
    for _ in 1:config.n_equilibration
        sweep!(lattice, beta; J=J)
    end

    # 測定
    avg = ThermalAverages()
    total_accepted = 0

    for sweep_i in 1:config.n_sweeps
        accepted = sweep!(lattice, beta; J=J)
        total_accepted += accepted

        if sweep_i % config.sample_interval == 0
            e = energy(lattice; J=J) / N
            m = magnetization(lattice)

            avg.sum_e += e
            avg.sum_e2 += e^2
            avg.sum_abs_m += abs(m)
            avg.sum_m2 += m^2
            avg.sum_m4 += m^4
            avg.n_samples += 1
        end
    end

    return (avg, total_accepted)
end

"""
    run_temperature_sweep(params::IsingParams, config::SimulationConfig,
                          temperatures::Vector{Float64}) -> TemperatureSweepResult

温度配列をループし、各温度でシミュレーションを実行する。
高温側から開始し、前の温度の格子状態を引き継ぐ。
"""
function run_temperature_sweep(params::IsingParams, config::SimulationConfig,
                               temperatures::Vector{Float64})
    N = params.L^2
    sorted_temps = sort(temperatures; rev=true)
    lattice = random_lattice(params.L)

    results = SingleTemperatureResult[]

    for (i, T) in enumerate(sorted_temps)
        beta = 1.0 / T
        println("  [$i/$(length(sorted_temps))] T = $(round(T, digits=4)) ...")

        avg, total_accepted = run_single_temperature(lattice, beta, N, config; J=params.J)
        thermo = compute_thermodynamics(avg, N, T)
        acceptance_rate = total_accepted / (config.n_sweeps * N)

        push!(results, SingleTemperatureResult(params, T, avg, thermo, acceptance_rate))
    end

    # 温度昇順に並べ替え
    sort!(results; by=r -> r.temperature)
    sorted_final_temps = sort(sorted_temps)

    return TemperatureSweepResult(params, config, sorted_final_temps, results)
end

"""
    quick_mc(T::Float64, L::Int; J::Float64=1.0, n_sweeps::Int=2000,
             n_equilibration::Int=500, sample_interval::Int=5) -> Dict{String, Float64}

温度 T と格子サイズ L を受け取り、短いMCシミュレーションを実行して
⟨E⟩, ⟨|M|⟩, ⟨M²⟩ の平均値を Dict で返す簡易インターフェース。
"""
function quick_mc(T::Float64, L::Int;
                  J::Float64=1.0, n_sweeps::Int=2000,
                  n_equilibration::Int=500, sample_interval::Int=5)
    T > 0 || throw(ArgumentError("T must be > 0, got $T"))
    L >= 2 || throw(ArgumentError("L must be ≥ 2, got $L"))

    config = SimulationConfig(n_sweeps, n_equilibration, sample_interval)
    lattice = random_lattice(L)
    N = L^2
    beta = 1.0 / T

    avg, _ = run_single_temperature(lattice, beta, N, config; J=J)

    n = avg.n_samples
    return Dict(
        "E"  => avg.sum_e / n,
        "M"  => avg.sum_abs_m / n,
        "M2" => avg.sum_m2 / n
    )
end
