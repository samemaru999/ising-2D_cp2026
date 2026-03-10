"""
Ising2D の全型定義。
data_structure.md に対応する。
"""

# =============================================================================
# パラメータ型
# =============================================================================

struct IsingParams
    L::Int
    J::Float64
    h::Float64

    function IsingParams(L::Int, J::Float64, h::Float64)
        L >= 2 || throw(ArgumentError("L must be ≥ 2, got $L"))
        J != 0.0 || throw(ArgumentError("J must be non-zero"))
        new(L, J, h)
    end
end

struct SimulationConfig
    n_sweeps::Int
    n_equilibration::Int
    sample_interval::Int
    seed::Union{Int, Nothing}

    function SimulationConfig(n_sweeps::Int, n_equilibration::Int,
                              sample_interval::Int, seed::Union{Int, Nothing}=nothing)
        n_sweeps > 0 || throw(ArgumentError("n_sweeps must be > 0, got $n_sweeps"))
        n_equilibration >= 0 || throw(ArgumentError("n_equilibration must be ≥ 0, got $n_equilibration"))
        sample_interval > 0 || throw(ArgumentError("sample_interval must be > 0, got $sample_interval"))
        new(n_sweeps, n_equilibration, sample_interval, seed)
    end
end

# =============================================================================
# 物理量の型
# =============================================================================

struct Observables
    energy_per_site::Float64
    magnetization_per_site::Float64
end

mutable struct ThermalAverages
    sum_e::Float64
    sum_e2::Float64
    sum_abs_m::Float64
    sum_m2::Float64
    sum_m4::Float64
    n_samples::Int
end

ThermalAverages() = ThermalAverages(0.0, 0.0, 0.0, 0.0, 0.0, 0)

struct ThermodynamicQuantities
    mean_energy::Float64
    mean_abs_magnetization::Float64
    specific_heat::Float64
    susceptibility::Float64
    binder_cumulant::Float64
end

# =============================================================================
# 結果型
# =============================================================================

struct SingleTemperatureResult
    params::IsingParams
    temperature::Float64
    averages::ThermalAverages
    thermodynamics::ThermodynamicQuantities
    acceptance_rate::Float64
end

struct TemperatureSweepResult
    params::IsingParams
    config::SimulationConfig
    temperatures::Vector{Float64}
    results::Vector{SingleTemperatureResult}
end

# =============================================================================
# 性能最適化型
# =============================================================================

struct BoltzmannTable
    factors::Dict{Int, Float64}
    beta::Float64
    J::Float64
end
