module Ising2D

# =============================================================================
# Types
# =============================================================================

include("Types.jl")

# =============================================================================
# Ising/ [Pure]
# =============================================================================

include("Ising/Lattice.jl")
include("Ising/Physics.jl")
include("Ising/Statistics.jl")

# =============================================================================
# Runtime/ [Impure]
# =============================================================================

include("Runtime/MonteCarlo.jl")
include("Runtime/IO.jl")
include("Runtime/Visualization.jl")

# =============================================================================
# Exports
# =============================================================================

export IsingParams, SimulationConfig
export Observables, ThermalAverages, ThermodynamicQuantities
export SingleTemperatureResult, TemperatureSweepResult
export BoltzmannTable
export random_lattice, uniform_lattice
export magnetization, energy, delta_energy
export metropolis_step!, sweep!
export save_result, load_result
export plot_spin_lattice, plot_thermodynamics, plot_binder_cumulant, plot_timeseries

end # module Ising2D
