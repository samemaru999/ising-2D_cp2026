module Ising2D

# =============================================================================
# Domain Layer [Pure]
# =============================================================================

# 格子生成
include("Ising/Lattice.jl")

# =============================================================================
# Exports
# =============================================================================

export random_lattice, uniform_lattice
export magnetization, energy, delta_energy
export metropolis_step!, sweep!

end # module Ising2D
