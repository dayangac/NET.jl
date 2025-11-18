"""
NET.jl - Non-Equilibrium Thermodynamics Library

A Julia library for working with Non-Equilibrium Thermodynamics (NET) principles,
including Onsager Relations, Entropy Production, and Stochastic Thermodynamics.
"""
module NET

# External dependencies
using Plots

# Utility modules (must be included first)
include("utils/ode_solver.jl")
include("utils/gradient.jl")

# Submodules - Onsager Relations
include("onsager/linear_onsager.jl")
include("onsager/non_linear_onsager.jl")
include("onsager/stochastic_onsager_relations.jl")
include("onsager/time_dependent_onsager.jl")

# Submodules - Entropy Production
include("entropy_production/local_entropy_production.jl")
include("entropy_production/entropy_production_rate.jl")
include("entropy_production/time_dependent_entropy.jl")
include("entropy_production/entropy_optimization.jl")

# Submodules - Stochastic Thermodynamics
include("stohcastic_thermodynamics/fluctuation_theorems.jl")
include("stohcastic_thermodynamics/fokker_planck.jl")
include("stohcastic_thermodynamics/langevin_dynamics.jl")
include("stohcastic_thermodynamics/noise_induced_bifurications.jl")
include("stohcastic_thermodynamics/non_equilibrium_potential.jl")
include("stohcastic_thermodynamics/path_integral_formation.jl")
include("stohcastic_thermodynamics/driven_stohcastic_systems.jl")

# Physical constants
export k_B
const k_B = 1.380649e-23  # Boltzmann constant (J/K)

# Re-export from submodules
# Onsager Relations
using .LinearOnsager
using .NonLinearOnsager
using .StochasticOnsagerRelations
using .TimeDependentOnsager

# Entropy Production
using .LocalEntropyProduction
using .EntropyProductionRate
using .TimeDependentEntropy
using .EntropyOptimization

# Stochastic Thermodynamics
using .FluctuationTheorems
using .FokkerPlanck
using .LangevinDynamics
using .NoiseInducedBifurications
using .NonEquilibriumPotential
using .PathIntegralFormation
# Note: DrivenStochasticSystems is included but not re-exported to avoid ambiguity
# with StochasticOnsagerRelations which exports the same types

end # module NET

