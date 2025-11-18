using NET
using Test

@testset "NET.jl" begin
    @testset "Utils" begin
        include("utils/ode_solver_tests.jl")
        include("utils/gradient_tests.jl")
    end

    @testset "Onsager Relations" begin
        include("onsager/linear_onsager_tests.jl")
        include("onsager/non_linear_onsager_tests.jl")
        include("onsager/stochastic_onsager_relations_tests.jl")
        include("onsager/time_dependent_onsager_tests.jl")
    end

    @testset "Entropy Production" begin
        include("entropy_production/local_entropy_production_tests.jl")
        include("entropy_production/entropy_production_rate_tests.jl")
        include("entropy_production/time_dependent_entropy_tests.jl")
        include("entropy_production/entropy_optimization_tests.jl")
    end

    @testset "Stochastic Thermodynamics" begin
        include("stohcastic_thermodynamics/fluctuation_theorems_tests.jl")
        include("stohcastic_thermodynamics/fokker_planck_tests.jl")
        include("stohcastic_thermodynamics/langevin_dynamics_tests.jl")
        include("stohcastic_thermodynamics/noise_induced_bifurications_tests.jl")
        include("stohcastic_thermodynamics/non_equilibrium_potential_tests.jl")
        include("stohcastic_thermodynamics/path_integral_formation_tests.jl")
        include("stohcastic_thermodynamics/driven_stohcastic_systems_tests.jl")
    end
end
