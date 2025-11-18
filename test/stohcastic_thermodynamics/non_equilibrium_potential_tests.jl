using NET
using Test

@testset "NonEquilibriumPotential" begin
    @testset "calculate_non_equilibrium_potential" begin
        grid = collect(0.0:0.1:1.0)
        drift(x) = -x  # Linear drift
        diffusion = 1.0
        
        potential = NET.calculate_non_equilibrium_potential(grid, drift, diffusion)
        
        @test length(potential) == length(grid)
        @test potential[1] == 0.0  # Initial potential
    end
    
    @testset "calculate_non_equilibrium_potential - constant drift" begin
        grid = collect(0.0:0.1:1.0)
        drift(x) = 1.0
        diffusion = 1.0
        
        potential = NET.calculate_non_equilibrium_potential(grid, drift, diffusion)
        
        @test potential[end] > potential[1]  # Should increase
    end
    
    @testset "calculate_non_equilibrium_potential - zero drift" begin
        grid = collect(0.0:0.1:1.0)
        drift(x) = 0.0
        diffusion = 1.0
        
        potential = NET.calculate_non_equilibrium_potential(grid, drift, diffusion)
        
        @test all(potential .== 0.0)  # Should remain zero
    end
end

