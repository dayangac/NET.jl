using NET
using Test

@testset "DrivenStohcasticSystems" begin
    # Note: This file appears to be a duplicate of stochastic_onsager_relations
    # Testing the same functionality
    @testset "StochasticOnsagerMatrix construction" begin
        L = zeros(2, 2, 3)
        L[:, :, 1] = [1.0 0.5; 0.5 2.0]
        matrix = StochasticOnsagerMatrix(L)
        @test matrix.L == L
    end
    
    @testset "compute_stochastic_fluxes" begin
        L = zeros(2, 2, 2)
        L[:, :, 1] = [1.0 0.0; 0.0 1.0]
        L[:, :, 2] = [1.0 0.0; 0.0 1.0]
        matrix = StochasticOnsagerMatrix(L)
        F = [1.0 2.0; 3.0 4.0]
        noise = [0.1 0.2; 0.3 0.4]
        
        J = compute_stochastic_fluxes(matrix, F, noise)
        @test size(J) == size(F)
    end
end

