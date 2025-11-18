using NET
using Test

@testset "StochasticOnsagerRelations" begin
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
        # J = L*F + noise, with identity L
        @test J[:, 1] ≈ [1.1, 3.3] atol=0.01
        @test J[:, 2] ≈ [2.2, 4.4] atol=0.01
    end
    
    @testset "validate_dimensions_stochastic" begin
        L = zeros(2, 2, 2)
        L[:, :, 1] = [1.0 0.0; 0.0 1.0]
        L[:, :, 2] = [1.0 0.0; 0.0 1.0]
        matrix = StochasticOnsagerMatrix(L)
        F = [1.0 2.0; 3.0 4.0]
        noise = [0.1 0.2; 0.3 0.4]
        
        @test validate_dimensions_stochastic(matrix, F, noise) === nothing
        
        # Wrong noise dimensions
        noise_wrong = [0.1; 0.2]
        @test_throws DimensionMismatch compute_stochastic_fluxes(matrix, F, noise_wrong)
    end
end

