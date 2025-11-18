using NET
using Test

@testset "TimeDependentOnsager" begin
    @testset "TimeDependentOnsagerMatrix construction" begin
        L = zeros(2, 2, 3)
        L[:, :, 1] = [1.0 0.5; 0.5 2.0]
        matrix = NET.TimeDependentOnsagerMatrix(L)
        @test matrix.L == L
    end
    
    @testset "compute_time_dependent_fluxes" begin
        L = zeros(2, 2, 2)
        L[:, :, 1] = [1.0 0.0; 0.0 1.0]
        L[:, :, 2] = [2.0 0.0; 0.0 2.0]
        matrix = NET.TimeDependentOnsagerMatrix(L)
        F = [1.0 2.0; 3.0 4.0]
        
        J = NET.compute_time_dependent_fluxes(matrix, F)
        @test size(J) == size(F)
        # At t=1: J = [1.0 0.0; 0.0 1.0] * [1.0; 3.0] = [1.0; 3.0]
        @test J[:, 1] ≈ [1.0, 3.0]
        # At t=2: J = [2.0 0.0; 0.0 2.0] * [2.0; 4.0] = [4.0; 8.0]
        @test J[:, 2] ≈ [4.0, 8.0]
    end
    
    @testset "validate_dimensions_time_dependent" begin
        L = zeros(2, 2, 2)
        L[:, :, 1] = [1.0 0.0; 0.0 1.0]
        L[:, :, 2] = [1.0 0.0; 0.0 1.0]
        matrix = NET.TimeDependentOnsagerMatrix(L)
        F = [1.0 2.0; 3.0 4.0]
        
        @test NET.validate_dimensions_time_dependent(matrix, F) === nothing
    end
end

