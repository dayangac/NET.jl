using NET
using Test

@testset "LinearOnsager" begin
    @testset "OnsagerMatrix construction" begin
        L_2d = [1.0 0.5; 0.5 2.0]
        matrix_2d = OnsagerMatrix(L_2d)
        @test matrix_2d.L == L_2d
        
        L_3d = zeros(2, 2, 3)
        L_3d[:, :, 1] = [1.0 0.5; 0.5 2.0]
        matrix_3d = OnsagerMatrix(L_3d)
        @test matrix_3d.L == L_3d
    end
    
    @testset "compute_fluxes - 2D matrix" begin
        L = OnsagerMatrix([1.0 0.5; 0.5 2.0])
        F = [1.0 2.0; 3.0 4.0]
        
        J = compute_fluxes(L, F)
        @test size(J) == size(F)
        @test J[:, 1] ≈ [2.5, 6.5]  # [1*1 + 0.5*3, 0.5*1 + 2*3]
    end
    
    @testset "compute_fluxes - 3D matrix" begin
        L_3d = zeros(2, 2, 2)
        L_3d[:, :, 1] = [1.0 0.0; 0.0 1.0]
        L_3d[:, :, 2] = [2.0 0.0; 0.0 2.0]
        L = OnsagerMatrix(L_3d)
        F = [1.0 2.0; 3.0 4.0]
        
        J = compute_fluxes(L, F)
        @test size(J) == size(F)
        @test J[:, 1] ≈ [1.0, 3.0]
        @test J[:, 2] ≈ [4.0, 8.0]
    end
    
    @testset "validate_dimensions" begin
        L = OnsagerMatrix([1.0 0.5; 0.5 2.0])
        F = [1.0 2.0; 3.0 4.0]
        
        # Should not throw
        @test validate_dimensions(L, F) === nothing
        
        # Should throw for mismatched dimensions
        F_wrong = [1.0; 2.0; 3.0]
        @test_throws DimensionMismatch compute_fluxes(L, F_wrong)
    end
end

