using NET
using Test

@testset "NonLinearOnsager" begin
    @testset "NonLinearOnsagerMatrix construction" begin
        L = zeros(2, 2, 3)
        L[:, :, 1] = [1.0 0.5; 0.5 2.0]
        matrix = NonLinearOnsagerMatrix(L)
        @test matrix.L == L
    end
    
    @testset "NonLinearFluxSystem construction" begin
        L = zeros(2, 2, 2)
        L[:, :, 1] = [1.0 0.0; 0.0 1.0]
        L[:, :, 2] = [1.0 0.0; 0.0 1.0]
        matrix = NonLinearOnsagerMatrix(L)
        F_linear = [1.0 2.0; 3.0 4.0]
        F_nl = [[0.1 0.2; 0.3 0.4]]
        
        system = NonLinearFluxSystem(matrix, F_linear, F_nl)
        @test system.matrix == matrix
        @test system.linear_forces == F_linear
        @test system.nonlinear_forces == F_nl
    end
    
    @testset "compute_total_fluxes" begin
        L = zeros(2, 2, 2)
        L[:, :, 1] = [1.0 0.0; 0.0 1.0]
        L[:, :, 2] = [1.0 0.0; 0.0 1.0]
        matrix = NonLinearOnsagerMatrix(L)
        F_linear = [1.0 2.0; 3.0 4.0]
        F_nl = [[0.1 0.2; 0.3 0.4]]
        system = NonLinearFluxSystem(matrix, F_linear, F_nl)
        
        J = compute_total_fluxes(system)
        @test size(J) == size(F_linear)
        # Linear part: J_linear = F_linear (identity matrix)
        # Nonlinear part: J_nl = F_nl
        # Total: J = F_linear + F_nl
        @test J[:, 1] ≈ [1.1, 3.3] atol=0.01
    end
    
    @testset "validate_dimensions_multinonlinear" begin
        L = zeros(2, 2, 2)
        L[:, :, 1] = [1.0 0.0; 0.0 1.0]
        L[:, :, 2] = [1.0 0.0; 0.0 1.0]
        matrix = NonLinearOnsagerMatrix(L)
        F_linear = [1.0 2.0; 3.0 4.0]
        F_nl = [[0.1 0.2; 0.3 0.4]]
        system = NonLinearFluxSystem(matrix, F_linear, F_nl)
        
        @test validate_dimensions_multinonlinear(system) === nothing
    end
end

