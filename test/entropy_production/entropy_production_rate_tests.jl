using NET
using Test

@testset "EntropyProductionRate" begin
    @testset "compute_entropy_rate" begin
        J1 = [1.0 2.0; 3.0 4.0]
        X1 = [1.0 2.0; 3.0 4.0]
        J2 = [2.0 3.0; 4.0 5.0]
        X2 = [2.0 3.0; 4.0 5.0]
        Δt = 1.0
        
        rate = NET.compute_entropy_rate(J1, X1, J2, X2, Δt)
        
        # Initial entropy: 10 + 20 = 30
        # Final entropy: 20 + 34 = 54
        # Rate: (54 - 30) / 1.0 = 24
        @test rate ≈ 24.0 atol=0.1
    end
    
    @testset "validate_dimensions_entropy" begin
        J = [1.0 2.0; 3.0 4.0]
        X = [1.0 2.0; 3.0 4.0]
        
        # Use LocalEntropyProduction version to avoid ambiguity
        @test NET.LocalEntropyProduction.validate_dimensions_entropy(J, X) === nothing
    end
end

