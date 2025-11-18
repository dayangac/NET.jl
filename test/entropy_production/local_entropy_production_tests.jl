using NET
using Test

@testset "LocalEntropyProduction" begin
    @testset "compute_local_entropy_prozduction" begin
        J = [1.0 2.0; 3.0 4.0]
        X = [1.0 2.0; 3.0 4.0]
        
        σ = NET.compute_local_entropy_prozduction(J, X)
        @test length(σ) == size(J, 2)
        @test σ[1] ≈ 10.0  # 1*1 + 3*3 = 10
        @test σ[2] ≈ 20.0  # 2*2 + 4*4 = 20
    end
    
    @testset "validate_dimensions_entropy" begin
        J = [1.0 2.0; 3.0 4.0]
        X = [1.0 2.0; 3.0 4.0]
        
        @test NET.validate_dimensions_entropy(J, X) === nothing
        
        # Wrong dimensions (different size)
        X_wrong = [1.0 2.0 3.0; 4.0 5.0 6.0]  # 2x3 instead of 2x2
        @test_throws DimensionMismatch NET.compute_local_entropy_prozduction(J, X_wrong)
    end
end

