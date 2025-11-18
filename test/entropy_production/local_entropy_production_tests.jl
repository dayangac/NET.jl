using NET
using Test

@testset "LocalEntropyProduction" begin
    @testset "compute_local_entropy_prozduction" begin
        J = [1.0 2.0; 3.0 4.0]
        X = [1.0 2.0; 3.0 4.0]
        
        σ = compute_local_entropy_prozduction(J, X)
        @test length(σ) == size(J, 2)
        @test σ[1] ≈ 10.0  # 1*1 + 3*3 = 10
        @test σ[2] ≈ 20.0  # 2*2 + 4*4 = 20
    end
    
    @testset "validate_dimensions_entropy" begin
        J = [1.0 2.0; 3.0 4.0]
        X = [1.0 2.0; 3.0 4.0]
        
        @test validate_dimensions_entropy(J, X) === nothing
        
        X_wrong = [1.0; 2.0; 3.0]
        @test_throws DimensionMismatch compute_local_entropy_prozduction(J, X_wrong)
    end
end

