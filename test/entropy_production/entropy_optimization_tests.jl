using NET
using Test

@testset "EntropyOptimization" begin
    @testset "compute_entropy" begin
        J = [1.0 2.0; 3.0 4.0]
        X = [1.0 2.0; 3.0 4.0]
        
        entropy = compute_entropy(J, X)
        @test entropy ≈ 30.0  # 10 + 20 = 30
    end
    
    @testset "entropy_gradient" begin
        J = [1.0 2.0; 3.0 4.0]
        X = [1.0 2.0; 3.0 4.0]
        
        ∇J, ∇X = entropy_gradient(J, X)
        @test ∇J == X
        @test ∇X == J
    end
    
    @testset "gradient_descent" begin
        J = [1.0 2.0; 3.0 4.0]
        X = [1.0 2.0; 3.0 4.0]
        learning_rate = 0.01
        max_iterations = 10
        
        J_opt, X_opt = gradient_descent(J, X, learning_rate, max_iterations)
        
        @test size(J_opt) == size(J)
        @test size(X_opt) == size(X)
        # Entropy should decrease (or at least change)
        @test compute_entropy(J_opt, X_opt) != compute_entropy(J, X)
    end
end

