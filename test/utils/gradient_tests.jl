using NET
using Test

@testset "Gradient" begin
    using NET.Gradient: gradient
    
    @testset "gradient - quadratic function" begin
        f(x) = x[1]^2 + 2*x[2]^2
        x = [1.0, 2.0]
        grad = gradient(f, x)
        
        @test length(grad) == length(x)
        @test grad[1] ≈ 2.0 atol=0.01
        @test grad[2] ≈ 8.0 atol=0.01
    end
    
    @testset "gradient - linear function" begin
        f(x) = 3*x[1] + 4*x[2]
        x = [1.0, 1.0]
        grad = gradient(f, x)
        
        @test grad[1] ≈ 3.0 atol=0.01
        @test grad[2] ≈ 4.0 atol=0.01
    end
    
    @testset "gradient - custom step size" begin
        f(x) = x[1]^2
        x = [2.0]
        grad1 = gradient(f, x; h=1e-6)
        grad2 = gradient(f, x; h=1e-4)
        
        @test abs(grad1[1] - 4.0) < abs(grad2[1] - 4.0)  # Smaller h should be more accurate
    end
end

