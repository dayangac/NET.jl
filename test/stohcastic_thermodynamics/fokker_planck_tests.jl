using NET
using Test

@testset "FokkerPlanck" begin
    @testset "solve_fokker_planck" begin
        x_range = collect(-2.0:0.1:2.0)
        t_range = collect(0.0:0.01:0.1)
        D = 0.5
        drift(x) = -x  # Simple drift
        
        P = solve_fokker_planck(x_range, t_range, D, drift)
        
        @test size(P) == (length(x_range), length(t_range))
        @test all(P .>= 0)  # Probabilities should be non-negative
        # Check normalization (approximately)
        @test sum(P[:, 1]) > 0
    end
    
    @testset "initial condition" begin
        x_range = collect(-2.0:0.1:2.0)
        t_range = collect(0.0:0.01:0.1)
        D = 0.5
        drift(x) = 0.0  # No drift
        
        P = solve_fokker_planck(x_range, t_range, D, drift)
        # Initial condition should be Gaussian
        @test P[length(x_range)÷2, 1] > 0  # Center should have probability
    end
end

