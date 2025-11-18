using NET
using Test

@testset "TimeDependentEntropy" begin
    @testset "analyze_time_dependent_entropy_production" begin
        J_func(t) = [1.0 * t; 2.0 * t]
        X_func(t) = [1.0; 1.0]
        t_span = (0.0, 1.0)
        
        sol = analyze_time_dependent_entropy_production(J_func, X_func, t_span)
        
        @test sol isa NET.ODESolver.ODESolution
        @test length(sol.t) > 0
        @test length(sol.u) == length(sol.t)
        @test sol.u[1] ≈ 0.0  # Initial entropy production
        @test sol.u[end] > 0.0  # Should increase
    end
    
    @testset "constant fluxes and forces" begin
        J_func(t) = [1.0; 1.0]
        X_func(t) = [1.0; 1.0]
        t_span = (0.0, 1.0)
        
        sol = analyze_time_dependent_entropy_production(J_func, X_func, t_span)
        @test sol.u[end] ≈ 2.0 atol=0.1  # Integral of 2 from 0 to 1
    end
end

