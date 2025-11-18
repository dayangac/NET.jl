using NET
using Test

@testset "ODESolver" begin
    using NET.ODESolver: rk4_solve, solve_ode, ODESolution
    
    @testset "rk4_solve - exponential decay" begin
        f(u, p, t) = -u
        u0 = 1.0
        tspan = (0.0, 1.0)
        
        t, u = rk4_solve(f, u0, tspan, 0.1)
        @test length(t) > 0
        @test length(u) == length(t)
        @test u[1] ≈ u0
        @test u[end] < u[1]
    end
    
    @testset "rk4_solve - linear growth" begin
        f(u, p, t) = 1.0
        u0 = 0.0
        tspan = (0.0, 2.0)
        
        t, u = rk4_solve(f, u0, tspan, 0.1)
        @test u[end] ≈ 2.0 atol=0.1
    end
    
    @testset "solve_ode" begin
        f(u, p, t) = -u
        u0 = 1.0
        tspan = (0.0, 1.0)
        
        sol = solve_ode(f, u0, tspan; dt=0.01)
        @test sol isa ODESolution
        @test length(sol.t) > 0
        @test length(sol.u) == length(sol.t)
        @test sol.u[1] ≈ u0
    end
    
    @testset "ODESolution struct" begin
        t = [0.0, 1.0, 2.0]
        u = [1.0, 0.5, 0.25]
        sol = ODESolution(t, u)
        @test sol.t == t
        @test sol.u == u
    end
end

