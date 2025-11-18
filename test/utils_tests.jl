using NET
using Test

@testset "Utils - ODE Solver" begin
    using NET.ODESolver: rk4_solve, solve_ode, ODESolution
    
    # Test simple exponential decay: du/dt = -u
    f_exp(u, p, t) = -u
    u0 = 1.0
    tspan = (0.0, 1.0)
    
    t, u = rk4_solve(f_exp, u0, tspan, 0.1)
    @test length(t) > 0
    @test length(u) == length(t)
    @test u[1] ≈ u0
    @test u[end] < u[1]  # Should decay
    
    # Test ODESolution struct
    sol = solve_ode(f_exp, u0, tspan; dt=0.01)
    @test sol isa ODESolution
    @test length(sol.t) > 0
    @test length(sol.u) == length(sol.t)
    
    # Test linear ODE: du/dt = 1
    f_linear(u, p, t) = 1.0
    sol_linear = solve_ode(f_linear, 0.0, (0.0, 2.0); dt=0.1)
    @test sol_linear.u[end] ≈ 2.0 atol=0.1
end

@testset "Utils - Gradient" begin
    using NET.Gradient: gradient
    
    # Test gradient of quadratic function: f(x) = x[1]^2 + 2*x[2]^2
    f_quad(x) = x[1]^2 + 2*x[2]^2
    x_test = [1.0, 2.0]
    grad = gradient(f_quad, x_test)
    
    @test length(grad) == length(x_test)
    @test grad[1] ≈ 2.0 atol=0.01  # df/dx1 = 2*x[1] = 2
    @test grad[2] ≈ 8.0 atol=0.01  # df/dx2 = 4*x[2] = 8
    
    # Test gradient of linear function: f(x) = 3*x[1] + 4*x[2]
    f_linear(x) = 3*x[1] + 4*x[2]
    x_linear = [1.0, 1.0]
    grad_linear = gradient(f_linear, x_linear)
    @test grad_linear[1] ≈ 3.0 atol=0.01
    @test grad_linear[2] ≈ 4.0 atol=0.01
end

