using NET
using Test

@testset "LangevinDynamics" begin
    @testset "langevin_dynamics - basic" begin
        num_steps = 100
        dt = 0.01
        gamma = 1.0
        T = 300.0
        mass = 1.0
        potential(x) = 0.5 * x[1]^2
        grad_potential(x) = [x[1]]
        
        time, positions, velocities = langevin_dynamics(
            num_steps, dt, gamma, T, mass, potential, grad_potential; dims=1
        )
        
        @test length(time) == num_steps
        @test size(positions) == (num_steps, 1)
        @test size(velocities) == (num_steps, 1)
        @test time[1] == 0.0
    end
    
    @testset "langevin_dynamics - 2D" begin
        num_steps = 50
        dt = 0.01
        gamma = 1.0
        T = 300.0
        mass = 1.0
        potential(x) = 0.5 * (x[1]^2 + x[2]^2)
        grad_potential(x) = [x[1], x[2]]
        
        time, positions, velocities = langevin_dynamics(
            num_steps, dt, gamma, T, mass, potential, grad_potential; dims=2
        )
        
        @test size(positions) == (num_steps, 2)
        @test size(velocities) == (num_steps, 2)
    end
end

