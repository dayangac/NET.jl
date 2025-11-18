using NET
using Test

@testset "PathIntegralFormation" begin
    @testset "calculate_action" begin
        path = [0.0, 0.1, 0.2, 0.3, 0.4]
        drift(x) = -x
        diffusion = 1.0
        
        action = NET.calculate_action(path, drift, diffusion)
        
        @test action >= 0
        @test isfinite(action)
    end
    
    @testset "calculate_action - linear path" begin
        path = [0.0, 1.0, 2.0, 3.0]
        drift(x) = 0.0
        diffusion = 1.0
        
        action = NET.calculate_action(path, drift, diffusion)
        
        @test action > 0
    end
    
    @testset "calculate_action - constant path" begin
        path = [1.0, 1.0, 1.0, 1.0]
        drift(x) = 0.0
        diffusion = 1.0
        
        action = NET.calculate_action(path, drift, diffusion)
        
        @test action >= 0
    end
end

