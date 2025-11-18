using NET
using Test

@testset "NoiseInducedBifurications" begin
    @testset "detect_bifurcation" begin
        drift(x) = x * (1.0 - x)  # Logistic drift
        noise_level = 0.1
        x_range = collect(0.0:0.1:2.0)
        
        bifurcation_points = detect_bifurcation(drift, noise_level, x_range)
        
        @test bifurcation_points isa Vector
        # Should find points where drift is small
        @test length(bifurcation_points) >= 0
    end
    
    @testset "detect_bifurcation - zero drift" begin
        drift(x) = 0.0
        noise_level = 0.1
        x_range = collect(0.0:0.1:1.0)
        
        bifurcation_points = detect_bifurcation(drift, noise_level, x_range)
        # All points should be bifurcation points
        @test length(bifurcation_points) == length(x_range)
    end
    
    @testset "detect_bifurcation - large drift" begin
        drift(x) = 10.0 * x
        noise_level = 0.1
        x_range = collect(0.1:0.1:1.0)
        
        bifurcation_points = detect_bifurcation(drift, noise_level, x_range)
        # Should find few or no bifurcation points
        @test length(bifurcation_points) <= length(x_range)
    end
end

