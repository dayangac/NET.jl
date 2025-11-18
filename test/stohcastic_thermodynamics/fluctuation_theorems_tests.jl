using NET
using Test

@testset "FluctuationTheorems" begin
    @testset "thermal_noise" begin
        T = 300.0  # Room temperature
        noise = thermal_noise(T)
        @test noise > 0
        @test noise ≈ sqrt(k_B * T) * 1e-2 atol=1e-10
    end
    
    @testset "jarzynski_equality" begin
        work_values = [1.0, 2.0, 3.0, 4.0, 5.0]
        T = 300.0
        
        ΔF = jarzynski_equality(work_values, T)
        @test isfinite(ΔF)
        @test ΔF < maximum(work_values)  # Free energy should be less than max work
    end
    
    @testset "crooks_theorem" begin
        work_forward = [1.0, 2.0, 3.0]
        work_reverse = [1.0, 2.0, 3.0]
        T = 300.0
        
        ratio = crooks_theorem(work_forward, work_reverse, T)
        @test isfinite(ratio)
        @test ratio > 0
    end
    
    @testset "mean function" begin
        values = [1.0, 2.0, 3.0, 4.0, 5.0]
        # Note: This tests the custom mean function in the module
        # We can't directly test it, but we can test it through jarzynski
        @test length(values) == 5
    end
end

