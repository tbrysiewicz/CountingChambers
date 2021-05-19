using Test
using CountingChambers

@testset "Tests" begin
@testset "Correctness-Tests" begin
    @testset "Resonance" begin
        ResBetti=[betti_numbers(resonance_hyperplanes(n); symmetry_group=symmetry_resonance(n)) for n in 2:6]
        @test ResBetti==[[1, 3, 2], [1, 7, 15, 9], [1, 15, 80, 170, 104], [1, 31, 375, 2130, 5270, 3485], [1, 63, 1652, 22435, 159460, 510524, 371909]]
    end
    @testset "Threshold" begin
        ThreshBetti=[betti_numbers(threshold_hyperplanes(n); symmetry_group=symmetry_threshold(n)) for n in 2:6]
        @test ThreshBetti==[[1, 4, 6, 3], [1, 8, 28, 44, 23], [1, 16, 120, 460, 820, 465], [1, 32, 496, 4240, 19660, 43014, 27129], [1, 64, 2016, 36848, 400400, 2453248, 7111650, 5023907]]
    end

end
end
