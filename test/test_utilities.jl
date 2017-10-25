include("../src/utilities.jl")
using Base.Test

@testset "UnitStep Test" begin
pairs = Dict(-1 => 0, 0 =>1, 1 => 1)
    for key in keys(pairs)
        expected = pairs[key]
        calculated = unit_step(key)
        @test expected == calculated
    end
end


@testset "UnitPulse Test" begin
    pairs = Dict(-4 => 0, -2 =>1, 0 =>1, 2 => 0, 4 => 0)
    for key in keys(pairs)
        expected = pairs[key]
        calculated = unit_pulse(key, tl=-2, th=2)
        @test expected == calculated
    end
end

@testset "FunctionNorms" begin
    f(t) = t^2
    t0, t1 = 0, 1
    @test fnorm1(f, t0, t1) ≈ 1 / 3
    @test fnormp(f, t0, t1, 2) ≈ (1 / 5)^(1 / 2)
    @test fnormInf(f, t0, t1) ≈ 1
end
