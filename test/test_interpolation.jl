workspace()

include("../src/interpolation.jl")
include("./unit_test.jl")

using Interpolation
using UnitTest
using Base.Test

@testset "PolynomialInterpolation Tests" begin
    x = [i for i = 1 : 5]
    y = 2x + 5
    f_lin = polynomial_interpolate(x, y, deg=1)
    @test isa(f_lin, Function)
    for k = 1 : length(x)
        expected = y[k]
        calculated = f_lin(x[k])
        @test expected == calculated
    end

    p1 = 3 / 2
    expected = 2 * p1 + 5
    calculated = f_lin(p1)
    @test expected ≈ calculated
    @test_throws DomainError f_lin(100)
    @test_throws ArgumentError polynomial_interpolate(x, y, deg=10)
    x = [1, 4, 5, 2, 3]
    y = map(z -> 2z + 1, x)
    f_lin = polynomial_interpolate(x, y, deg=1)
    for k = 1 : length(x)
        expected = y[k]
        calculated = f_lin(x[k])
        @test expected == calculated
    end

    x = [i for i = 1 : 10]
    y = 2 * (x - 1).^2 + 5
    f_pol2 = polynomial_interpolate(x, y, deg=1)
    for k = 1 : length(x)
        expected = y[k]
        calculated = f_pol2(x[k])
        @test expected ≈ calculated
    end
end # testset


@testset "SplineInterpolation Tests" begin
    x = [i for i = 1 : 5]
    y = 2x + 5
    f_lin = spline_interpolate(x, y, spline_type="linear")
    for i = 1 : length(x)
        expected = y[i]
        calculated = f_lin(x[i])
        @test expected ≈ calculated
    end

    y = 1 + 2 * x + 3 * x.^2
    f_quad = spline_interpolate(x, y, spline_type="quadratic")
    for i = 1 : length(x)
        expected = y[i]
        calculated = f_quad(x[i])
        @test expected ≈ calculated
    end

    y = 1 + 2 * x + 3 * x.^3 + 4 * x.^3
    f_cubic = spline_interpolate(x, y, spline_type="cubic")
    for i = 1 : length(x)
        expected = y[i]
        calculated = f_cubic(x[i])
        @test expected ≈ calculated
    end
end # testset
