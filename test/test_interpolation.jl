include("../src/interpolation.jl")

import Interpolation
using Base.Test

@testset "PolynomialInterpolation Tests" begin
    x = [i for i = 1 : 5]
    y = 2x + 5
    f_lin = Interpolation.polynomial_interpolate(x, y, deg=1)
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
    @test_throws ArgumentError Interpolation.polynomial_interpolate(x, y, deg=10)
    x = [1, 4, 5, 2, 3]
    y = map(z -> 2z + 1, x)
    f_lin = Interpolation.polynomial_interpolate(x, y, deg=1)
    for k = 1 : length(x)
        expected = y[k]
        calculated = f_lin(x[k])
        @test expected == calculated
    end

    x = [i for i = 1 : 10]
    y = 2 * (x - 1).^2 + 5
    f_pol2 = Interpolation.polynomial_interpolate(x, y, deg=1)
    for k = 1 : length(x)
        expected = y[k]
        calculated = f_pol2(x[k])
        @test expected ≈ calculated
    end
end # testset


@testset "SplineInterpolation Tests" begin
    x = [i for i = 1 : 5]
    y = 2x + 5
    f_lin = Interpolation.spline_interpolate(x, y, spline_type="linear")
    for i = 1 : length(x)
        expected = y[i]
        calculated = f_lin(x[i])
        @test expected ≈ calculated
    end

    y = 1 + 2 * x + 3 * x.^2
    f_quad = Interpolation.spline_interpolate(x, y, spline_type="quadratic")
    for i = 1 : length(x)
        expected = y[i]
        calculated = f_quad(x[i])
        @test expected ≈ calculated
    end

    y = 1 + 2 * x + 3 * x.^3 + 4 * x.^3
    f_cubic = Interpolation.spline_interpolate(x, y, spline_type="cubic")
    for i = 1 : length(x)
        expected = y[i]
        calculated = f_cubic(x[i])
        @test expected ≈ calculated
    end
end # testset

@testset "FractalInterpolation Test" begin
    x = collect(linspace(0, 1, 11))
    y = x.^2
    d = 0.5 * ones(length(x) - 1)
    func0(t) = (y[end] - y[1]) / (x[end] - x[1]) * (t - x[end]) + y[end]
    func = Interpolation.fractal_interpolate(x, y, d, func0, tol=1e-2)
    for i = 1 : length(x)
        @test func(x[i]) ≈ y[i]
    end
end
