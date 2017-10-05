workspace()

include("../src/interpolation.jl")
include("./unit_test.jl")

using Interpolation
using UnitTest

println("Test: unit_step")
pairs = Dict(-1 => 0, 0 =>1, 1 => 1)
for key in keys(pairs)
    expected = pairs[key]
    calculated = unit_step(key)
    equal(expected, calculated) ? pass() : fail("Expected $expected, got $calculated")
end
println("Passed...\n")


println("Test: unit_pulse")
pairs = Dict(-4 => 0, -2 =>1, 0 =>1, 2 => 0, 4 => 0)
for key in keys(pairs)
    expected = pairs[key]
    calculated = unit_pulse(key, tl=-2, th=2)
    equal(expected, calculated) ? pass() : fail("Expected $expected, got $calculated")
end
println("Passed...\n")


println("Test: polynomial_interpolate")
x = [i for i = 1 : 5]
y = 2x + 5
# Linear interpolation
f_lin = polynomial_interpolate(x, y, deg=1)
# Check f_lin is a Function
isa(f_lin, Function) ? pass() : fail("Expected type $Function, got $(typeof(f_lin))")
# Check if f_lin passes through interpolaton points
for k = 1 : length(x)
    expected = y[k]
    calculated = f_lin(x[k])
    equal(expected, calculated) ? pass() : fail("Expected $expected, got $calculated")
end
# Check if f_lin interpolates an intermediate  value
p1 = 3 / 2
expected = 2 * p1 + 5
calculated = f_lin(p1)
almost_equal(expected, calculated, tol=EPSILON) ? pass() : fail("Expected $expected, got $calculated")
# Check if DomainError is raised when called outside of interpolation domain
raises(f_lin, DomainError, 100) ? pass() : fail("$DomainError does not raised")
# Check if ArgumentError is raised when polynomial degree is greater than number of points
raises(polynomial_interpolate, ArgumentError, x, y, deg=10)
# Check if interpolation works properly when the points are mixed.
x = [1, 4, 5, 2, 3]
y = map(z -> 2z + 1, x)
f_lin = polynomial_interpolate(x, y, deg=1)
for k = 1 : length(x)
    expected = y[k]
    calculated = f_lin(x[k])
    equal(expected, calculated) ? pass() : fail("Expected $expected, calculated $calculated")
end
# Polynomial interpolate to the data
x = [i for i = 1 : 10]
y = 2 * (x - 1).^2 + 5
f_pol2 = polynomial_interpolate(x, y, deg=1)
for k = 1 : length(x)
    expected = y[k]
    calculated = f_pol2(x[k])
    almost_equal(expected, calculated, tol=1e-3) ? pass() : fail("Expected $expected, got $calculated")
end
# # Plotting the results
# using PyPlot
# f_pol1 = polynomial_interpolate(x, y, deg=1)
# f_pol2 = polynomial_interpolate(x, y, deg=2)
# f_pol9 = polynomial_interpolate(x, y, deg=9)
# domain = linspace(1, 10, 100)
# plot(domain, map(f_pol1, domain), label="1st degree")
# plot(domain, map(f_pol2, domain), label="2nd degree")
# plot(domain, map(f_pol9, domain), label="9th degree")
# stem(x, y, label="Points")
# legend()
# println("Passed...\n")


println("Test: spline_interpolate...")
x = [i for i = 1 : 5]
y = 2x + 5
# Perform linear interpolation.
f_lin = spline_interpolate(x, y, spline_type="linear")
# Check if f_lin passes through the interpolation data
for i = 1 : length(x)
    expected = y[i]
    calculated = f_lin(x[i])
    almost_equal(expected, calculated) ? pass() : fail("Expected $expected, got $calculated instead.")
end
# Perform quadratic interpolation
y = 1 + 2 * x + 3 * x.^2
f_quad = spline_interpolate(x, y, spline_type="quadratic")
# Check if f_quad passes through the interpolation data
for i = 1 : length(x)
    expected = y[i]
    calculated = f_quad(x[i])
    almost_equal(expected, calculated) ? pass() : fail("Expected $expected, got $calculated instead.")
end
# Perform qubic interpolation
y = 1 + 2 * x + 3 * x.^3 + 4 * x.^3
f_cubic = spline_interpolate(x, y, spline_type="cubic")
# Check if f_quad passes through the interpolation data
for i = 1 : length(x)
    expected = y[i]
    calculated = f_cubic(x[i])
    almost_equal(expected, calculated) ? pass() : fail("Expected $expected, got $calculated instead.")
end

# using PyPlot
# println("Test: spline_interpolate...")
# x = [i for i = 1 : 100]
# y = 1 + 2 * x + 3 * x.^3
# f = spline_interpolate(x, y, spline_type="quadratic")
# plot(x, y, ".")
# plot(x, map(f, x))
# println("Passed...\n")

# using PyPlot
# println("Test: spline_interpolate...")
# x = [i for i in  linspace(-4, 4, 100)]
# y = 1 ./ (1 + x.^2) .* sin.(x)
# f = spline_interpolate(x, y, spline_type="cubic")
# plot(x, y, ".", label="Points")
# plot(x, map(f, x), label="Cubic Interpolation")
# legend()
# xlabel("x")
# ylabel("y")
# println("Passed...\n")

# Animated plot
using PyPlot
println("Test: spline_interpolate...")
x = [i for i in  linspace(-4, 4, 100)]
y = 1 ./ (1 + x.^2) .* sin.(x)
f = spline_interpolate(x, y, spline_type="cubic")

fig = figure()
ax = fig[:add_subplot](111)
ax[:plot](x, y, ".")
curve, = ax[:plot]([], [])
for xi in x
    curve[:set_xdata](xi)
    curve[:set_ydata](f(xi))
    fig[:canvas][:draw]()
    fig[:canvas][:flush_events]()
end
println("Passed...\n")
