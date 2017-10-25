
include("../src/interpolation.jl")

using PyPlot
import Interpolation

x = [i for i in  linspace(-4, 4, 100)]
y = 1 ./ (1 + x.^2) .* sin.(x)
f = Interpolation.spline_interpolate(x, y, spline_type="cubic")

fig = figure()
ax = fig[:add_subplot](111)
ax[:plot](x, y, ".", label="Interpolation Data")
ax[:set_xlabel]("x")
ax[:set_ylabel]("y")
x_data = []
y_data = []
curve, = ax[:plot]([], [], label="Cubic Spline Interpolation")
ax[:legend]()
ax[:set_title]("Cubic spline interplation to f(x) = 1 / (1 + x^2) sin(x)")
for xi in x
    x_data = [x_data; xi]
    y_data = [y_data; f(xi)]
    curve[:set_xdata](x_data)
    curve[:set_ydata](y_data)
    fig[:canvas][:flush_events]()
    ax[:relim]()
    sleep(0.25)
end
