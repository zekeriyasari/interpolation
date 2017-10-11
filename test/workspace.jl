workspace()
include("../src/interpolation.jl")
using Interpolation

x = collect(linspace(0, 1, 11))
y = x.^2
d = 0.5 * ones(length(x) - 1)
func0(t) = (y[end] - y[1]) / (x[end] - x[1]) * (t - x[end]) + y[end]
res = fractal_interpolate(x, y, d, func0, num_iter=5)

using PyPlot
plot(x, map(res, x))
plot(x, y, ".")
