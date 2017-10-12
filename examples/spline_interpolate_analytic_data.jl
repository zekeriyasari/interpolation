workspace()
include("../src/interpolation.jl")

using PyPlot
using Interpolation

x = [i for i in  linspace(-4, 4, 100)]
y = 1 ./ (1 + x.^2) .* sin.(x)
f = spline_interpolate(x, y, spline_type="cubic")
plot(x, y, ".", label="Points")
plot(x, map(f, x), label="Cubic Interpolation")
legend()
xlabel("x")
ylabel("y")
