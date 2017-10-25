# workspace()
include("../src/interpolation.jl")

# using PyPlot
using Plots
using PlotThemes
using Interpolation

x = [i for i in  linspace(-4, 4, 100)]
y = 1 ./ (1 + x.^2) .* sin.(x)
f = spline_interpolate(x, y, spline_type="cubic")

# theme(:dark)
# plot(x, y)
plt1 = scatter(x, y)
plot!(x, map(f, x), xlabel="x", ylabel="y")
display(plt1)
