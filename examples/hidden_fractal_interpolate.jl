# include("../src/interpolation.jl")
# using Interpolation
# using PyPlot
#
# x = [0; 30; 60; 100]
# y = [0; 50; 40; 10]
# H = [0; 30; 60; 100]
# d = [0.3; 0.3; 0.3]
# h = [0.2; 0.2; 0.2]
# l = [-0.1; -0.1; -0.1]
# m = [0.3; 0.; -0.1]
#
# include("../src/fractal_interpolation.jl")
# using FractalInterpolation
# using PyPlot
#
# f0(xi) = [(y[end] - y[1]) / (x[end] - x[1]) * (xi - x[N]) + y[N];
#           (H[end] - H[1]) / (x[end] - x[1]) * (xi - x[N]) + H[end]]
# x = linspace(0, pi, 1000)
# y = sin.(x)
# N = length(x)
# d = ones(N - 1) * 0.5
# plot(x, y, "o",linewidth=5, label="Data set")
# title("Fractal Interpolation")
# legend(loc="upper right")
# xlabel("x")
# ylabel("y")
#
# # f_interpolated = interpolate(f0, x, y, d, 1e-3)
# # domain = linspace(x[1], x[end], 1000)
# # range = f_interpolated.(domain)
# # plot(domain, range , label="Interpolate")
# # legend(loc="upper right")
# f_interpolated = interpolate(x, y, d, 1e-6)
# domain = linspace(x[1], x[end], 1000)
# range = f_interpolated.(domain)
# plot(domain, range , label="Interpolate")
