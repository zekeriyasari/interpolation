workspace()
include("../src/interpolation.jl")

using PyPlot
using Interpolation

x = [i for i = 1 : 10]
y = 2 * (x - 1).^2 + 5
f_pol1 = polynomial_interpolate(x, y, deg=1)
f_pol2 = polynomial_interpolate(x, y, deg=2)
f_pol9 = polynomial_interpolate(x, y, deg=9)
domain = linspace(1, 10, 100)
plot(domain, map(f_pol1, domain), label="1st degree")
plot(domain, map(f_pol2, domain), label="2nd degree")
plot(domain, map(f_pol9, domain), label="9th degree")
stem(x, y, label="Points")
legend()
