# workspace()

using Plots

theme(:dark)
x = collect(1:10)
y = x.^2
plt = plot(x, y)
display(plt)
gui()
