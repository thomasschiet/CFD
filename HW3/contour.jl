using Gadfly
using Contour

R(λτ) = 1 + λτ + 1/2*(λτ)^2 + 1/6 * (λτ)^3 + 1/24 * (λτ)^4 + 1/144 * (λτ)^5

xs = linspace(-10, 10, 150)
ys = linspace(-10, 10, 150)
zs = [abs(R(x + im*y)) for x in xs, y in ys]
c = contour(xs, ys, zs, 1.0)
plot_x, plot_y = coordinates(lines(c)[1])
plot(x = plot_x, y = plot_y, Geom.line)
