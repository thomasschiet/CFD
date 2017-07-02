using Gadfly
using Cairo
using Contour
using Fontconfig

u = 1
h = 1
τ = 1

c = abs(u)*τ/h
d(Pe) = 2τ/(Pe*h^2)

s(θ) = sinpi(θ/2)^2

γ_1(θ, κ) = 2(1-κ)c*s(θ)^2
γ_2(θ, κ) = c*((1-κ)s(θ) + 1) * sinpi(θ)

δ(θ, Pe) = 2d(Pe)*s(θ)

Chat(θ, κ) = γ_1(θ, κ) + im * γ_2(θ, κ)
Dhat(θ, Pe) = δ(θ, Pe)
Lhat(θ, κ, Pe) = Chat(θ, κ) + Dhat(θ, Pe)

layers = Any[]
for κ in 0:0.1:1
  xs = linspace(-4, 1, 1500)
  ys = linspace(-4, 4, 1500)
  zs = [abs(Lhat(x + im*y, κ, 100000)) for x in xs, y in ys]
  cont = contour(xs, ys, zs, 1)
  plot_x, plot_y = coordinates(lines(cont)[1])
  push!(layers, layer(x = plot_x, y = plot_y, Geom.point))
end

R(λτ) = 1 + λτ + 1/2*(λτ)^2 + 1/6 * (λτ)^3 + 1/24 * (λτ)^4 + 1/144 * (λτ)^5

xs = linspace(-4, 1, 1500)
ys = linspace(-4, 4, 1500)
zs = [abs(R(x + im*y)) for x in xs, y in ys]
cont = contour(xs, ys, zs, 1)
plot_x, plot_y = coordinates(lines(cont)[1])
spaceship = layer(x = plot_x, y = plot_y, Geom.point)

plot(layer_0)
plot(layers..., spaceship,
  Theme(
    panel_fill = colorant"white",
    default_point_size = 0.5pt,
    default_color = colorant"black",
    highlight_width = 0pt
  ))
