using Gadfly
using Cairo
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

Lhat(1, 1/3, 100)
# go
for θ = [0, 0.5, 1, 1.5]
  f = (x) -> Lhat(θ, 1/3, x)

  xs = linspace(0, 1e5, 1000)
  ys = map(f, xs)

  plt = plot(
    layer(y = real(ys), x = xs, Geom.line),
    layer(y = imag(ys), x = xs, Geom.line, Theme(default_color = colorant"green")),
    Guide.ylabel("L_h"),
    Guide.xlabel("Pe"),
    Guide.title(string("L_h, Θ = ", θ, "π, κ = 1/3")),
    Guide.manual_color_key("Legend", ["Real(L_h)", "Imag(L_h)"], ["green", "deepskyblue"])
  )

  draw(PNG(string("LvsPE", round(Int, 10θ), ".pdf"), 4inch, 3inch), plt)
# Lhat(1, 1/3, 100)
  # go
    # θ = 1
    f = (x) -> Lhat(θ, x, 10000)

    xs = linspace(0, 1, 1000)
    ys = map(f, xs)

    plt = plot(
      layer(y = real(ys), x = xs, Geom.line),
      layer(y = imag(ys), x = xs, Geom.line, Theme(default_color = colorant"green")),
      Guide.ylabel("L_h"),
      Guide.xlabel("κ"),
      Guide.title(string("L_h, Θ = ", θ, "π, Pe = 10000")),
      Guide.manual_color_key("Legend", ["Real(L_h)", "Imag(L_h)"], ["green", "deepskyblue"])
    )
  draw(PNG(string("LvsKappa", round(Int, 10θ), ".pdf"), 4inch, 3inch), plt)
end
