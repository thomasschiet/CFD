using Gadfly
using Contour

R(λτ) = 1 + λτ + 1/2*(λτ)^2 + 1/6 * (λτ)^3 + 1/24 * (λτ)^4 + 1/144 * (λτ)^5

xs = linspace(-4, 1, 1500)
ys = linspace(-4, 4, 1500)
zs = [abs(R(x + im*y)) for x in xs, y in ys]
c = contour(xs, ys, zs, 1)
plot_x, plot_y = coordinates(lines(c)[1])

  plot(x = plot_x, y = plot_y, Geom.point,Theme(
            # panel_fill = colorant"black",
             default_point_size = 0.5pt,
            #  default_color = colorant"red",
             highlight_width = 0pt
       ),
       Guide.ylabel("Im(τλ)"),
       Guide.xlabel("Re(τλ)")
  )
