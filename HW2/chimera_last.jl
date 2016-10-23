using Gadfly
using DataFrames

include("functions.jl")

K = 20
a = 0
b = 1
pe = 1

J = round(Int, collect(logspace(4, 1, 400)))
L2 = zeros(length(J))
i = 1
for j in J
  (x, ϕ_n) = chimera_solution(a = a,
  b = b,
  h = 0.1/(j-1),
  J = j,
  H = 0.9/(K-1),
  K = K
  )
  ϕ_e = sinpi(x) + a + (b-a) * x
  L2[i] = L2Error(ϕ_n, ϕ_e)
  i = i + 1
end

p = plot(x = 1./J, y = L2,
 Scale.y_log10, Scale.x_log10,
 Geom.point, Geom.line,
 Guide.XLabel("h"), Guide.YLabel("L₂ error")
 )

draw(PNG("HW2/plots/chimera_L2error.png", 6inch, 3inch), p)
