using Gadfly
using DataFrames

include("functions.jl")

a=0.1
b=1
pe = 40
scheme="central1"
refinement = true
m1 = 20
m2 = 40

q(x) = -π/2 * a * sinpi(x/2) + π^2/(4 * pe) * a * cospi(x/2) + (1/pe - 1)*b/(1-e)*exp(x)


# (ϕ_n, x) = numerical_solution(a=a, b=b, pe=pe, scheme=scheme, refinement=refinement, m1=m1, m2=m2, q = q)
# df_n = DataFrame(x=x, y=ϕ_n, label="Numerical")
#
# ϕ_e = exact_solution(x, pe, a, b)
# df_e = DataFrame(x = x, y = ϕ_e, label="Exact")
#
# df = vcat(df_n)
#
# p = plot(df, x="x", y="y", color="label", Geom.line,
#          Scale.discrete_color_manual("blue","red"))





ϕ_e = map(x_1 -> a * cospi(x_1/2) + b / (1-e) * (1 - exp(x_1)), x)
df_e = DataFrame(x = x, y = ϕ_e, label="Exact")

<<<<<<< HEAD
df = vcat(df_n, df_e)
=======
# Determine the L2Error
>>>>>>> 4503243ca5fc2f61e52e66a5e605c9490514d69c

# J = round(Int, collect(logspace(3, 1, 100)))
# L2 = zeros(length(J))
# i = 1
# for j in J
#   (ϕ_n, x) = numerical_solution(a=a, b=b, pe=pe, scheme=scheme, refinement=refinement, m1=j, m2=j, q = q)
#   ϕ_e = exact_solution(x, pe, a, b)
#   L2[i] = L2Error(ϕ_n, ϕ_e)
#   i = i + 1
# end
#
#
# plot(x = 1./J, y = L2)
