using Gadfly
using DataFrames

include("functions.jl")

a=0.1
b=1
pe=40
scheme="central1"
refinement = true
m1 = 200
m2 = 500

q(x) = -2/π * a * sinpi(x/2) + 4/(π^2 * pe) * a * cospi(x/2) + (1/pe - 1)*b/(1-e)*exp(x)

(ϕ_n, x) = numerical_solution(a=a, b=b, pe=pe, scheme=scheme, refinement=refinement, m1=m1, m2=m2, q = q)
df_n = DataFrame(x=x, y=ϕ_n, label="Numerical")

ϕ_e = exact_solution(x, pe, a, b)
ϕ_e = map(x_1 -> a * cospi(x_1/2) + b / (1-e) * (1 - exp(x_1)), x)
df_e = DataFrame(x = x, y = ϕ_e, label="Exact")

df = vcat(df_n, df_e)

p = plot(df, x="x", y="y", color="label", Geom.line,
         Scale.discrete_color_manual("blue","red"))
