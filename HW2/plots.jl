using Gadfly
using DataFrames

include("functions.jl")

a=0.5
b=1
pe=10
scheme="central1"
refinement = true
m1 = 20
m2 = 50

(ϕ_n, x) = numerical_solution(a=a, b=b, pe=pe, scheme=scheme, refinement=refinement, m1=m1, m2=m2)
df_n = DataFrame(x=x, y=ϕ_n, label="Numerical")

ϕ_e = exact_solution(x, pe, a, b)
df_e = DataFrame(x=x, y=ϕ_e, label="Exact")

df = vcat(df_n, df_e)

p = plot(df, x="x", y="y", color="label", Geom.line,
         Scale.discrete_color_manual("blue","red"))
