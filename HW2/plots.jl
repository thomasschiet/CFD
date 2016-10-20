using Gadfly
using DataFrames

include("functions.jl")

a=0.1
b=1
pe=10
scheme="central1"
refinement = false
m1 = 100
m2 = 100

q(x) = 10000cos(20π*x)

(ϕ_n, x) = numerical_solution(a=a, b=b, pe=pe, scheme=scheme, refinement=refinement, m1=m1, m2=m2, q = q)
df_n = DataFrame(x=x, y=ϕ_n, label="Numerical")

ϕ_e = exact_solution(x, pe, a, b)
df_e = DataFrame(x = x, y = ϕ_e, label="Exact")

df = vcat(df_n)

p = plot(df, x="x", y="y", color="label", Geom.line,
         Scale.discrete_color_manual("blue","red"))
