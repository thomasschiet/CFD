using Gadfly
using DataFrames

include("functions.jl")

J = 100
K = 200
(x, y) = chimera_solution(J = J, K = K, h = 0.5/(J-1), H = 0.5/K)

plot(x = x, y = y)
