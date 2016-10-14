using Gadfly

function numerical_solution(;
    J::Union{Int, Void} = nothing,
    h::Union{Number, Void} = nothing,
    a::Number = 0,
    b::Number = 1,
    boundarytype::AbstractString = "Dirichlet",
    q::Union{Vector, Number} = 0,
    u::Number = 1,
    ϵ::Number = 1e-4
  )

  @assert ((typeof(J) == Void) $ (typeof(h) == Void)) || h ≈ 1/J "J or h should be provided, but not both"
  @assert boundarytype == "Dirichlet" || boundarytype == "Neumann" "Boundary type should be either Dirichlet or Neumann"

  if typeof(J) ≠ Void
    # J given, so determine h
    h = 1/J
  else
    # h given, so determine J
    J = round(Int, 1/h)
  end

  if typeof(q) <: Number
    q = fill(q, J)
    println(q)
  end

  if boundarytype == "Neumann"
    γ_0 = (u + 2ϵ/h)a
    γ_1 = (u*h/2 - ϵ)b

    β_0J = u
  else # Dirichlet
    γ_0 = (u + 2ϵ/h)a
    γ_1 = (u - 2ϵ/h)b

    β_0J = 2ϵ/h
  end

  q_tilde = [h*q[1] + γ_0; h*q[2:end-1]; h*q[end] - γ_1]

  β_11 = -2ϵ/h

  β_0 = zeros(J)
  β_0[1:J-1] = u/2 + ϵ/h
  β_0[J] = β_0J

  β_1 = zeros(J)
  β_1[1] = β_11
  β_1[2:J] = u/2 - ϵ/h

  α_minus1 = [0; -β_0[1:end-1]]
  α_0 = β_0 - β_1
  α_1 = [β_0[2:end]; 0]

  A = diagm(α_minus1[2:end], -1) + diagm(α_0, 0) + diagm(α_1[1:end-1], 1)
  b = q_tilde

  Pe_m = u*h/(ϵ)
  println("Mesh Peclet = ", Pe_m)

  return b\A
end
J = 100
φ = numerical_solution(
  u = 1,
  ϵ = 1,
  J = J,
  q = 0,
  a = 0,
  b = 1,
  boundarytype = "Dirichlet")
plot(x = 1:J, y = φ)
