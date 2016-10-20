exact_solution(x, pe, a, b) = a + (b-a)*(exp(x*pe) - 1)/(exp(pe) - 1)

function L2Error(ϕ_e, ϕ_n)
  1/sqrt(length(ϕ_e)) * sqrt(sum((ϕ_e-ϕ_n).^2))
end

function numerical_solution(;
  a=0.2,
  b=1.0,
  pe=100,
  scheme="central1",
  refinement = true,
  m1 = 20,
  m2 = 50,
  q = 0
  )

  if refinement
    del = 6/pe
  else
    del = m2/(m1+m2)
  end

  x = collect(zeros(m1+m2+1)) # Contains coords of cell bounds
  j=collect(1:(m1+1)) # indices of first m1 cells
  x[j] = (j-1)*(1-del)/m1 # set the coords of first m1 cells
  j = collect((m1+2):(m1+m2+1))
  x[j] = x[m1+1] + (j-1-m1)*del/m2

  dx = diff(x)
  n = length(x) - 1
  y = (x[1:n] + x[2:n+1])/2
  dy = [dx[1]/2; diff(y); dx[n]/2];

  β_diffusion = 1./(pe*dy)

  if scheme == "central1"
    β_0 = 1/2 + β_diffusion[2:n+1]
    β_1 = 1/2 - β_diffusion[1:n]
    β_0[n] = β_diffusion[n+1]
    β_1[1] = - β_diffusion[1]
    γ_0 = (1+(2/pe)/dx[1])*a
    γ_1 = (1-(2/pe)/dx[n])*b
  end

  if typeof(q) == Function
    f = map(x -> x[2] * q(x[1]), zip(y, dx))
  elseif typeof(q) <: Vector
    f = q
  elseif typeof(q) <: Number
    f = fill(float(q), length(y))
  else
    f = zeros(length(y))
  end

  f[1] = f[1] + γ_0
  f[n] = f[n] - γ_1

  A = spdiagm((-β_0[1:end-1], β_0-β_1, β_1[2:end]), (-1, 0, 1))
  full(A)
  return (inv(full(A))*f, y)
end
