exact_solution(x, pe, a, b) = a + (b-a)*(exp(x*pe) - 1)/(exp(pe) - 1)

manufactured_exact(x, a, b) = a * cospi(x/2) + b / (1-e) * (1 - exp(x))

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

  if typeof(q) <: Function
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


function chimera_solution(;
  a = 0,
  b = 1.0,
  h = 0.02,
  J = 20,
  H = 0.1,
  K = 7
  )

  pe = 1

  @assert J * h + K * H ≥ 1 "Cells don't overlap"

  q(x) = sinpi(x)

  # fine grid
  x_fine = zeros(J + 1) # Contains coords of cell bounds
  indices = collect(1:J+1) # indices of first m1 cells
  x_fine[indices] = (indices - 1) * h # set the coords of first m1 cells

  x_coarse = zeros(K + 1) # Contains coords of cell bounds
  indices = 1:K+1 # indices of first m1 cells
  x_coarse[indices] = (indices - 1)*H # set the coords of first m1 cells
  x_coarse = reverse(1 - x_coarse)

  # create fine grid matrix
  dx_fine = diff(x_fine) # width of cells
  y_fine = (x_fine[1:J] + x_fine[2:J+1])/2 # cell centers
  dy_fine = [dx_fine[1]/2; diff(y_fine); dx_fine[end]/2]

  dx_coarse = diff(x_coarse) # width of cells
  y_coarse = (x_coarse[1:K] + x_coarse[2:K+1])/2 # cell centers
  dy_coarse = [dx_coarse[1]/2; diff(y_coarse); dx_coarse[end]/2]

  β_diffusion = 1./(pe*dy_fine)
  β_0 = 0 + β_diffusion[2:J+1]
  β_1 = 0 - β_diffusion[1:J]
  β_0[J] = β_diffusion[J+1]
  β_1[1] = - β_diffusion[1]
  γ_0 = (0+(2/pe)/dx_fine[1])*a
  γ_1 = 0

  if typeof(q) <: Function
    f_1 = map(x -> x[2] * q(x[1]), zip(y_fine, dx_fine))
  elseif typeof(q) <: Vector
    f_1 = q
  elseif typeof(q) <: Number
    f_1 = fill(float(q), J)
  else
    f_1 = zeros(J)
  end

  y_fine_virtual = y_fine[end] + h

  sort!(y_fine)
  sort!(y_coarse)
  y_coarse_left_index  = maximum(find(x -> x < y_fine_virtual, y_coarse))
  y_coarse_right_index = minimum(find(x -> x ≥ y_fine_virtual, y_coarse))

  f_1[1] = f_1[1] + γ_0

  A_11 = full(spdiagm((-β_0[1:end-1], β_0-β_1, β_1[2:end]), (-1, 0, 1)))
  A_12 = zeros(J, K)
  A_12[J, y_coarse_left_index] = β_1[end] * (y_fine_virtual - y_coarse[y_coarse_left_index])/H
  A_12[J, y_coarse_right_index] = β_1[end] * (1 - (y_fine_virtual - y_coarse[y_coarse_left_index])/H)

  β_diffusion = 1./(pe*dy_coarse)
  β_0 = 0 + β_diffusion[2:K+1]
  β_1 = 0 - β_diffusion[1:K]
  β_0[K] = β_diffusion[K+1]
  β_1[1] = - β_diffusion[1]
  γ_0 = 0
  γ_1 = (0-(2/pe)/dx_coarse[K])*b

  if typeof(q) <: Function
    f_2 = map(x -> x[2] * q(x[1]), zip(y_coarse, dx_coarse))
  elseif typeof(q) <: Vector
    f_2 = q
  elseif typeof(q) <: Number
    f_2 = fill(float(q), K)
  else
    f_2 = zeros(K)
  end

  y_coarse_virtual = y_coarse[1] - H

  sort!(y_fine)
  sort!(y_coarse)
  y_fine_left_index  = maximum(find(x -> x < y_coarse_virtual, y_fine))
  y_fine_right_index = minimum(find(x -> x ≥ y_coarse_virtual, y_fine))

  f_2[end] = f_2[end] - γ_1

  f = [f_1; f_2]

  A_22 = full(spdiagm((-β_0[1:end-1], β_0-β_1, β_1[2:end]), (-1, 0, 1)))
  A_21 = zeros(K, J)
  A_21[1, y_fine_left_index] = -β_0[1] * (y_coarse_virtual - y_fine[y_fine_left_index])/h;
  A_21[1, y_fine_right_index] = -β_0[1] * (1 - (y_coarse_virtual - y_fine[y_fine_left_index])/h);
  A = [A_11 A_12; A_21 A_22]
  return ([y_fine; y_coarse], inv(A) * f)
end
