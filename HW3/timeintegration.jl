using Gadfly

α = 4
β = 2
exact_solution(t, x) = cospi(β * (x-u*t)) + exp(-α^2*π^2 * ϵ * t) * cospi(α * (x-u*t))

# Merson's method
## Butcher array
s = 5
c = zeros(s)
a = zeros(s, s)
b = zeros(s)

n_t = 100
τ = 1 / (n_t)

n_x = 100
h = 1 / (n_x)
x = linspace(0, 1, n_x)

c[2] = 1/3; a[2, 1] = 1/3;
c[3] = 1/3; a[3, 1] = 1/6; a[3, 2] = 1/6;
c[4] = 1/2; a[4, 1] = 1/8; a[4, 2] =   0; a[4, 3] =  3/8;
c[5] = 1;   a[5, 1] = 1/2; a[5, 2] =   0; a[5, 3] = -3/2; a[5, 4] =   2;
            b[   1] = 1/6; b[   2] =   0; b[   3] =    0; b[   4] = 2/3; b[5] = 1/6;

Pe = 400
ϵ = 1/Pe
u = 1
κ = 1/2
const_c = abs(u)*τ/h
const_d = 2ϵ*τ/h^2

q(0, 0)
q(t, x::Number) = (2π)^2 * ϵ * cospi(2*(x - u*t))
q(t, y::AbstractArray) = map((x) -> q(t, x), y)

# cyclicly define A
A_min_2 = const_c/4 * (1 - κ)
A_min_1 = -const_c/4 * (5 - 3κ) - const_d/2
A_id = const_d + (3 - 3κ)const_c/4
A_plus_1 = (1 + κ)const_c/4 - const_d/2
A = spzeros(n_x, n_x)
for i in 3:n_x-1
    A[i, i - 2] = A_min_2
    A[i, i - 1] = A_min_1
    A[i, i] = A_id
    A[i, i + 1] = A_plus_1
end

# -2 term
A[1, end - 1] = A[2, end] = A[end, end-2] = A_min_2
# -1 term
A[1, end] = A[2, 1] = A[end, end-1] = A_min_1
# 0 term
A[1, 1] = A[2, 2] = A[end, end] = A_id
# +1 term
A[1, 2] = A[2, 3] = A[end, 1] = A_plus_1

f(t, y) = -A*y/τ + q(t, x)

y = zeros(n_t, n_x)

# initial y
y[1, :] = map((x) -> exact_solution(0, x), linspace(0, 1, n_x))

# iterate
for n = 1:(n_t - 1)
  t_n = n * τ

  # determine k's
  k = zeros(s, n_x)
  for i in 1:s
    # f's argument in t
    f_t = t_n + c[i] * τ

    # f's argument in y
    # f_y = y_n + τ * Σa_ij k_j
    f_y = y[n, :]
    for j in 1:s
      if a[i, j] ≠ 0
        f_y += τ * a[i, j] * k[j, :]
      end
    end

    k[i, :] = f(f_t, f_y)
  end

  # perform one step
  y[n+1, :] = y[n, :] + τ * (sum([ b[i] * k[i, :] for i in 1:s ]))
end

# gp
plot_n = n_t
  plot(
    layer(Geom.point, y = y[plot_n, :], x = linspace(0, 1, n_x)),
    layer(Geom.line, y = map((x) -> exact_solution((plot_n-1)*τ, x), linspace(0, 1, n_x)), x = linspace(0, 1, n_x)))
