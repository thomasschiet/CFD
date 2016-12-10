using Gadfly

# Merson's method
## Butcher array
s = 5
c = zeros(s)
a = zeros(s, s)
b = zeros(s)

n_t = 100
τ = 0.1

n_x = 100
h = 1 / n_x

c[2] = 1/3; a[2, 1] = 1/3;
c[3] = 1/3; a[3, 1] = 1/6; a[3, 2] = 1/6;
c[4] = 1/2; a[4, 1] = 1/8; a[4, 2] =   0; a[4, 3] =  3/8;
c[5] = 1;   a[5, 1] = 1/2; a[5, 2] =   0; a[5, 3] = -2/3; a[5, 4] =   2;
            b[   1] = 1/6; b[   2] =   0; b[   3] =    0; b[   4] = 2/3; b[5] = 1/6;

Pe = 200
ϵ = 1/Pe
u = 1
κ = 1/3
const_c = abs(u)*τ/h
const_d = 2ϵ*τ/h^2

q(t, x::Number) = (2π)^2 * ϵ * cospi(2*(x - u*t))
q(t, y::AbstractArray) = map((x) -> q(t, x), y)

# cyclicly define A
A = spzeros(n_x, n_x)
for i in 3:n_x-1
    A[i, i - 2] = const_c/4 * (1 - κ)
    A[i, i - 1] = -const_c/4 * (5 - 3κ) - const_d/2
    A[i, i] = const_d + (3 - 3κ)const_c/4
    A[i, i + 1] = (1 + κ)const_c/4 - const_d/2
end
# -2 term
A[1, end - 1] = A[2, end] = A[end, end-2] = const_c/4 * (1 - κ)
# -1 term
A[1, end] = A[2, 1] = A[end, end-1] = -const_c/4 * (5 - 3κ) - const_d/2
# 0 term
A[1, 1] = A[2, 2] = A[end, end] = const_d + (3 - 3κ)const_c/4
# +1 term
A[1, 2] = A[2, 3] = A[end, 1] = (1 + κ)const_c/4 - const_d/2

# Af = factorize(A)

f(t, y) = -A*y + q(t, y)

y = zeros(n_t, n_x)
rand(3)
# initial y
y[1, :] = abs(rand(n_x))

# iterate
for n = 1:(n_t - 1)
  t_n = n_t * τ

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
  y[n+1, :] = y[n, :] + τ * sum([ b[i] * k[i, :] for i in 1:s ])
end

plot(y = y[10, :], Geom.line)
