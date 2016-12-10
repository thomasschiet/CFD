using Gadfly

# Merson's method
## Butcher array
s = 5
c = zeros(s)
a = zeros(s, s)
b = zeros(s)

n_t = 100
n_x = 100

c[2] = 1/3; a[2, 1] = 1/3;
c[3] = 1/3; a[3, 1] = 1/6; a[3, 2] = 1/6;
c[4] = 1/2; a[4, 1] = 1/8; a[4, 2] =   0; a[4, 3] =  3/8;
c[5] = 1;   a[5, 1] = 1/2; a[5, 2] =   0; a[5, 3] = -2/3; a[5, 4] =   2;
            b[   1] = 1/6; b[   2] =   0; b[   3] =    0; b[   4] = 2/3; b[5] = 1/6;

Pe = 200
ϵ = 1/Pe
u = 1
q(t, x) = (2π)^2 * ϵ * cospi(2*(x - u*t))
f(t, y) = map((x) -> q(t, x), y)


y = zeros(n_t, n_x)
rand(3)
# initial y
y[1, :] = abs(rand(n_x))
τ = 0.1
# iterate
for n = 1:(n_t - 1)
  t_n = n_t * τ

  # determine k's
  k = zeros(s, n_x)
  for i in 1:s
    f_t = t_n + c[i] * τ

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

plot(y = y[3, :], Geom.line)
