%% settings
n_r = 10;
n_theta = 10;
U_infty = 1;
R_1 = 1;
R_2 = 5;

%% 
A = sparse(n_r, n_theta);

d_r = (R_2-R_1) / n_r;
d_theta = pi / n_theta;

for i = 2:n_r-1
    for j = 2:n_theta-1
        r_i = R_1 + d_r * i;
        a = d_r^2 / (2*d_r^2 + r_i*d_r*d_theta^2 + 2*r_i^2 * d_theta^2);
        b = d_r^2 / (2*d_r^2 + r_i*d_r*d_theta^2 + 2*r_i^2 * d_theta^2);
        c = (r_i*d_r*d_theta^2 + r_i^2*d_theta^2) / (2*d_r^2 + r_i*d_r*d_theta^2 + 2*r_i^2 * d_theta^2);
        d = r_i^2*d_theta^2 / (2*d_r^2 + r_i*d_r*d_theta^2 + 2*r_i^2 * d_theta^2);
        
        row = (i-1)*n_r + j;
        A(row, row + 1) = a;
        A(row, row - 1) = b;
        A(row, row + n_r) = c;
        A(row, row - n_r) = d;
    end
end