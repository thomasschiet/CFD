%% settings
n_r = 10;
n_theta = 10;
U_infty = 1;
R_1 = 1;
R_2 = 5;
Delta_eta_1 = 1/(n_theta-1);
Delta_eta_2 = 1/(n_r-1);

%% 
A = sparse(n_r, n_theta);

d_r = (R_2-R_1) / n_r;
d_theta = pi / n_theta;

ind = @(i, j) (i-1)*n_r + j;

for i = 1:n_r
    for j = 1:n_theta
        row = ind(i, j);
        if i ~= 1 && i ~= n_r && j ~= 1 && j ~= n_theta
            eta_1 = (j - 1) / (n_theta - 1);
            eta_2 = (i - 1) / (n_r - 1);
            lambda = (2*Delta_eta_2/eta_2/Delta_eta_1 + (eta_2^2+(eta_2 + 1/(n_theta-1))^2)/Delta_eta_2)^-1;
            a = lambda * eta_2 * Delta_eta_2 / Delta_eta_1;
            b = b;
            c = lambda * (eta_2^2 + (eta_2 + 1/(n_theta-1))^2);
            d = lambda * (eta_2^2 + (eta_2 - 1/(n_theta-1))^2);

            A(row, row + 1) = a;
            A(row, row - 1) = b;
            A(row, row + n_r) = c;
            A(row, row - n_r) = d;
        end
        
        if i > 1 && i < n_r && j ~= 1 && j ~= n_theta
            lambda = 
        end
    end
end

