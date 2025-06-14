function [sigma] = custom_power_iteration_noops(Ahat_1, Ahat_2, Bhat_minus, Bhat_plus, Shat, lx1_inv, lx2_inv, x, xi, u, s)
%CUSTOM_POWER_ITERATION Computes max eigenvalue squared of dynamics matrix
%   Input variables are problem variables obtained after L2-hypersphere
%   preconditioning
arguments
    Ahat_1
    Ahat_2
    Bhat_minus
    Bhat_plus
    Shat
    lx1_inv
    lx2_inv
    x
    xi
    u
    s
end

tol_abs = 1e-3;
tol_rel = 1e-3;
eps_buff = 0.1;
j_max = 5000;

N = size(x, 2);

sigma = s ^ 2 + sum(vecnorm(x, 2, 1) .^ 2 + vecnorm(xi, 2, 1) .^ 2 + vecnorm(u, 2, 1) .^ 2);
sigma = sqrt(sigma);

for j = 1 : j_max
    w_k = 1 / sigma * (pagemtimes(Ahat_1, reshape(x(:, 1 : (N - 1)), [], 1, N - 1)) ...
                     + pagemtimes(Ahat_2, reshape(xi(:, 1 : (N - 1)), [], 1, N - 1)) ...
                     + pagemtimes(Bhat_minus, reshape(u(:, 1 : (N - 1)), [], 1, N - 1)) ...
                     + pagemtimes(Bhat_plus, reshape(u(:, 2 : N), [], 1, N - 1)) ...
                     + Shat * s ...
                     - lx1_inv * reshape(x(:, 2 : N), [], 1, N - 1) - lx2_inv * reshape(xi(:, 2 : N), [], 1, N - 1));

    x(:, 1) = Ahat_1(:, :, 1)' * w_k(:, :, 1);
    xi(:, 1) = Ahat_2(:, :, 1)' * w_k(:, :, 1);
    u(:, 1) = Bhat_minus(:, :, 1)' * w_k(:, :, 1);
    s(:, 1) = Shat(:, :, 1)' * w_k(:, :, 1);

    x(:, 2 : (N - 1)) = pagemtimes(pagetranspose(Ahat_1(:, :, 2 : (N - 1))), w_k(:, :, 2 : (N - 1))) - lx1_inv * w_k(:, :, 1 : (N - 2));
    xi(:, 2 : (N - 1)) = pagemtimes(pagetranspose(Ahat_2(:, :, 2 : (N - 1))), w_k(:, :, 2 : (N - 1))) - lx2_inv * w_k(:, :, 1 : (N - 2));
    u(:, 2 : (N - 1)) = pagemtimes(pagetranspose(Bhat_minus(:, :, 2 : (N - 1))), w_k(:, :, 2 : (N - 1))) + pagemtimes(pagetranspose(Bhat_plus(:, :, 1 : (N - 2))), w_k(:, :, 1 : (N - 2)));
    s(:, 1) = s + sum(pagemtimes(pagetranspose(Shat(:, :, 2 : (N - 1))), w_k(:, :, 2 : (N - 1))), 3);

    x(:, N) = -lx1_inv * w_k(:, :, N - 1);
    xi(:, N) = -lx2_inv * w_k(:, :, N - 1);
    u(:, N) = Bhat_plus(:, :, N - 1)' * w_k(:, :, N - 1);

    sigma_star = s ^ 2 + sum(vecnorm(x, 2, 1) .^ 2 + vecnorm(xi, 2, 1) .^ 2 + vecnorm(u, 2, 1) .^ 2);
    sigma_star = sqrt(sigma_star);

    if abs(sigma_star - sigma) <= tol_abs + tol_rel * max(sigma_star, sigma)
        break
    elseif j < j_max
        sigma = sigma_star;
    end
end

sigma = (1 + eps_buff) * sigma_star;

end

