function [z_star, what_star, sol_info] = PIPG_new_explicit_nostop(qhat, Hhat, h, L, L_inv, lambda, sigma, omega, rho, tol_abs, tol_rel, tol_infeas, j_check, j_max, z_ref, what_ref)
%PIPG Proportional Integral Projected Gradient convex conic solver
%   Extrapolated Proportional Integral Projected Gradient (xPIPG)
% rho in [1.5, 1.9] is usually a good choice - should converge about twice
% as fast as not doing extrapolation (rho = 1)
% sigma should be >= norm(H)
% lambda should be >= norm(Q)

sol_info = struct;
sol_info.primal_infeasible = false;
sol_info.dual_infeasible = false;
sol_info.solution_status = "Unsolved";
sol_info.iterations = -1;
sol_info.zhat_inf_dj = zeros([j_max, 1]);
sol_info.what_inf_dj = zeros([j_max, 1]);
sol_info.zhis = zeros([size(z_ref, j_max)]);

N = size(qhat, 2);

t1 = tic;

%% Transform previous primal solution
zhat_ref = L * z_ref;

%% Initialize transformed primal variables
zeta_j = zhat_ref;

%% Intialize transformed dual variable
eta_j = what_ref;
what_j = what_ref;

%% Calculate step-size (sigma from power iteration)
alpha = 2 / (lambda + sqrt(lambda^2 + 4 * omega * sigma));
beta = omega * alpha;

%% Initialize
zhat_j = zhat_ref;

zhat_jp1 = zeros(size(zhat_j));
what_jp1 = zeros(size(what_j));

zhis = zeros([size(zhat_j), j_max]);

for j = 1 : j_max
    %% Projected gradient step
    zhat_jp1 = ball_project(1 : N, L_inv * (zeta_j - alpha * (lambda * zeta_j + qhat + Hhat' * eta_j)), 80/5, 20:1000);
    zhat_jp1 = L * singleton_project(1 : N, zhat_jp1, [1; 2; 3], 10:12);
 
    %aff_violation(j) = norm(Hhat * zhat_jp1 - h);
    zhis(:, j) = zhat_j;

    %% Proportional-integral feedback of affine equality constraint violation
    what_jp1 = eta_j + beta * (Hhat * (2 * zhat_jp1 - zeta_j) - h);

    %% Extrapolate transformed primal variables
    zeta_jp1 = (1 - rho) * zeta_j + rho * zhat_jp1;

    %% Extrapolate transformed dual variable
    eta_jp1 = (1 - rho) * eta_j + rho * what_jp1;
    
    %% Update values for next iteration
    zhat_j = zhat_jp1;
    
    what_j = what_jp1;
    
    zeta_j = zeta_jp1;
    
    eta_j = eta_jp1;
end

%% Recover original primal variables
z_star = L_inv * zhat_jp1;

%% Retain transformed dual variable
what_star = what_jp1;

t2 = toc(t1);
sol_info.time = t2;
% % 
%fprintf("PIPG Time: %.3f ms\n", t2 * 1000 )

% 
% figure
% plot(aff_violation)
% yscale("log")
% 
% 
% figure
% zhis = L_inv * zhis;
% scatter(zhis(25 * 7 + 1, :), zhis(25 * 7 + 2, :))
% axis equal
% grid on
% 
% figure
% zhis = L_inv * zhis;
% plot(vecnorm(zhis(1:2, :) - zhis(1:2, end)))
% yscale("log")
% grid on
% % 

sol_info.zhis = zhis;

end

function x = ball_project(k, z, r, indices)
    %PROJECT Summary of this method goes here
    %   Detailed explanation goes here
    x = z;

    shape = size(z(indices, :));
    z1 = reshape(z(indices, :), [], size(r, 2));

    z1 = r ./ max(vecnorm(z1, 2, 1), r) .* z1;

    x(indices, :) = reshape(z1, shape);
end

function x = singleton_project(k, z, y, indices)
    %PROJECT Summary of this method goes here
    %   Detailed explanation goes here

    x = z;

    x(indices, :) = y;
end