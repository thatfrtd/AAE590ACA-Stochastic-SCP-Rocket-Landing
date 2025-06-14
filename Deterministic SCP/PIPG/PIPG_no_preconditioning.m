function [z_star, w_star, sol_info] = PIPG_no_preconditioning(Q, q, H, h, D, lambda, sigma, rho, tol_abs, tol_rel, j_check, j_restart, j_max, z_ref, w_ref)
%PIPG Proportional Integral Projected Gradient convex conic solver
%   Extrapolated Proportional Integral Projected Gradient (xPIPG)
% rho in [1.5, 1.9] is usually a good choice - should converge about twice
% as fast as not doing extrapolation (rho = 1)
% sigma should be >= norm(H)
% lambda should be >= norm(Q)

sol_info = [];
check_infeasibility = true;
sol_info.primal_infeasible = false;
sol_info.dual_infeasible = false;
sol_info.solution_status = "Unsolved";

tol_infeas = tol_abs;

N = size(q, 2);

t1 = tic;

%% Initialize primal variables
zeta_j = z_ref;

%% Intialize dual variable
eta_j = w_ref * 0 + 1;
w_j = w_ref * 0 + 1;

%% Calculate step-size (sigma from power iteration)
alpha = 2 / (lambda + sqrt(lambda^2 + 4 * sigma));

%% Initialize
z_j = z_ref;

z_jp1 = zeros(size(z_j));

for j = 1 : j_max
    %% Projected gradient step
    z_jp1 = D.project_onto_cones(1 : N, (zeta_j - alpha * (Q * zeta_j + q + H' * eta_j)));
 
    aff_violation(j) = norm(H * z_jp1 - h);
    zhis(:, j) = z_j;

    %% Proportional-integral feedback of affine equality constraint violation
    w_jp1 = eta_j + alpha * (H * (2 * z_jp1 - zeta_j) - h);

    %% Extrapolate transformed primal variables
    zeta_jp1 = (1 - rho) * zeta_j + rho * z_jp1;

    %% Extrapolate transformed dual variable
    eta_jp1 = (1 - rho) * eta_j + rho * w_jp1;

    %% Check stopping criterion every j_check iterations
    if mod(j, j_check) == 0
        [terminate, zhat_inf_dj, what_inf_dj] = stopping(z_jp1, w_jp1, z_j, w_j, tol_abs, tol_rel);
        sol_info.zhat_inf_dj(j / j_check) = zhat_inf_dj;
        sol_info.what_inf_dj(j / j_check) = what_inf_dj;

       if terminate % stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel)
           sol_info.iterations = j;
           sol_info.solution_status = "Optimal";
           %fprintf("Opt at %d\n", j)
           break;
       elseif j == j_max && check_infeasibility % Try to check infeasibility - should be test which can be applied before the end...
           if vecnorm(w_jp1 - w_j, 2, 1) > tol_infeas
               sol_info.primal_infeasible = true;
               sol_info.iterations = j;
               sol_info.solution_status = "Infeasible";
              %fprintf("Prim infeas at %d\n", j)
           end
           if vecnorm(z_jp1 - z_j, 2, 1) > tol_infeas 
               sol_info.dual_infeasible = true;
               sol_info.iterations = j;
               sol_info.solution_status = "Infeasible";
               %fprintf("Dual infeas at %d\n", j)
           end
       end
    end

    %% Update values for next iteration
    z_j = z_jp1;

    w_j = w_jp1;

    zeta_j = zeta_jp1;

    eta_j = eta_jp1;

    %% Restart
    if mod(j, j_restart) == 0
        zeta_j = ones(size(zeta_j));
        eta_j = ones(size(eta_j));
    end
end

%% Recover original primal variables
z_star = z_jp1;

%% Retain transformed dual variable
w_star = w_jp1;

t2 = toc(t1);

fprintf("PIPG Time: %.3f ms\n", t2 * 1000 )


figure
plot(aff_violation)
yscale("log")


figure
zhis =  zhis;
scatter(zhis(1, :), zhis(2, :))
axis equal
grid on

figure
zhis = zhis;
plot(vecnorm(zhis - zhis(:, end)))
yscale("log")
grid on


end

