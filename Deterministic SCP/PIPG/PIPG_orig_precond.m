function [z_star, what_star, sol_info] = PIPG_orig_precond(qhat, Hhat, h, D, L, L_inv, lambda, sigma, rho, tol_abs, tol_rel, tol_infeas, j_check, j_restart, j_max, z_ref, what_ref, S_z, c_z)
%PIPG Proportional Integral Projected Gradient convex conic solver
%   Extrapolated Proportional Integral Projected Gradient (xPIPG)
% rho in [1.5, 1.9] is usually a good choice - should converge about twice
% as fast as not doing extrapolation (rho = 1)
% sigma should be >= norm(H)
% lambda should be >= norm(Q)

sol_info = [];
check_infeasibility = true;
sol_info.iterations = j_max;
sol_info.primal_infeasible = false;
sol_info.dual_infeasible = false;
sol_info.solution_status = "Unsolved";

N = size(qhat, 2);

t1 = tic;

%% Intialize variables
vhat_j = what_ref;
what_j = what_ref;

%% Calculate step-size (sigma from power iteration)
%alpha = 2 / (lambda + sqrt(lambda^2 + 4 * sigma));
mu = 1;
alpha = @(j) 2 / (2 * lambda + (j + 1) * mu);

beta = @(j) (j + 1) * mu / (2 * sigma);

%% Initialize
zhat_j = L * z_ref;

zhat_jp1 = zeros(size(zhat_j));

j_r = 1;

zhat_inf_dj2 = 1e5;
what_inf_dj2 = 1e5;

for j = 1 : j_max
    alpha_j = alpha(j_r);
    beta_j = beta(j_r);

    %% Proportional-integral feedback of affine equality constraint violation
    what_jp1 = vhat_j + beta_j * (Hhat * zhat_j - h);

    %% Projected gradient step
    zhat_jp1 = L * (S_z \ (D.project_onto_cones(1 : N, S_z * (L_inv * (zhat_j - alpha_j * (lambda * zhat_j + qhat + Hhat' * what_jp1))) + c_z) - c_z));
 
    aff_violation(j) = norm(Hhat * zhat_jp1 - h);
    zhis(:, j) = zhat_j;

    %%
    vhat_jp1 = what_jp1 + beta_j * Hhat * (zhat_jp1 - zhat_j);


    %% Extrapolate variables
    what_jp1 = (1 - rho) * what_j + rho * what_jp1;
    zhat_jp1 = (1 - rho) * zhat_j + rho * zhat_jp1;
    vhat_jp1 = (1 - rho) * vhat_j + rho * vhat_jp1;
        
    %% Check stopping criterion every j_check iterations
    if mod(j, j_check) == 0
        [terminate, zhat_inf_dj, what_inf_dj] = stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel);
        sol_info.zhat_inf_dj(j / j_check) = zhat_inf_dj;
        sol_info.what_inf_dj(j / j_check) = what_inf_dj;
            
        if terminate % stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel)
            sol_info.iterations = j;
            sol_info.solution_status = "Optimal";
            %fprintf("Opt at %d\n", j)
            break;
        elseif check_infeasibility % Try to check infeasibility - should be test which can be applied before the end...
            terminate = false;
            if what_inf_dj < tol_abs && zhat_inf_dj > tol_infeas
                sol_info.primal_infeasible = true;
                sol_info.iterations = j;
                sol_info.solution_status = "Infeasible";
                terminate = true;
                %fprintf("Prim infeas at %d\n", j)
            end
            if zhat_inf_dj < tol_abs && what_inf_dj > tol_infeas
                sol_info.dual_infeasible = true;
                sol_info.iterations = j;
                sol_info.solution_status = "Infeasible";
                terminate = true;
                %fprintf("Dual infeas at %d\n", j)
            end
            if terminate
                %break;
            end
        end
    end

    %% Restart
    [~, zhat_inf_djp1, what_inf_djp1] = stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel);
%    if (zhat_inf_djp1 - zhat_inf_dj2  <= 0) %(vhat_jp1 - zhat_jp1)' * (zhat_jp1 - zhat_j) <= 0 && j_restart == -1 || mod(j, j_restart) ~= 0
    if ([vhat_jp1; zeros([numel(zhat_j) - numel(vhat_jp1), 1])] - zhat_jp1)' * (zhat_jp1 - zhat_j) <= 0 && j_restart == -1 || mod(j, j_restart) ~= 0
        j_r = j_r + 1;
    else
        j_r = 1;
        fprintf("Restarting at j = %d\n", j)
    end
    zhat_inf_dj2 = zhat_inf_djp1;
    what_inf_dj2 = what_inf_djp1;

    %% Update values for next iteration
    zhat_j = zhat_jp1;

    what_j = what_jp1;

    vhat_j = vhat_jp1;

end

%% Recover original primal variables
z_star = S_z * L_inv * zhat_jp1 + c_z;

%% Retain transformed dual variable
what_star = what_jp1;

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
% 
% figure
% zhis = zhis;
% plot(vecnorm(zhis - zhis(:, end)))
% yscale("log")
% grid on

sol_info.zhis = zhis;

end

