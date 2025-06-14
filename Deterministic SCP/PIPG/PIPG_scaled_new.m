function [z_star, what_star, sol_info] = PIPG_scaled_new(qhat, Hhat, h, D, L, L_inv, lambda, sigma, omega, rho, tol_abs, tol_rel, tol_infeas, j_check, j_max, z_ref, what_ref, S_z, c_z, t_k, glideslope_angle_max, gimbal_max, T_min, T_max)
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

num = size(qhat, 2);
% 
% nx = 7;
% nu = 2;
% np = 1;
% 
% N = 15;

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

%%
%z_star_j = S_z * (L_inv * (zeta_j - alpha * (lambda * zeta_j + qhat + Hhat' * eta_j))) + c_z; %S_z * L_inv * zhat_j + c_z;
% z_star_j = S_z * L_inv * zeta_j + c_z;
% 
% x = reshape(z_star_j(1 : (N * nx)), nx, N);
% xi = reshape(z_star_j((N * nx + 1) : (N * nx * 2)), nx, N); 
% u_p = reshape(z_star_j((N * nx * 2 + 1) : (N * nx * 2 + N * nu)), nu, N); 
% u = u_p;
% for k = 1 : N
%     u(:, k) = [cos(u_p(2, k)) -sin(u_p(2, k)); sin(u_p(2, k)) cos(u_p(2, k))] * [u_p(1, k); 0];
% end
% u = [u; u_p(1, :)];
% 
% s = z_star_j(N * nx * 2 + N * nu + 1); 
% 
% figure
% plot_3DoF_trajectory(t_k, x, u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);
% 
% figure
% %plot_3DoF_time_histories(t_k, x, u);
% comparison_plot_3DoF_time_histories({t_k, t_k}, {x, xi}, {u, u}, ["x", "\xi"], linestyle = ["-", "--"]);

for j = 1 : j_max
    %% Projected gradient step
    zhat_jp1 = L * (S_z \ (D.project_onto_cones(1 : num, S_z * (L_inv * (zeta_j - alpha * (lambda * zeta_j + qhat + Hhat' * eta_j))) + c_z) - c_z));
 
    aff_violation(j) = norm(Hhat * zhat_jp1 - h);
    zhis(:, j) = zhat_j;

    %% Proportional-integral feedback of affine equality constraint violation
    what_jp1 = eta_j + beta * (Hhat * (2 * zhat_jp1 - zeta_j) - h);

    %% Extrapolate transformed primal variables
    zeta_jp1 = (1 - rho) * zeta_j + rho * zhat_jp1;

    %% Extrapolate transformed dual variable
    eta_jp1 = (1 - rho) * eta_j + rho * what_jp1;

    %% Check stopping criterion every j_check iterations
    if mod(j, j_check) == 0
        [terminate, zhat_inf_dj, what_inf_dj] = stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel);
        sol_info.zhat_inf_dj(j / j_check) = zhat_inf_dj;
        sol_info.what_inf_dj(j / j_check) = what_inf_dj;

       if terminate % stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel)
           sol_info.iterations = j;
           sol_info.solution_status = "Optimal";
           %fprintf("Opt at %d\n", j)
           %break;
       elseif check_infeasibility % Try to check infeasibility - should be test which can be applied before end...
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

    %% Update values for next iteration
    zhat_j = zhat_jp1;

    what_j = what_jp1;

    zeta_j = zeta_jp1;

    eta_j = eta_jp1;
    fprintf("Iter: %d\n", j)

    % %%
    % z_star_j = S_z * L_inv * zhat_jp1 + c_z;
    % 
    % x = reshape(z_star_j(1 : (N * nx)), nx, N);
    % xi = reshape(z_star_j((N * nx + 1) : (N * nx * 2)), nx, N); 
    % u_p = reshape(z_star_j((N * nx * 2 + 1) : (N * nx * 2 + N * nu)), nu, N); 
    % u = u_p;
    % for k = 1 : N
    %     u(:, k) = [cos(u_p(2, k)) -sin(u_p(2, k)); sin(u_p(2, k)) cos(u_p(2, k))] * [u_p(1, k); 0];
    % end
    % u = [u; u_p(1, :)];
    % 
    % s = z_star_j(N * nx * 2 + N * nu + 1); 
    % 
    % figure
    % plot_3DoF_trajectory(t_k, x, u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);
    % 
    % figure
    % %plot_3DoF_time_histories(t_k, x, u);
    % comparison_plot_3DoF_time_histories({t_k, t_k}, {x, xi}, {u, u}, ["x", "\xi"], linestyle = ["-", "--"]);
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
zhis = L_inv * zhis;
scatter(zhis(25 * 7 + 1, :), zhis(25 * 7 + 2, :))
axis equal
grid on
% 
% figure
% zhis = L_inv * zhis;
% plot(vecnorm(zhis(1:2, :) - zhis(1:2, end)))
% yscale("log")
% grid on
% 

end

