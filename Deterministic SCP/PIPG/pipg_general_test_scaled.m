%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA
% AC Optimal Control
% Author: Travis Hastreiter 
% Created On: 16 May, 2025
% Description: Test of parts of PIPG conic convex solver on a general
% problem and compare to CVX. 
% 
% PIPG solves the following optimization problem
%   min_z    1/2 z' Q z + <q, z>
%    s.t.    H z - h = 0
%            z in D
%
% D is a convex cone (which can be composed of many simple cones)
% Q is a positive definite matrix
%
% Most Recent Change: 17 May, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct General Problem
n = 1000;

% Objective
Q = full(sprandsym(n, 0.5, abs(randn(n, 1) * 100)));
q = randn(n, 1) * 100;
q(abs(q) < 80) = 0; % Make q relatively sparse too;

c = 0;

% Initial guess
z_ref = randn(n, 1);

% Affine constraints
z = sym("z", [n, 1]);
g_1 = z(3) + 2 * z(4) - 20;
g_2 = z(5) + z(8) - z(3) - 60;
g_3 = z(3) + z(4);
g_4 = z(3) - z(4) - 1;

g = {g_1, g_2, g_3};

h = zeros(numel(g), 1);
H = zeros(numel(g), n);
for i = 1:numel(g)
    h(i, 1) = -subs(g{i}, z, zeros(size(z)));
    H(i, :) = jacobian(g{i}, z);
end
H = sparse(H);

% H = sprand(n, n, 0.00001); % sparse(H);
% h = zeros(n, 1);
% for i = 1:n
%     if sum(H(i, :)) ~= 0
%         h(i) = randn(1) * 10;
%     end
% end

% Equivalent convex cones for constraints
% u_1 = [sqrt(2) / 2; sqrt(2) / 2];
% zeta_1 = 3;
% u_2 = [sqrt(2) / 2; -sqrt(2) / 2];
% zeta_2 = 3;
% D = Halfspace2(1:2, u_1, u_2, zeta_1, zeta_2);
%%
D = [Ball(1:2, 160), Singleton(10:12, [100;200;300])];

% Convex constraints (for CVX but also to validate satisfaction)
convex_constraints = D.get_constraint_functions;

%%
S_z = diag(ones(size(z_ref)) * 1000);
c_z = -ones(size(z_ref)) * 1000 * 0;
z_guess_hat = S_z \ (z_ref - c_z);

Qtilde = S_z' * Q * S_z;
qtilde = S_z' * (Q' * c_z + q);
Htilde = H * S_z;
htilde = -H * c_z + h;
ctilde = (1/2 * c_z' * Q + q') * c_z + c;

%%
% PIPG hyperparameters
lambda = 1.08 * norm(Qtilde); % lambda should be >= norm(Q) - 1.08 seems optimal?
rho = 1.9;

tol_abs = 1e-5;
tol_rel = 1e-5;
tol_infeas = 1e-3;

j_check = 10;
j_max = 100000;

what_ref = zeros(size(h));

%% Test Customized Preconditioning
t1 = tic;
[qhat, Hhat, L, L_inv] = L2_hypersphere_preconditioning(Qtilde, qtilde, Htilde, lambda);
t2 = toc(t1);
fprintf("L2 Hypersphere Preconditioning Time: %.3f ms\n", t2 * 1000 )

Hhat = sparse(Hhat);
%L = sparse(L);
%L_inv = sparse(L_inv);

%% Test Power Iteration
t1 = tic;
sigma = power_iteration(Hhat, z_ref, eps_buff = 0.1);
t2 = toc(t1);
fprintf("Power Iteration Time: %.3f ms\n", t2 * 1000 )

t1 = tic;
sigma_ck = norm(full(Hhat)) .^ 2;
t2 = toc(t1);
fprintf("Exact max(eig(H))^2 Time: %.3f ms\n", t2 * 1000 )

%assert(abs(abs(sigma - sigma_ck) / sigma_ck - 0.1) < 1e-4, "Power iteration lossy")

%% Solve with CVX First
[z_star_cvx, sol_info_cvx] = cvx_PIPG(Q, q, c, H, h, convex_constraints, {});

1/2 * z_star_cvx' * Q * z_star_cvx + q' * z_star_cvx + c

%% Solve with CVX First
[z_star_cvx, sol_info_cvx] = cvx_PIPG_scaled(Q, q, c, H, h, convex_constraints, {}, S_z, c_z);

%%
[z_star_cvx_sc, sol_info_cvx_sc] = cvx_PIPG_scaled_constraints_only(Qtilde, qtilde, ctilde, Htilde, htilde, convex_constraints, {}, S_z, c_z);

z_tilde = S_z \ (z_star_cvx_sc - c_z);
z_hat = L * z_tilde;

%fprintf("CVX_pc Optimal value: %.5f, solution %.3f%s off from MOSEK\n", 1/2 * z_star_cvx_sc' * Q * z_star_cvx_sc + q' * z_star_cvx_sc, norm(z_star_cvx - z_star_cvx_sc) / norm(z_star_cvx) * 100, "%")
%%
1/2 * z_tilde' * Qtilde * z_tilde + qtilde' * z_tilde + ctilde 

(lambda/2 * z_hat' * z_hat + qhat' * z_hat + ctilde * lambda) / lambda


%%
[z_star_cvx_pc, sol_info_cvx_pc] = cvx_PIPG_pc(lambda, qhat, Hhat, htilde, L_inv, convex_constraints, {}, S_z, c_z);

fprintf("CVX_pc Optimal value: %.5f, solution %.3f%s off from MOSEK\n", 1/2 * z_star_cvx_pc' * Q * z_star_cvx_pc + q' * z_star_cvx_pc, norm(z_star_cvx - z_star_cvx_pc) / norm(z_star_cvx) * 100, "%")

1/2 * z_star_cvx_pc' * Qtilde * z_star_cvx_pc + qtilde' * z_star_cvx_pc + ctilde

%z_ref = z_star_cvx;

%% Test PIPG
rho = 1;
j_check = 1;
j_max = 1000;
[z_star_np, w_star_np, sol_info_np] = PIPG_orig_precond(qhat, Hhat, htilde, D, L, L_inv, lambda, sigma, rho, tol_abs, tol_rel, tol_infeas, j_check, 0, j_max, D.project_onto_cones(1, z_ref), what_ref, S_z, c_z);
% %rho = 1.9;
% %[z_star_np, w_star_np, sol_info_np] = PIPG_no_preconditioning(Q, q, H, h, D, lambda, sigma, rho, tol_abs, tol_rel, j_check, 5000000, j_max, D.project_onto_cones(1, z_ref), what_ref);
%% 
figure
plot(1:j_check:sol_info_np.iterations, abs(sol_info_np.what_inf_dj)); hold on
plot(1:j_check:sol_info_np.iterations, abs(sol_info_np.zhat_inf_dj)); hold off
title("PIPG Iterations")
subtitle(sprintf("Status: %s, Primal Infeasible: %s, Dual Infeasible: %s", sol_info_np.solution_status, string(sol_info_np.primal_infeasible), string(sol_info_np.dual_infeasible)))
xlabel("Iteration")
ylabel("Convergence Metric")
yscale("log")
legend("w_{\infty\Delta j}", "z_{\infty\Delta j}")
grid on

    %obs(l) = 1/2 * z_star' * Q * z_star + q' * z_star;
    %mos_dis(l) = norm(z_star_cvx - z_star);

fprintf("PIPG Optimal value: %.5f, solution %.3f%s off from MOSEK\n", 1/2 * z_star_np' * Q * z_star_np + q' * z_star_np, norm(z_star_cvx - z_star_np) / norm(z_star_cvx) * 100, "%")
%end
% 
% figure
% plot(vecnorm(sol_info_np.zhis - z_star_cvx))
% yscale("log")
% grid on

%%


%lambda_gains = 1.079:0.001:1.081;

%obs = zeros(size(lambda_gains));
%mos_dis = zeros(size(lambda_gains));

%for l = 1:numel(lambda_gains)
rho = 1.9;
%lambda = lambda_gains(l) * norm(Q); % lambda should be >= norm(Q)
j_check = 1;
[z_star, what_star, sol_info] = PIPG_scaled(qhat, Hhat, htilde, D, L, L_inv, lambda, sigma, rho, tol_abs, tol_rel, tol_infeas, j_check, 1000, L_inv * S_z \ (z_star_cvx - c_z), what_ref, S_z, c_z);
% 
figure
plot((1:numel(sol_info.what_inf_dj)) * j_check, abs(sol_info.what_inf_dj)); hold on
plot((1:numel(sol_info.zhat_inf_dj)) * j_check, abs(sol_info.zhat_inf_dj)); hold off
title("PIPG Iterations")
subtitle(sprintf("Status: %s, Primal Infeasible: %s, Dual Infeasible: %s", sol_info.solution_status, string(sol_info.primal_infeasible), string(sol_info.dual_infeasible)))
xlabel("Iteration")
ylabel("Convergence Metric")
yscale("log")
legend("w_{\infty\Delta j}", "z_{\infty\Delta j}")
grid on

    %obs(l) = 1/2 * z_star' * Q * z_star + q' * z_star;
    %mos_dis(l) = norm(z_star_cvx - z_star);

fprintf("PIPG Optimal value: %.5f compared to %.5f from MOSEK, solution %.3f%s off from MOSEK\n", 1/2 * z_star' * Q * z_star + q' * z_star, 1/2 * z_star_cvx' * Q * z_star_cvx + q' * z_star_cvx, norm(z_star_cvx - z_star) / norm(z_star_cvx) * 100, "%")
%end
%%



%lambda_gains = 1.079:0.001:1.081;

%obs = zeros(size(lambda_gains));
%mos_dis = zeros(size(lambda_gains));

%for l = 1:numel(lambda_gains)
rho = 1.9;
%lambda = lambda_gains(l) * norm(Q); % lambda should be >= norm(Q)
j_check = 1;
omega = 1e-5;
[z_star, what_star, sol_info] = PIPG_scaled_new(qhat, Hhat, htilde, D, L, L_inv, lambda, sigma, omega, rho, tol_abs, tol_rel, tol_infeas, j_check, 1000, L_inv * S_z \ (z_star_cvx - c_z), what_ref, S_z, c_z);
% 
figure
plot((1:numel(sol_info.what_inf_dj)) * j_check, abs(sol_info.what_inf_dj)); hold on
plot((1:numel(sol_info.zhat_inf_dj)) * j_check, abs(sol_info.zhat_inf_dj)); hold off
title("PIPG Iterations")
subtitle(sprintf("Status: %s, Primal Infeasible: %s, Dual Infeasible: %s", sol_info.solution_status, string(sol_info.primal_infeasible), string(sol_info.dual_infeasible)))
xlabel("Iteration")
ylabel("Convergence Metric")
yscale("log")
legend("w_{\infty\Delta j}", "z_{\infty\Delta j}")
grid on

    %obs(l) = 1/2 * z_star' * Q * z_star + q' * z_star;
    %mos_dis(l) = norm(z_star_cvx - z_star);

fprintf("PIPG Optimal value: %.5f compared to %.5f from MOSEK, solution %.3f%s off from MOSEK\n", 1/2 * z_star' * Q * z_star + q' * z_star, 1/2 * z_star_cvx' * Q * z_star_cvx + q' * z_star_cvx, norm(z_star_cvx - z_star) / norm(z_star_cvx) * 100, "%")
%end


%%

D.eval_constraint_functions(z_star)
D.eval_constraint_functions(z_star_cvx)

%%
sparse(H * z_star - h)
% sparse(H * z_star_cvx - h)
% sparse(H * z_star_cvx_pc - h)