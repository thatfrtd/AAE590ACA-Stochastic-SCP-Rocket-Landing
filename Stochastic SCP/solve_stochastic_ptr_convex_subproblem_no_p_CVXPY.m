function [x_sol, u_sol, X_sol, S_sol, sol_info] = solve_stochastic_ptr_convex_subproblem_no_p_CVXPY(prob, ptr_ops, x_ref, u_ref, X_k_ref, S_k_ref)
%SOLVE_PTR_CONVEX_SUBPROBLEM Summary of this function goes here
%   Detailed explanation goes here

t1 = tic;

P_yk_sqrt = sqrtm_array(prob.disc.Ptilde_minus_k);

L_kp1_times_P_yk_sqrt = pagemtimes(prob.disc.L_k(:, :, 2:(prob.N)), P_yk_sqrt(:, :, 1:prob.N - 1));

% Solve
if prob.n.x == 7
    [X_py, U_py, X_C_py, S_k_py, eta_py, V_py, v_prime_py, v_0_py, v_N_py, v_NP_py, solve_status_py] = pyrunfile("Stochastic3DoF_clarabel.py", ["X_sol", "U_sol", "X_C_sol", "S_k_sol", "eta", "V", "v_prime", "v_0", "v_N", "v_NP", "solve_status"], x_ref = x_ref, u_ref = u_ref, X_k_ref = X_k_ref, S_k_ref = S_k_ref, x_0 = prob.x0, x_f = prob.xf, Phat0 = prob.Phat0, Phatf = prob.Pf - prob.disc.Ptilde_k(:, :, end), A_k = prob.disc.A_k, B_k = prob.disc.B_k, c_k = prob.disc.c_k, params = prob.params, N = prob.N - 1, delta_t = prob.tf / (prob.N - 1), w_vc = ptr_ops.w_vc, w_tr = ptr_ops.w_tr, L_kp1_times_P_yk_sqrt = L_kp1_times_P_yk_sqrt);
end


% Convert solution to Matlab friendly types
X = double(X_py);
U = double(U_py);
X_C = double(X_C_py);
S_k = double(S_k_py);
eta = double(eta_py);
V = double(V_py);
v_prime = double(v_prime_py);
v_0 = double(v_0_py);
v_N = double(v_N_py);
v_NP = double(v_NP_py);

solve_status = string(solve_status_py);

t2 = toc(t1);

% Check constraint satisfaction
nc_ck = zeros([prob.n.ncvx, prob.Nu]);
nc_ck_noref = zeros([prob.n.ncvx, prob.Nu]);
for k = 1:prob.Nu
    % Nonconvex Constraints
    for nc = 1:prob.n.ncvx
        nc_ck(nc, k) = prob.nonconvex_constraints{nc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, X_C(:, (tri(k - 1, prob.n.x) + 1):tri(k, prob.n.x)), S_k(:, (tri(k - 1, prob.n.x) + 1):tri(k, prob.n.x)), x_ref, u_ref, 0, X_k_ref, S_k_ref, k);
            %- v_prime(nc);
        nc_ck_noref(nc, k) = prob.nonconvex_constraints{nc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, X_C(:, (tri(k - 1, prob.n.x) + 1):tri(k, prob.n.x)), S_k(:, (tri(k - 1, prob.n.x) + 1):tri(k, prob.n.x)), prob.unscale_x(X), prob.unscale_u(U), 0, X_C, S_k, k);
    end

    min_test(k) = 33.9916 / exp(X(7, k)) - (norm(U(:, k)) - 3.89894920704084 * norm(S_k(:, (tri(k - 1, prob.n.x) + 1):tri(k, prob.n.x))));
end

figure; plot(min_test)
% 
% theta = zeros(1, size(U, 2));
% for k = 1:prob.nu
%     theta = (1, k) / norm(U(1:2, k));
% end
x_sol = X;
u_sol = U;
X_sol = X_C;
S_sol = S_k;

fprintf("CVXPY Time: %.3f ms\n", t2 * 1000 )

sol_info.status = solve_status;
sol_info.vd = V;
sol_info.vs = v_prime;
sol_info.vbc_0 = v_0;
sol_info.vbc_N = v_N;
sol_info.vbc_NP = v_NP;
sol_info.J = prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0, X_C, S_k);
sol_info.J_tr = trust_region_cost(eta, 0, ptr_ops.w_tr, 0);
sol_info.J_vc = virtual_control_cost(V, v_prime, v_0, v_N, v_NP, ptr_ops.w_vc);
sol_info.dJ = 100 * (prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0, X_C, S_k) - prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0, 0, S_k_ref)) / prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0, 0, S_k_ref);
sol_info.dx = vecnorm(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), ptr_ops.q, 1);
sol_info.du = vecnorm(U - u_ref, ptr_ops.q, 1);
sol_info.dp = 0;
sol_info.eta = eta;
sol_info.eta_x = 0;
sol_info.eta_u = 0;
sol_info.eta_p = 0;
end

function [J_tr] = trust_region_cost(eta, eta_p, w_tr, w_tr_p)
    J_tr = w_tr * eta' + w_tr_p * eta_p;
end

function [J_vc] = virtual_control_cost(V, v_prime, v_0, v_N, v_NP, w_vc)
    J_vc = w_vc * (norm(v_prime, 1) + v_NP + sum(norms(V, 1, 1)) + norm(v_0, 1) + norm(v_N, 1));
end