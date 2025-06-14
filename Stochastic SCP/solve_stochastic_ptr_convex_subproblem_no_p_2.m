function [x_sol, u_sol, X_sol, S_sol, sol_info] = solve_stochastic_ptr_convex_subproblem_no_p_2(prob, ptr_ops, x_ref, u_ref, X_k_ref, S_k_ref)
%SOLVE_PTR_CONVEX_SUBPROBLEM Summary of this function goes here
%   Detailed explanation goes here

t1 = tic;

P_yk_sqrt = sqrtm_array(prob.disc.Ptilde_minus_k);

L_kp1_times_P_yk_sqrt = pagemtimes(prob.disc.L_k(:, :, 2:(prob.N)), P_yk_sqrt(:, :, 1:prob.N - 1));

tri = @(k) k * (k + 1) / 2 * prob.n.x;

cvx_begin
    variable X(prob.n.x, prob.N)
    variable U(prob.n.u, prob.Nu)
    variable eta(1, prob.Nu)
    variable V(prob.n.x, prob.N - 1)
    variable v_prime(prob.n.ncvx)
    variable v_0(prob.n.x, 1)
    variable v_N(prob.n.x, 1)
    variable v_NP(1,1)
    variable X_C(prob.n.x, tri(prob.N)) % This is for the covariance propoagation. This mf changes dimensions every iteration
    variable S_k(prob.n.u, tri(prob.N))
    minimize( prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0, X_C, S_k) ...
        + virtual_control_cost(V, v_prime, v_0, v_N, v_NP, ptr_ops.w_vc) ...
        + trust_region_cost(eta, 0, ptr_ops.w_tr, 0) )
    subject to
        % Dynamics
        for k = 1:(prob.N - 1)
            X(:, k + 1) == prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                         + prob.disc.B_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                         + prob.disc.c_k(:, :, k) ...
                         + V(:, k);

            %X_C(:, (tri(k) + 1):tri(k + 1)) == [prob.disc.A_k(:, :, k) * X_C(:, (tri(k - 1) + 1):tri(k)) + prob.disc.B_k(:, :, k) * S_k(:, (tri(k - 1) + 1):tri(k)), prob.disc.L_k(:, :, k + 1) * P_yk_sqrt(:, :, k)];
            X_C(:, (tri(k) + 1):tri(k + 1)) == [prob.disc.A_k(:, :, k) * X_C(:, (tri(k - 1) + 1):tri(k)) + prob.disc.B_k(:, :, k) * S_k(:, (tri(k - 1) + 1):tri(k)), L_kp1_times_P_yk_sqrt(:, :, k)];
        end
        % Constraints
        for k = 1:prob.Nu
            % Convex Constraints
            for cc = 1:prob.n.cvx
                prob.convex_constraints{cc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, X_C(:, (tri(k - 1) + 1):tri(k)), S_k(:, (tri(k - 1) + 1):tri(k))) <= 0;
            end
            % Nonconvex Constraints
            for nc = 1:prob.n.ncvx
                prob.nonconvex_constraints{nc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, X_C(:, (tri(k - 1) + 1):tri(k)), S_k(:, (tri(k - 1) + 1):tri(k)), x_ref, u_ref, 0, X_k_ref, S_k_ref, k) ...
                    - v_prime(nc) <= 0;
            end
        end
        v_prime >= 0;

        % Boundary Conditions
        prob.initial_bc(prob.unscale_x(X(:, 1)), 0) + v_0 == 0; % Initial mean state
        sqrtm(prob.Phat0) == X_C(:, 1:prob.n.x); % Initial estimated state covariance
        prob.terminal_bc(prob.unscale_x(X(:, prob.N)), 0) + v_N == 0; % Final mean state
        %norm(sqrtm(prob.Pf - prob.disc.Ptilde_k(:, :, end)) \ X_C(:, (tri(prob.N - 1) + 1):tri(prob.N))) - 1 - v_NP <= 0; % Final estimated state covariance
        norm(inv(sqrtm(prob.Pf - prob.disc.Ptilde_k(:, :, end))) * X_C(:, (tri(prob.N - 1) + 1):tri(prob.N))) - 1 - v_NP <= 0; % Final estimated state covariance

        v_NP >= 0;

        % Trust Region Constraints
        ptr_ops.alpha_x * norms(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), ptr_ops.q, 1) + ptr_ops.alpha_u * norms(U - u_ref, ptr_ops.q, 1) <= eta;
cvx_end

X_k_err = zeros([1, prob.Nu]);
X_k_ref_err = zeros([1, prob.Nu]);
for k = 1:prob.Nu
    X_k_err(k) = norm(X_C(:, (tri(k) + 1):tri(k + 1)) - [prob.disc.A_k(:, :, k) * X_C(:, (tri(k - 1) + 1):tri(k)) + prob.disc.B_k(:, :, k) * S_k(:, (tri(k - 1) + 1):tri(k)), L_kp1_times_P_yk_sqrt(:, :, k)]);
    X_k_ref_err(k) = norm(X_k_ref(:, (tri(k) + 1):tri(k + 1)) - [prob.disc.A_k(:, :, k) * X_k_ref(:, (tri(k - 1) + 1):tri(k)) + prob.disc.B_k(:, :, k) * S_k_ref(:, (tri(k - 1) + 1):tri(k)), L_kp1_times_P_yk_sqrt(:, :, k)]);
end

t2 = toc(t1);

nc_ck = zeros([prob.n.ncvx, prob.Nu]);
nc_ck_noref = zeros([prob.n.ncvx, prob.Nu]);
for k = 1:prob.Nu
    % Nonconvex Constraints
    for nc = 1:prob.n.ncvx
        nc_ck(nc, k) = prob.nonconvex_constraints{nc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, X_C(:, (tri(k - 1) + 1):tri(k)), S_k(:, (tri(k - 1) + 1):tri(k)), x_ref, u_ref, 0, X_k_ref, S_k_ref, k);
            %- v_prime(nc);
        nc_ck_noref(nc, k) = prob.nonconvex_constraints{nc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, X_C(:, (tri(k - 1) + 1):tri(k)), S_k(:, (tri(k - 1) + 1):tri(k)), prob.unscale_x(X), prob.unscale_u(U), 0, X_C, S_k, k);
    end

    min_test(k) = 33.9916 / exp(X(7, k)) - (norm(U(:, k)) - 3.89894920704084 * norm(S_k(:, (tri(k - 1) + 1):tri(k))));
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

fprintf("MOSEK Time: %.3f ms\n", t2 * 1000 )

sol_info.status = cvx_status;
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