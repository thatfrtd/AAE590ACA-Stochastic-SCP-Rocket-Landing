function [x_sol, u_sol, sol_info] = solve_stochastic_ptr_convex_subproblem_no_p(prob, ptr_ops, x_ref, u_ref)
%SOLVE_PTR_CONVEX_SUBPROBLEM Summary of this function goes here
%   Detailed explanation goes here
P_yk_sqrt = get_P_k_root(prob.disc.P, prob.disc.C, prob.disc.D);

cvx_begin quiet
    variable X(prob.n.x, prob.N)
    variable U(prob.n.u, prob.Nu)
    variable eta(1, prob.Nu)
    variable V(prob.n.x, prob.N - 1)
    variable v_prime(prob.n.ncvx)
    variable v_0(prob.n.x, 1)
    variable v_N(prob.n.x, 1)
    variable X_C(prob.n.x, ) % This is for the covariance propoagation. This mf changes dimensions every iteration
    variable S_k(prob.n.u, (prob.N^2 + prob.N)/2 *prob.n.x)
    minimize( prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0) ...
        + virtual_control_cost(V, v_prime, v_0, v_N, ptr_ops.w_vc) ...
        + trust_region_cost(eta, 0, ptr_ops.w_tr, 0) ) % Need to update this to include minimization of S
    subject to
        
        % Stochastic/Kalman Implementaion
        if prob.u_hold == "ZOH"
            for k = 1:(prob.N - 1)
                X(:, k + 1) == prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                             + prob.disc.B_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                             + prob.disc.E_k(:, :, k) * prob.unscale_p(p) ...
                             + prob.disc.c_k(:, :, k) ...
                             + V(:, k);
 
            
                X_C(:, 1:((k+1)*(k+2)/2)*prob.n.x) == [A*X + B*S_k, ...
                                                   prob.disc.L(:,:,k) * P_yk_sqrt(:,:,k)]; % Dimensions don't work here, will have to replace with new method of doing it

                P_hat_N = X_C * X_C'; % Need to fix dimensions of X

                norm(sqrtm(inv(prob.Pf - P_hat_N))*X_C, 2) - 1 <= 0; % Covariance path constraint
            end
        elseif prob.u_hold == "FOH"
            for k = 1:(prob.N - 1)
                X(:, k + 1) == prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                             + prob.disc.B_minus_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                             + prob.disc.B_plus_k(:, :, k) * prob.unscale_u(U(:, k + 1)) ...
                             + prob.disc.E_k(:, :, k) * prob.unscale_p(p) ...
                             + prob.disc.c_k(:, :, k) ...
                             + V(:, k);
            end
        end
        % Constraints
        for k = 1:prob.Nu
            % Convex Constraints
            for cc = 1:prob.n.cvx
                prob.convex_constraints{cc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0) <= 0;
            end
            % Nonconvex Constraints
            for nc = 1:prob.n.ncvx
                prob.nonconvex_constraints{nc}(prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), 0, x_ref, u_ref, 0) ...
                    - v_prime(nc) <= 0;
            end
        end
        v_prime >= 0;

        % Boundary Conditions
        prob.initial_bc(prob.unscale_x(X(:, 1)), 0) + v_0 == 0;
        prob.terminal_bc(prob.unscale_x(X(:, prob.N)), 0) + v_N == 0;

        % Trust Region Constraints
        ptr_ops.alpha_x * norms(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), ptr_ops.q, 1) + ptr_ops.alpha_u * norms(U - u_ref, ptr_ops.q, 1) <= eta;
cvx_end

x_sol = X;
u_sol = U;

sol_info.status = cvx_status;
sol_info.vd = V;
sol_info.vs = v_prime;
sol_info.vbc_0 = v_0;
sol_info.vbc_N = v_N;
sol_info.J = prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0);
sol_info.J_tr = trust_region_cost(eta, 0, ptr_ops.w_tr, 0);
sol_info.J_vc = virtual_control_cost(V, v_prime, v_0, v_N, ptr_ops.w_vc);
sol_info.dJ = 100 * (prob.objective(prob.unscale_x(X), prob.unscale_u(U), 0) - prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0)) / prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), 0);
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

function [J_vc] = virtual_control_cost(V, v_prime, v_0, v_N, w_vc)
    J_vc = w_vc * (norm(v_prime, 1) + sum(norms(V, 1, 1)) + norm(v_0, 1) + norm(v_N, 1));
end