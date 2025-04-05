function [prob] = solve_ptr_convex_subproblem(prob, i)
%SOLVE_PTR_CONVEX_SUBPROBLEM Summary of this function goes here
%   Detailed explanation goes here

cvx_begin
    variable X(prob.n_x, prob.N)
    variable U(prob.n_u, prob.N_u)
    variable p(prob.n_p, 1)
    variable eta(prob.N)
    variable eta_p(prob.np, 1)
    variable V(prob.nx, prob.N - 1)
    variable v_prime(prob.n_ncvx)
    variable v_0(prob.n_x, 1)
    variable v_N(prob.n_x, 1)
    minimize( prob.objective(X, U, p) ...
        + virtual_control_cost(V, v_prime, v_0, v_N, w_vc) ...
        + trust_region_cost(eta, eta_p, ptr_ops.w_tr, ptr_ops.w_tr_p) )
    subject to
        % Dynamics
        if prob.u_hold == "ZOH"
            for k = 1:(prob.N - 1)
                X(:, k + 1) == prob.disc.A_k(:, :, k) * X(:, k) ...
                             + prob.disc.B_k(:, :, k) * U(:, k) ...
                             + prob.disc.E_k(:, :, k) * p ...
                             + prob.disc.c_k(:, :, k) ...
                             + V(:, k);
            end
        elseif prob.u_hold == "FOH"
            for k = 1:(prob.N - 1)
                X(:, k + 1) == prob.disc.A_k(:, :, k) * X(:, k) ...
                             + prob.disc.B_minus_k(:, :, k) * U(:, k) ...
                             + prob.disc.B_plus_k(:, :, k) * U(:, k + 1) ...
                             + prob.disc.E_k(:, :, k) * p ...
                             + prob.disc.c_k(:, :, k) ...
                             + V(:, k);
            end
        end

        % Constraints
        for k = 1:prob.N
            % Convex Constraints
            for cc = 1:prob.n_cvx
                prob.convex_constraints{cc}(X(:, k), U(:, k), p) <= 0;
            end
            % Nonconvex Constraints
            for nc = 1:prob.n_ncvx
                prob.nonconvex_constraints{nc}(X(:, k), U(:, k), p, prob.x_ref(:, :, i), prob.u_ref(:, :, i), prob.p_ref(:, i)) ...
                    - v_prime(nc) <= 0;
            end
        end
        v_prime >= 0;

        % Boundary Conditions
        prob.initial_bc(X(:, 1), p) + v_0 == 0;
        prob.terminal_bc(X(:, prob.N), p) + v_N == 0;

        % Trust Region Constraints
        trust_region_constraints(X, U, p, prob.x_ref(:, :, i + 1), prob.u_ref(:, :, i + 1), prob.p_ref(:, i + 1), ptr_ops.q, ptr_ops.alpha_x, ptr_ops.alpha_u, ptr_ops.alpha_p, eta, eta_p);
cvx_end

end

function [J_tr] = trust_region_cost(eta, eta_p, w_tr, w_tr_p)
    J_tr = w_tr' * eta + w_tr_p * eta_p;
end

function [J_vc] = virtual_control_cost(V, v_prime, v_0, v_N, w_vc)
    J_vc = w_vc * (norm(v_prime, 1) + sum(norms(V, 1, 2)) + norm(v_0, 1) + norm(v_N, 1));
end

function [] = trust_region_constraints(X, U, p, x_ref, u_ref, p_ref, q, alpha_x, alpha_u, alpha_p, eta, eta_p)
    % Limit artificial unboundedness
    alpha_x * norms(X - x_ref, q, 2) ...
        + alpha_u * vecnorm(U - u_ref, q, 2) <= eta;
    alpha_p * norms(p - p_ref, q, 2) <= eta_p;
end

