function [x_sol, u_sol, p_sol, sol_info] = solve_ptr_convex_subproblem_phase(prob, ptr_ops, x_ref, u_ref, p_ref)
%SOLVE_PTR_CONVEX_SUBPROBLEM Summary of this function goes here
%   Detailed explanation goes here

cvx_begin quiet
    variable X(prob.n.x, prob.N)
    variable U(prob.n.u, prob.Nu)
    variable p(prob.n.p, 1)
    variable eta(1, prob.Nu)
    variable eta_p(1, prob.n.p)
    variable V(prob.n.x, prob.N - 1)
    variable v_prime(prob.n.ncvx)
    variable v_trans(prob.n.pncvx)
    variable v_0(prob.n.x, 1)
    variable v_N(prob.n.x, 1)
    minimize( prob.objective(prob.unscale_x(X), prob.unscale_u(U), prob.unscale_p(p), prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref)) ...
        + virtual_control_cost(V, [v_prime v_trans], v_0, v_N, ptr_ops.w_vc) ...
        + trust_region_cost(eta, eta_p, ptr_ops.w_tr, ptr_ops.w_tr_p) )
    subject to
        % Dynamics
        if prob.u_hold == "ZOH"
            for k = 1:(prob.N - 1)
                if ~ismember(k, prob.phase_transition_k)
                    X(:, k + 1) == prob.scale_x(prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                                 + prob.disc.B_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                                 + prob.disc.E_k(:, :, k) * prob.unscale_p(p) ...
                                 + prob.disc.c_k(:, :, k) ...
                                 + V(:, k));
                end
            end
        elseif prob.u_hold == "FOH"
            for k = 1:(prob.N - 1)
                if ~ismember(k, prob.phase_transition_k)
                    X(:, k + 1) == prob.scale_x(prob.disc.A_k(:, :, k) * prob.unscale_x(X(:, k)) ...
                                 + prob.disc.B_minus_k(:, :, k) * prob.unscale_u(U(:, k)) ...
                                 + prob.disc.B_plus_k(:, :, k) * prob.unscale_u(U(:, k + 1)) ...
                                 + prob.disc.E_k(:, :, k) * prob.unscale_p(p) ...
                                 + prob.disc.c_k(:, :, k) ...
                                 + V(:, k));
                end
            end
        end

        % Constraints
        for k = 1:prob.Nu
            % Convex Constraints
            for cc = 1:prob.n.cvx
                cc_k = prob.convex_constraints{cc}{1};
                if ismember(k, cc_k)
                    cvx_constraint_func = prob.convex_constraints{cc}{2};
                    cvx_constraint_func(prob.t_k(k), prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), prob.unscale_p(p)) <= 0;
                end
            end
            % Nonconvex Constraints
            for nc = 1:prob.n.ncvx
                nc_k = prob.nonconvex_constraints{nc}{1};
                if ismember(k, nc_k)
                    ncvx_constraint_func = prob.nonconvex_constraints{nc}{2};
                    ncvx_constraint_func(prob.t_k(k), prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), prob.unscale_p(p), prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref), k) ...
                        - v_prime(nc) <= 0;
                end
            end
        end

        % Phase transition constraints
        pncvx_count = 0;
        for e = 1:(prob.n.phase - 1)
            ph_k = prob.phase_transition_k(e);
            % Convex Constraints
            for pcc = 1:prob.phase_transition{e}.n.cvx
                cvx_constraint_func = prob.phase_transition{e}.convex_constraints{pcc}{1};
                cvx_constraint_type = prob.phase_transition{e}.convex_constraints{pcc}{2};
                if cvx_constraint_type == "<="
                    cvx_constraint_func(prob.t_k(ph_k), prob.unscale_x(X(:, [ph_k, ph_k + 1])), prob.unscale_u(U(:, [ph_k, ph_k + 1])), prob.unscale_p(p)) <= 0;
                elseif cvx_constraint_type == "=="
                    cvx_constraint_func(prob.t_k(ph_k), prob.unscale_x(X(:, [ph_k, ph_k + 1])), prob.unscale_u(U(:, [ph_k, ph_k + 1])), prob.unscale_p(p)) == 0;
                end
            end
            % Nonconvex Constraints
            for pnc = 1:prob.phase_transition{e}.n.ncvx
                pncvx_count = pncvx_count + 1;
                ncvx_constraint_func = prob.phase_transition{e}.nonconvex_constraints{pnc}{1};
                ncvx_constraint_type = prob.phase_transition{e}.nonconvex_constraints{pnc}{2};
                if ncvx_constraint_type == "<="
                    ncvx_constraint_func(prob.t_k(ph_k), prob.unscale_x(X(:, [ph_k, ph_k + 1])), prob.unscale_u(U(:, [ph_k, ph_k + 1])), prob.unscale_p(p), prob.unscale_x(x_ref(:, [ph_k, ph_k + 1])), prob.unscale_u(u_ref(:, [ph_k, ph_k + 1])), prob.unscale_p(p_ref)) ...
                        - v_trans(pncvx_count) <= 0;
                    v_trans(pncvx_count) >= 0;
                elseif ncvx_constraint_type == "=="
                    ncvx_constraint_func(prob.t_k(ph_k), prob.unscale_x(X(:, [ph_k, ph_k + 1])), prob.unscale_u(U(:, [ph_k, ph_k + 1])), prob.unscale_p(p), prob.unscale_x(x_ref(:, [ph_k, ph_k + 1])), prob.unscale_u(u_ref(:, [ph_k, ph_k + 1])), prob.unscale_p(p_ref)) ...
                        - v_trans(pncvx_count) == 0;
                end
            end
        end

        %p(1) == 2.62;
        %p(2) == 10;

        % Boundary Conditions
        prob.initial_bc(prob.unscale_x(X(:, 1)), prob.unscale_p(p), prob.unscale_x(x_ref(:, 1)), prob.unscale_p(p_ref)) + v_0 == 0;
        prob.terminal_bc(prob.unscale_x(X(:, prob.N)), prob.unscale_p(p), prob.unscale_x(x_ref(:, prob.N)), prob.unscale_p(p_ref)) + v_N == 0;

        % Trust Region Constraints
        ptr_ops.alpha_x * sum(sum_square(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu))) + ptr_ops.alpha_u * sum(sum_square(U - u_ref)) <= eta;
        ptr_ops.alpha_p * norms(p - p_ref, ptr_ops.q, 1) <= eta_p;
cvx_end
% 
% ck = zeros([prob.n.cvx, prob.Nu]);
% 
% for k = 1:prob.Nu
%     % Convex Constraints
%     for cc = 1:prob.n.cvx
%         ck(cc, k) = prob.convex_constraints{cc}(prob.t_k(k), prob.unscale_x(X(:, k)), prob.unscale_u(U(:, k)), prob.unscale_p(p));
%     end
% end
% 
% figure
% plot(ck(1, :)); hold on
% plot(ck(2, :)); hold on
% plot(ck(3, :)); hold off

x_sol = X;
u_sol = U;
p_sol = p;

sol_info.status = cvx_status;
sol_info.vd = V;
sol_info.vs = [v_prime v_trans];
sol_info.vbc_0 = v_0;
sol_info.vbc_N = v_N;
sol_info.J = prob.objective(prob.unscale_x(X), prob.unscale_u(U), prob.unscale_p(p), prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref));
sol_info.J_tr = trust_region_cost(eta, eta_p, ptr_ops.w_tr, ptr_ops.w_tr_p);
sol_info.J_vc = virtual_control_cost(V, [v_prime v_trans], v_0, v_N, ptr_ops.w_vc);
sol_info.dJ = 100 * (prob.objective(prob.unscale_x(X), prob.unscale_u(U), prob.unscale_p(p), prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref)) ...
    - prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref), prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref))) ...
    / prob.objective(prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref), prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref));
sol_info.dx = vecnorm(X(:, 1:prob.Nu) - x_ref(:, 1:prob.Nu), ptr_ops.q, 1);
sol_info.du = vecnorm(U - u_ref, ptr_ops.q, 1);
sol_info.dp = vecnorm(p - p_ref, ptr_ops.q, 1);
sol_info.eta = eta;
sol_info.eta_x = 0;
sol_info.eta_u = 0;
sol_info.eta_p = eta_p;
end

function [J_tr] = trust_region_cost(eta, eta_p, w_tr, w_tr_p)
    J_tr = w_tr * eta' + w_tr_p * eta_p';
end

function [J_vc] = virtual_control_cost(V, v_prime, v_0, v_N, w_vc)
    J_vc = w_vc * (norm(v_prime, 1) + sum(norms(V, 1, 1)) + norm(v_0, 1) + norm(v_N, 1));
end