function [ptr_sol] = ptr(prob, ptr_ops)
%PTR Sequential Convex Programming algorithm
%   If converged, solution satisfies the nonlinear continuous-time equations of motion
% to within a tolerance on the order of eps_feasible feasible, satisfies all algebraic constraints at each
% temporal node, and approximates (local) optimality of the original optimal control problem.

x_ref = zeros([prob.n.x, prob.N, ptr_ops.iter_max + 1]);
u_ref = zeros([prob.n.u, prob.Nu, ptr_ops.iter_max + 1]);
p_ref = zeros([prob.n.p, ptr_ops.iter_max + 1]);

x_ref(:, :, 1) = prob.scale_x(prob.guess.x);
u_ref(:, :, 1) = prob.scale_u(prob.guess.u);
p_ref(:, 1) = prob.scale_p(prob.guess.p);

ptr_sol.converged = false;
ptr_sol.objective = zeros([1, ptr_ops.iter_max + 1]);
ptr_sol.Delta = zeros([prob.Nu, ptr_ops.iter_max + 1]);
ptr_sol.delta_xp = zeros([ptr_ops.iter_max]);

% Convexify along initial guess
[prob, ptr_sol.Delta(:, 1)] = convexify_along_reference(prob, prob.guess.x, prob.guess.u, prob.guess.p);

ptr_sol.objective(1) = prob.objective(prob.guess.x, prob.guess.u, prob.guess.p);

for i = 1:(ptr_ops.iter_max)
    % Solve convex subproblem and update reference
    if prob.n.p == 0
        [x_ref(:, :, i + 1), u_ref(:, :, i + 1), ptr_sol.objective(i + 1)] = solve_ptr_convex_subproblem_no_p(prob, ptr_ops, x_ref(:, :, i), u_ref(:, :, i));
    else
        [x_ref(:, :, i + 1), u_ref(:, :, i + 1), p_ref(:, i + 1), ptr_sol.objective(i)] = solve_ptr_convex_subproblem(prob, ptr_ops, x_ref(:, :, i), u_ref(:, :, i), p_ref(:, i));
    end

    % Convexify along reference trajectory
    [prob, ptr_sol.Delta(:, i + 1)] = convexify_along_reference(prob, prob.unscale_x(x_ref(:, :, i + 1)), prob.unscale_u(u_ref(:, :, i + 1)), prob.unscale_p(p_ref(:, i + 1)));

    % Update algorithm weights (4.24)
    ptr_ops.w_tr = update_trust_region_weights(ptr_sol.Delta(:, i + 1)', ptr_ops.update_w_tr, ptr_ops.w_tr, ptr_ops.Delta_min);

    % Check stopping criteria (4.30)
    ptr_sol.delta_xp(i) = ptr_stopping(x_ref(:, i + 1), p_ref(:, i + 1), x_ref(:, i), p_ref(:, i), ptr_ops.q);
    
    if ptr_sol.delta_xp(i) < ptr_ops.delta_tol && ~ptr_sol.converged
        ptr_sol.converged = true;
        ptr_sol.converged_i = i;
    end
end

if ptr_sol.converged == false
    warning("PTR did not converge after %g iterations. delta_xp = %.3f. norm(Delta) = %.3f\n", i, ptr_sol.delta_xp(i), norm(ptr_sol.Delta(:, end)))
end

ptr_sol.x = prob.unscale_x(x_ref);
ptr_sol.u = prob.unscale_u(u_ref);
ptr_sol.p = prob.unscale_p(p_ref);
end

function [w_tr] = update_trust_region_weights(Delta, update, w_tr, Delta_min)
    w_tr_update = 1 ./ abs(Delta) .* (abs(Delta) >= Delta_min) ...
             + 1 / Delta_min * (abs(Delta_min) < Delta_min);

    w_tr = update .* w_tr_update + ~update .* w_tr;
end

function [delta_xp] = ptr_stopping(X, p, x_ref, p_ref, q)
    delta_xp = zero_if_empty(vecnorm(p - p_ref, q, 2)) + max(vecnorm(X - x_ref, q, 2), [], 1);
end