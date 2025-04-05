function [ptr_sol] = ptr(prob, ptr_ops)
%PTR Sequential Convex Programming algorithm
%   If converged, solution satisfies the nonlinear continuous-time equations of motion
% to within a tolerance on the order of eps_feasible feasible, satisfies all algebraic constraints at each
% temporal node, and approximates (local) optimality of the original optimal control problem.

x_ref = zeros([prob.N, ptr_ops.iter_max]);
u_ref = zeros([prob.Nu, ptr_ops.iter_max]);
p_ref = zeros([size(p_guess), ptr_ops.iter_max]);

x_ref(:, 1) = prob.guess.x;
u_ref(:, 1) = prob.guess.u;
p_ref(:, 1) = prob.guess.p;

% Convexify along initial guess
prob = convexify_along_reference(prob, 1);

ptr_sol.converged = false;

for i = 1:ptr_ops.iter_max
    % Solve convex subproblem and update reference
    prob = solve_ptr_convex_subproblem(prob, i);
    
    % Convexify along reference trajectory
    prob = convexify_along_reference(prob, i);

    % Determine feasibility of reference trajectory
    ptr_sol.Delta = calculate_defect(prob, i + 1);

    % Update algorithm weights (4.24)
    ptr_ops.w_tr = update_trust_region_weights(ptr_sol.Delta, update, ptr_ops.w_tr, ptr_ops.Delta_min);

    % Check stopping criteria (4.30)
    ptr_sol.delta_xp = ptr_stopping(x_ref(:, i + 1), p_ref(:, i + 1), x_ref(:, i), p_ref(:, i), ptr_ops.q, prob.S_x, prob.S_p);
    
    if ptr_sol.delta_xp < ptr_ops.delta_tol
        ptr_sol.converged = true;
        ptr_sol.converged_i = i;
    end
end

if ptr_sol.converged == false
    warning("PTR did not converge after %g iterations. delta_xp = %.3f. norm(Delta) = %.3f\n", i, ptr_sol.delta_xp, ptr_sol.Delta)
end

ptr_sol.x = x_ref;
ptr_sol.u = u_ref;
ptr_sol.p = p_ref;
end

function [w_tr] = update_trust_region_weights(Delta, update, w_tr, Delta_min)
    w_tr_update = 1 ./ abs(Delta) * (abs(Delta) >= Delta_min) ...
             + 1 / Delta_min * (abs(Delta_min) < Delta_min);

    w_tr = update .* w_tr_update + ~update .* w_tr;
end

function [delta_xp] = ptr_stopping(X, p, x_ref, p_ref, q, S_x, S_p)
    delta_xp = vecnorm(S_p \ (p - p_ref), q, 2) + max(vecnorm(S_x \ (X - x_ref), q, 2), [], 1);
end