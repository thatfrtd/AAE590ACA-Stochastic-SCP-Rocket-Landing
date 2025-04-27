function [ptr_sol] = Stochastic_ptr_memoryefficient(prob, ptr_ops)
%PTR Sequential Convex Programming algorithm
%   If converged, solution satisfies the nonlinear continuous-time equations of motion
% to within a tolerance on the order of eps_feasible feasible, satisfies all algebraic constraints at each
% temporal node, and approximates (local) optimality of the original optimal control problem.

x_sol = zeros([prob.n.x, prob.N]);
u_sol = zeros([prob.n.u, prob.Nu]);
p_sol = zeros([prob.n.p]);
X_sol = zeros([prob.n.x, prob.n.x * prob.N * (prob.N + 1) / 2]);
S_sol = zeros([prob.n.u, prob.n.x * prob.N * (prob.N + 1) / 2]);

ptr_sol.converged = false;
ptr_sol.objective = zeros([1, ptr_ops.iter_max + 1]);
ptr_sol.Delta = zeros([prob.Nu, ptr_ops.iter_max + 1]);
ptr_sol.delta_xp = zeros([1, ptr_ops.iter_max + 1]);

% Convexify along initial guess
[prob, ptr_sol.Delta(:, 1)] = convexify_along_reference(prob, prob.guess.x, prob.guess.u, prob.guess.p);

disp(" k |       status      |   vd  |   vs  |  vbc_NP |  vbc_N |    J    |   J_tr  |   J_vc   |   dJ %  |   dx   |   du   |   dp   | delta |  dyn  |  eta  | eta_x | eta_u | eta_p")

for i = 1:(ptr_ops.iter_max)
    % Solve convex subproblem and update reference
    if prob.n.p == 0
        [x_sol, u_sol, X_sol, S_sol, sol_info] = solve_stochastic_ptr_convex_subproblem_no_p_2(prob, ptr_ops, prob.scale_x(prob.guess.x), prob.scale_u(prob.guess.u), prob.guess.X, prob.guess.S);
    else
        [x_sol, u_sol, p_sol, sol_info] = solve_stochastic_ptr_convex_subproblem(prob, ptr_ops, prob.scale_x(prob.guess.x), prob.scale_u(prob.guess.u), prob.scale_p(prob.guess.p));
    end

    % Convexify along reference trajectory
    [prob, ptr_sol.Delta(:, i + 1)] = convexify_along_reference(prob, prob.unscale_x(x_sol), prob.unscale_u(u_sol), prob.unscale_p(p_sol));

    % Update algorithm weights (4.24)
    ptr_ops.w_tr = update_trust_region_weights(ptr_sol.Delta(:, i + 1)', ptr_ops.update_w_tr, ptr_ops.w_tr, ptr_ops.Delta_min);

    % Check stopping criteria (4.30)
    ptr_sol.delta_xp(i) = ptr_stopping(x_sol, p_sol, prob.scale_x(prob.guess.x), prob.scale_p(prob.guess.p), ptr_ops.q);

    % Display results of iteration
    if i ~= 1
        ptr_sol.info(i).dJ = ptr_sol.info(i - 1);
    end
    sol_info.delta = ptr_sol.delta_xp(i);
    sol_info.dyn = norm(ptr_sol.Delta(:, i + 1)) <= ptr_ops.Delta_min;

    ptr_sol.info(i) = sol_info;

    if i == 1
        fprintf("%2.f | %17s | %5.1g | %5.1g | %5.1g | %5.1g | %5.3g | %5.3g | %5.3g |         | %6.1g | %6.1g | %5.1g | %5.1g | %5s | %5.1g | %5.1g | %5.1g | %5.1g\n", i, ptr_sol.info(i).status, norm(ptr_sol.info(i).vd), norm(ptr_sol.info(i).vs), norm(ptr_sol.info(i).vbc_NP), norm(ptr_sol.info(i).vbc_N), ptr_sol.info(i).J, ptr_sol.info(i).J_tr, ptr_sol.info(i).J_vc, sum(ptr_sol.info(i).dx), sum(ptr_sol.info(i).du), sum(ptr_sol.info(i).dp), ptr_sol.info(i).delta, string(ptr_sol.info(i).dyn), norm(ptr_sol.info(i).eta), norm(ptr_sol.info(i).eta_x), norm(ptr_sol.info(i).eta_u), norm(ptr_sol.info(i).eta_p))
    else
        fprintf("%2.f | %17s | %5.1g | %5.1g | %5.1g | %5.1g | %5.3g | %5.3g | %5.3g | %6.3f | %6.1g | %6.1g | %5.1g | %5.1g | %5s | %5.1g | %5.1g | %5.1g | %5.1g\n", i, ptr_sol.info(i).status, norm(ptr_sol.info(i).vd), norm(ptr_sol.info(i).vs), norm(ptr_sol.info(i).vbc_NP), norm(ptr_sol.info(i).vbc_N), ptr_sol.info(i).J, ptr_sol.info(i).J_tr, ptr_sol.info(i).J_vc, ptr_sol.info(i).dJ, sum(ptr_sol.info(i).dx), sum(ptr_sol.info(i).du), sum(ptr_sol.info(i).dp), ptr_sol.info(i).delta, string(ptr_sol.info(i).dyn), norm(ptr_sol.info(i).eta), norm(ptr_sol.info(i).eta_x), norm(ptr_sol.info(i).eta_u), norm(ptr_sol.info(i).eta_p))
    end

    if i >= ptr_ops.iter_min
        if ptr_sol.delta_xp(i) < ptr_ops.delta_tol && ~ptr_sol.converged && sol_info.dyn
            ptr_sol.converged = true;
            ptr_sol.converged_i = i + 1;
    
            break
        end
    end

    prob.guess.x = prob.unscale_x(x_sol);
    prob.guess.u = prob.unscale_u(u_sol);
    prob.guess.p = prob.unscale_p(p_sol);
    prob.guess.X = X_sol;
    prob.guess.S = S_sol;
end

if ptr_sol.converged == false
    warning("PTR did not converge after %g iterations. delta_xp = %.3f. norm(Delta) = %.3f\n", i, ptr_sol.delta_xp(i), norm(ptr_sol.Delta(:, end)))
end

ptr_sol.x = prob.unscale_x(x_sol);
ptr_sol.u = prob.unscale_u(u_sol);
ptr_sol.p = prob.unscale_p(p_sol);
ptr_sol.X = X_sol;
ptr_sol.S = S_sol;
end

function [w_tr] = update_trust_region_weights(Delta, update, w_tr, Delta_min)
    w_tr_update = 1 ./ abs(Delta) .* (abs(Delta) >= Delta_min) ...
             + 1 / Delta_min * (abs(Delta) < Delta_min);

    w_tr = update .* w_tr_update + ~update .* w_tr;
end

function [delta_xp] = ptr_stopping(X, p, x_ref, p_ref, q)
    delta_xp = zero_if_empty(vecnorm(p - p_ref, q, 1)) + max(vecnorm(X - x_ref, q, 1));
end