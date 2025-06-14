function [x_star, xi_star, u_star, s_star, sol_info] = most_general_cvx_PIPG(x_ref, u_ref, s_ref, A_k, B_k_minus, B_k_plus, S_k, d_k, c_x_i, c_x_f, c_u_i, convex_constraints, nx, nu, np, N, w_tr, w_tr_s, w_vse, w_m)
%CVX_PIPG Solve same problem form as PIPG does but using CVX
%   Used to validate solution of PIPG algorithm implementations.
% Should probably return solver time too

cvx_begin
    variable x(nx, N)
    variable xi(nx, N)
    variable u(nu, N)
    variable s(np, 1)
    minimize( obj_function(x, xi, u, s, x_ref, u_ref, s_ref, w_tr, w_tr_s, w_vse, w_m) )
    subject to
        % Dynamics constraints
        for k = 1 : (N - 1)
            x(:, k + 1) == A_k(:, :, k) * x(:, k) + B_k_minus(:, :, k) * u(:, k) + B_k_plus(:, :, k) * u(:, k + 1) + S_k(:, k) * s + d_k(:, k);
        end

        % Cone constraints
        for cc = 1:numel(convex_constraints.x)
            convex_constraints.x{cc}(xi(:, 2 : (N - 1))) <= 0;
        end
        for cc = 1:numel(convex_constraints.u)
            convex_constraints.u{cc}(u(:, 2 : N)) <= 0;
        end
        for cc = 1:numel(convex_constraints.s)
            convex_constraints.s{cc}(s) <= 0;
        end

        %c_x_i(x(:, 1)) == 0;
        c_x_i(xi(:, 1)) == 0;
        c_x_f(xi(:, N)) == 0;
        c_u_i(u(:, 1)) == 0;
cvx_end

x_star = x;
xi_star = xi;
u_star = u;
s_star = s;

sol_info.status = cvx_status;
sol_info.J = cvx_optval;

end

function [val] = obj_function(x, xi, u, s, x_ref, u_ref, s_ref, w_tr, w_tr_s, w_vse, w_m)
    val = 0;
    for k = 1 : size(x, 2)
        val = val + w_tr * (x(:, k) - x_ref(:, k))' * (x(:, k) - x_ref(:, k));
    end
    for k = 1 : size(u, 2) 
        val = val + w_tr * (u(:, k) - u_ref(:, k))' * (u(:, k) - u_ref(:, k));
    end
    for k = 1 : numel(s)
        val = val + w_tr_s * (s(k) - s_ref(k))' * (s(k) - s_ref(k));
    end
    for k = 1 : size(xi, 2)
        val = val + w_vse * (x(:, k) - xi(:, k))' * (x(:, k) - xi(:, k));
    end

    val = val / 2;

    val = val - w_m * x(7, end);
end
