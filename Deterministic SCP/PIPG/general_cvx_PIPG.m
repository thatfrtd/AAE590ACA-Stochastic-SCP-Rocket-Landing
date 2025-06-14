function [z_star, sol_info] = general_cvx_PIPG(z_ref, H, h, convex_constraints, convex_constraints_eq, nx, nu, np, N, w_tr, w_tr_s, w_vse, w_m)
%CVX_PIPG Solve same problem form as PIPG does but using CVX
%   Used to validate solution of PIPG algorithm implementations.
% Should probably return solver time too

x_i_0 = 1;
x_i_f = N * nx;
x_i = x_i_0 : x_i_f;

xi_i_0 = x_i_f + 1;
xi_i_f = x_i_f + N * nx;
xi_i = xi_i_0 : xi_i_f;

u_i_0 = xi_i_f + 1;
u_i_f = xi_i_f + N * nu;
u_i = u_i_0 : u_i_f;

s_i_0 = u_i_f + 1;
s_i_f = u_i_f + np;
s_i = s_i_0 : s_i_f;

cvx_begin
    variable z(size(z_ref))
    minimize( obj_function(z, z_ref, x_i, xi_i, u_i, s_i, w_tr, w_tr_s, w_vse, w_m) )
    subject to
        % Affine constraints
        H * z - h == 0;

        % Cone constraints
        for cc = 1:numel(convex_constraints)
            convex_constraints{cc}(z) <= 0;
        end
        for cq = 1:numel(convex_constraints_eq)
            convex_constraints_eq{cq}(z) == 0;
        end
cvx_end

z_star = z;

sol_info.status = cvx_status;

end

function [val] = obj_function(z, z_ref, x_i, xi_i, u_i, s_i, w_tr, w_tr_s, w_vse, w_m)
    val = 0;
    for k = 1 : numel(x_i)
        val = val + w_tr * (z(x_i(k)) - z_ref(x_i(k))) * (z(x_i(k)) - z_ref(x_i(k)));
    end
    for k = 1 : numel(u_i) 
        val = val + w_tr * (z(u_i(k)) - z_ref(u_i(k))) * (z(u_i(k)) - z_ref(u_i(k)));
    end
    for k = 1 : numel(s_i)
        val = val + w_tr_s * (z(s_i(k)) - z_ref(s_i(k))) * (z(s_i(k)) - z_ref(s_i(k)));
    end
    for k = 1 : numel(xi_i)
        val = val + w_vse * (z(x_i(k)) - z(xi_i(k))) * (z(x_i(k)) - z(xi_i(k)));
    end

    val = val / 2;

    val = val - w_m * z(x_i(end));
end
