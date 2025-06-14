function [z_star, sol_info] = cvx_PIPG_pc(lambda, qhat, Hhat, h, L_inv, convex_constraints, convex_constraints_eq, S_z, c_z)
%CVX_PIPG Solve same problem form as PIPG does but using CVX
%   Used to validate solution of PIPG algorithm implementations.
% Should probably return solver time too

cvx_begin
    variable zhat(size(qhat))
    minimize( lambda / 2 * zhat' * zhat + qhat' * zhat )
    subject to
        % Affine constraints
        Hhat * zhat - h == 0;

        % Cone constraints
        for cc = 1:numel(convex_constraints)
            convex_constraints{cc}(S_z * L_inv * zhat + c_z) <= 0;
        end
        for cq = 1:numel(convex_constraints_eq)
            convex_constraints_eq{cq}(S_z * L_inv * zhat + c_z) == 0;
        end
cvx_end

z_star = S_z * L_inv * zhat + c_z;

sol_info.status = cvx_status;

end

