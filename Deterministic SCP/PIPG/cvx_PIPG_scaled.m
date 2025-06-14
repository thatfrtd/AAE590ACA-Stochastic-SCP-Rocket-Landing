function [z_star, sol_info] = cvx_PIPG_scaled(Q, q, c, H, h, convex_constraints, convex_constraints_eq, S_z, c_z)
%CVX_PIPG Solve same problem form as PIPG does but using CVX
%   Used to validate solution of PIPG algorithm implementations.
% Should probably return solver time too

cvx_begin
    variable z(size(q))
    minimize( 1/2 * (S_z * z + c_z)' * Q * (S_z * z + c_z) + q' * (S_z * z + c_z) + c )
    subject to
        % Affine constraints
        H * (S_z * z + c_z) - h == 0;

        % Cone constraints
        for cc = 1:numel(convex_constraints)
            convex_constraints{cc}(S_z * z + c_z) <= 0;
        end
        for cq = 1:numel(convex_constraints_eq)
            convex_constraints_eq{cq}(S_z * z + c_z) == 0;
        end
cvx_end

z_star = S_z * z + c_z;

sol_info.status = cvx_status;

end

