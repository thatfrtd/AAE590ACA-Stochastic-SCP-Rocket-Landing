function [z_star, sol_info] = cvx_PIPG(Q, q, c, H, h, convex_constraints, convex_constraints_eq)
%CVX_PIPG Solve same problem form as PIPG does but using CVX
%   Used to validate solution of PIPG algorithm implementations.
% Should probably return solver time too

cvx_begin
    variable z(size(q))
    minimize( 1/2 * z' * Q * z + q' * z + c )
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

