function [g_eq, g_ineq, g_eq_relax, g_ineq_relax] = eval_nonconvex_constraints(prob, x, u, p, x_ref, u_ref, p_ref)
%EVAL_NONCONVEX_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

t_k = prob.t_k;

g_ineq = [];
g_ineq_relax = [];

ncvx_ineq_index = 0;
for k = 1:prob.Nu
    for nc = 1:prob.n.ncvx
        nc_k = prob.nonconvex_constraints{nc}{1};
        if ismember(k, nc_k)
            ncvx_ineq_index = ncvx_ineq_index + 1;
    
            nc_ineq_func = prob.nonconvex_constraints{nc}{3};
            nc_ineq_relax_func = prob.nonconvex_constraints{nc}{2};
    
            g_ineq(ncvx_ineq_index) = nc_ineq_func(t_k(k), prob.unscale_x(x(:, k)), prob.unscale_u(u(:, k)), prob.unscale_p(p), k);
            g_ineq_relax(ncvx_ineq_index) = nc_ineq_relax_func(t_k(k), ...
                                                               prob.unscale_x(x(:, k)), prob.unscale_u(u(:, k)), prob.unscale_p(p), ...
                                                               prob.unscale_x(x_ref), prob.unscale_u(u_ref), prob.unscale_p(p_ref), ...
                                                               k);
        end
    end
end

% NEED TO ADD NON DIFFERENTIAL NONCONVEX INEQUALITY CONSTRAINTS
g_eq = [];
g_eq_relax = [];

end