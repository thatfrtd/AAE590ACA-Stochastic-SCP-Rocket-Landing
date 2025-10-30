function [convex_constraints, nonconvex_constraints, original_constraints] = construct_GA_constraints_fixedtf(transition_k, t_k, mu, r_p_min, x_p, nr, v_scale_ref)
%CONSTRUCT_GA_CONSTRAINTS Summary of this function goes here
%   Detailed explanation goes here

r_i = 1:nr;
r_select = [eye(nr), zeros([nr, nr])];

v_i = nr + (1:nr);
v_select = [zeros([nr, nr]), eye(nr)];

% Phase convex state constraints
GA_position_constraint = {@(t, x, u, p) x(r_i, 1) - r_select * x_p(t_k(transition_k)), "=="};
GA_pos_mass_constraint = {@(t, x, u, p) x([r_i, 5], 1) - x([r_i, 5], 2), "=="};

convex_constraints = {GA_position_constraint, GA_pos_mass_constraint};

% Phase nonconvex state constraints
v_k_sym = sym("v_k", [nr, 1]);
v_kp1_sym = sym("v_kp1", [nr, 1]);

%GA_vinfinity_constraint = @(v_k, v_kp1, v_k_ref, v_kp1_ref) (v_kp1_ref - v_select * x_p(t_k(transition_k)))' * (v_kp1 - v_select * x_p(t_k(transition_k))) ...
%                                                          - (v_k_ref - v_select * x_p(t_k(transition_k)))' * (v_k - v_select * x_p(t_k(transition_k)));
GA_vinfinity_constraint = @(v_k, v_kp1) sum_square(v_kp1 - v_select * x_p(t_k(transition_k))) - sum_square(v_k - v_select * x_p(t_k(transition_k)));
%GA_vinfinity_constraint = @(v_k, v_kp1) dnorm(v_kp1 - v_select * x_p(t_k(transition_k))) - dnorm(v_k - v_select * x_p(t_k(transition_k)));


GA_vinfinity_constraint_partial_v_k = matlabFunction(jacobian(GA_vinfinity_constraint(v_k_sym, v_kp1_sym), [v_k_sym]),"Vars", [{v_k_sym}, {v_kp1_sym}]);
GA_vinfinity_constraint_partial_v_kp1 = matlabFunction(jacobian(GA_vinfinity_constraint(v_k_sym, v_kp1_sym), [v_kp1_sym]),"Vars", [{v_k_sym}, {v_kp1_sym}]);

GA_vinfinity_constraint_linearized = {@(t, x, u, p, x_ref, u_ref, p_ref) GA_vinfinity_constraint(x_ref(v_i, 1), x_ref(v_i, 2))  ... 
                                                                       + GA_vinfinity_constraint_partial_v_k(x_ref(v_i, 1), x_ref(v_i, 2)) * (x(v_i, 1) - x_ref(v_i, 1)) ...
                                                                       + GA_vinfinity_constraint_partial_v_kp1(x_ref(v_i, 1), x_ref(v_i, 2)) * (x(v_i, 2) - x_ref(v_i, 2)), "=="};
%GA_vinfinity_constraint_linearized = {@(t, x, u, p, x_ref, u_ref, p_ref) GA_vinfinity_constraint(x(v_i, 1), x(v_i, 2), x_ref(v_i, 1), x_ref(v_i, 2)), "=="}

eta = @(v_k, v_kp1, v_p) (v_k - v_p)' * (v_kp1 - v_p) / (dnorm(v_k - v_p) * dnorm(v_kp1 - v_p));
GA_minperiapsis_constraint = @(v_k, v_kp1) dnorm(v_k - v_select * x_p(t_k(transition_k))) * v_scale_ref - sqrt(mu / r_p_min * (1 / sin(1/2 * acos(eta(v_k, v_kp1, v_select * x_p(t_k(transition_k))))) - 1));

GA_minperiapsis_constraint_partial_v_k = matlabFunction(jacobian(GA_minperiapsis_constraint(v_k_sym, v_kp1_sym), [v_k_sym]),"Vars", [{v_k_sym}, {v_kp1_sym}]);
GA_minperiapsis_constraint_partial_v_kp1 = matlabFunction(jacobian(GA_minperiapsis_constraint(v_k_sym, v_kp1_sym), [v_kp1_sym]),"Vars", [{v_k_sym}, {v_kp1_sym}]);

GA_minperiapsis_constraint_linearized = {@(t, x, u, p, x_ref, u_ref, p_ref) GA_minperiapsis_constraint(x_ref(v_i, 1), x_ref(v_i, 2)) ...
                                                                          + GA_minperiapsis_constraint_partial_v_k(x_ref(v_i, 1), x_ref(v_i, 2)) * (x(v_i, 1) - x_ref(v_i, 1)) ...
                                                                          + GA_minperiapsis_constraint_partial_v_kp1(x_ref(v_i, 1), x_ref(v_i, 2)) * (x(v_i, 2) - x_ref(v_i, 2)), "<="};

nonconvex_constraints = {GA_vinfinity_constraint_linearized, GA_minperiapsis_constraint_linearized};

%nonconvex_constraints = [nonconvex_constraints, {{@(t, x, u, p, x_ref, u_ref, p_ref) GA_position_constraint{1}(t, x, u, p), "=="}, {@(t, x, u, p, x_ref, u_ref, p_ref) GA_pos_mass_constraint{1}(t, x, u, p), "=="}}];

%convex_constraints = {};

eval_constraints = @(t, x, u, p) [GA_position_constraint{1}(t, [x(:, transition_k), x(:, transition_k + 1)], u, p); GA_pos_mass_constraint{1}(t, [x(:, transition_k), x(:, transition_k + 1)], u, p); GA_vinfinity_constraint(x(v_i, transition_k), x(v_i, transition_k + 1)); GA_minperiapsis_constraint(x(v_i, transition_k), x(v_i, transition_k + 1))];
constraint_types = {GA_position_constraint{2}, GA_pos_mass_constraint{2}, GA_vinfinity_constraint_linearized{2}, GA_minperiapsis_constraint_linearized{2}};

original_constraints = {eval_constraints, constraint_types};

end

function [val] = dnorm(vec)
    val = sqrt(sum_square(vec));
end