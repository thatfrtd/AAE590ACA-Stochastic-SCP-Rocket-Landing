function [CasADi_sol] = CasADi_solve_6DoFdqg(x_initial, x_final, x_guess, u_guess, g, vehicle, N, delta_t, alpha_ME, alpha_RCS, tau_max)
%CASADI_SOLVE Summary of this function goes here
%   Detailed explanation goes here
% TODO:
% - Add FOH (don't forget to update objective)

% Define the optimization problem
opti = casadi.Opti();

max_iter = 400;

% Generate the array of state and control vectors

% States: x, y, x_dot, y_dot, theta, theta_dot
x = opti.variable(15, N);  % N x 6 matrix
% Controls: thrust (percent), thrust_angle (rad)
u = opti.variable(6, N - 1);

% Initial and final conditions
opti.subject_to(x(2:end, N) == x_final); % Final state

% Cost function to minimize effort and angular velocity
cost = -x(1, N);
opti.minimize(cost);

%constraints
for i = 1:(N-1)
    % Current state
    x_current = x(:, i);
    u_current = u(:, i);
    
    % Define the state derivatives
    
    % Runge Kutta 4 integration for dynamics constraints
    x_next = rk4(@(t, x, u, p) SymDynamicsDualQuat6DoF(x, u, 1, g, vehicle.L, vehicle.I, alpha_ME, alpha_RCS), delta_t * (i - 1), x_current, @(t) u_current, 0, delta_t);
    %x_next = x_current + delta_t * SymDynamics3DoF(delta_t * (i - 1), x_current, u_current, vehicle.m, vehicle.L, vehicle.I(2));

    % Impose the dynamics constraint
    opti.subject_to(x(:, i+1) == x_next);
end

% Thrust bounds
opti.subject_to(u(1, :) >= vehicle.min_thrust);
opti.subject_to(u(1, :) <= vehicle.max_thrust);
opti.subject_to(u(2, :) <= vehicle.max_gimbal);
opti.subject_to(u(2, :) >= 0);

opti.subject_to(abs(u(4:6)) <= tau_max)

% Solver options
%p_opts = struct('expand',true,'error_on_fail',false,'verbose',false);
%s_opts = struct('max_iter',max_iter);
%opti.solver('ipopt', p_opts, s_opts);
opti.solver('ipopt');

p = opti.parameter(15, 1);
opti.set_value(p, x_initial);
opti.subject_to(x(:, 1) == p); % Initial state

opti.set_initial(x, x_guess);
opti.set_initial(u, u_guess(:, 1:(N - 1)));

% Solve the optimization problem
sol = opti.solve();

CasADi_sol.u = sol.value(u);
CasADi_sol.x = sol.value(x);
CasADi_sol.objective = sol.value(cost);
end

