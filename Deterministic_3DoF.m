%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Physical Parameters
g = [-3.7114e-3; 0; 0]; % [km / s2]
alpha = 0.5086; % [s / km]
T_min = 4.97; % [kg km / s2]
T_max = 13.26; % [kg km / s2]
I = 15000 * (1e-3) ^ 2; % [kg m2]

% Problem Parameters
tf = 80; % [s]
N = 80; % []
r_0 = [1.5; 2.0]; % [km]
v_0 = [-0.075; 0.1]; % [km / s]
m_0 = 2000; % [kg]
gimbal_max = deg2rad(8); % [rad]
glideslope_angle_max = deg2rad(45); % [rad]

x_0 = [r_0; v_0];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

% PTR algorithm parameters
ptr_ops.iter_max = 15;
ptr_ops.Delta_min = 5e-3;
ptr_ops.w_vc = 1e1;
ptr_ops.w_tr = 1e-1;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.delta_tol = 1e-5;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

%% Specify Constraints
% Convex state path constraints
glideslope_constraint = @(x, u, p) norms(x, 2, 2) - x(2, :) / cos(glideslope_angle_max);
state_convex_constraints = {glideslope_constraint};

% Convex control constraints
max_thrust_constraint = @(x, u, p) u(3, :) - T_max;
min_thrust_constraint = @(x, u, p) T_min - u(3, :);
max_gimbal_constraint = @(x, u, p) norms(u, 2, 2) - u(1, :) / cos(gimbal_max);
control_convex_constraints = {max_thrust_constraint, min_thrust_constraint, max_gimbal_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

%% Specify Objective
min_fuel_objective = @(x, u, p) sum(u(3,:)) * delta_t;

%% Create Guess
guess = straight_line_interpolate_3DoF_guess(x0, xf, N);

%% Construct Problem Object
prob_3DoF = DeterministicProblem(x0, N, "ZOH", tf, f, guess, convex_constraints, min_fuel_objective);

%% Solve Problem with PTR
ptr_sol = ptr(prob_3DoF, ptr_ops);