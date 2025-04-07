%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Vehicle Parameters
alpha = 0.5086; % [s / km]
T_min = 4.97; % [kg km / s2]
T_max = 13.26; % [kg km / s2]
I = 15000 * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 3e-3; % [km] Distance from CoM to nozzle
m_dry = 2000; % [kg]
m_wet = 600; % [kg]
m_0 = m_dry + m_wet;
gimbal_max = deg2rad(8); % [rad]

vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I);

% Problem Parameters
tf = 80; % [s]
N = 80; % []
r_0 = [1.5; 2.0]; % [km]
v_0 = [-0.0185; -0.0247]; % [km / s]
theta_0 = deg2rad(90); % [rad]
w_0 = deg2rad(0); % [rad / s]
glideslope_angle_max = deg2rad(45); % [rad]

x_0 = [r_0; v_0; theta_0; w_0];
x_f = [zeros(2, 1); zeros(2, 1); pi / 2; 0];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

% PTR algorithm parameters
ptr_ops.iter_max = 4;
ptr_ops.Delta_min = 5e-3;
ptr_ops.w_vc = 1e1;
ptr_ops.w_tr = 1e-1 * ones(1, Nu);
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = true;
ptr_ops.delta_tol = 1e-5;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

%% Get Dynamics
f = @(t, x, u, p) SymDynamics3DoF(t, x, u, vehicle.m, vehicle.L, vehicle.I(2));

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
guess = guess_3DoF(x_0, x_f, N, Nu, delta_t, vehicle);

% figure
% plot_3DoF_trajectory(t_k, guess.x, guess.u, glideslope_angle_max, gimbal_max, T_min, T_max)
% 
% figure
% plot_3DoF_time_histories(t_k, guess.x, guess.u)

%% Construct Problem Object
prob_3DoF = DeterministicProblem(x_0, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective);

%% Test Scaling
% guess_scaled.x = prob_3DoF.scale_x(guess.x);
% guess_scaled.u = prob_3DoF.scale_u(guess.u);
% guess_scaled.p = prob_3DoF.scale_p(guess.p);
% 
% figure
% plot_3DoF_trajectory(t_k, guess_scaled.x, guess_scaled.u, glideslope_angle_max, gimbal_max, 0, 0)
% 
% figure
% plot_3DoF_time_histories(t_k, guess_scaled.x, guess_scaled.u)

%%
Delta = calculate_defect(prob_3DoF, guess.x, guess.u, guess.p);

%% Test Discretization on Initial Guess

guess.u(1, :) = T_min * 1.5;
%guess.u(2, :) = T_min / 2000;
guess.u(3, :) = vecnorm(guess.u(1:2, :));

prob_3DoF = prob_3DoF.discretize(guess.x, guess.u, guess.p);

x_disc = prob_3DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_3DoF.cont_prop(guess.u, guess.p);

figure
comparison_plot_3DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

figure
comparison_plot_3DoF_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess.u, u_cont, guess.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%% Solve Problem with PTR
ptr_sol = ptr(prob_3DoF, ptr_ops);

%%
tiledlayout(1, 3)

nexttile
plot(0:ptr_ops.iter_max, ptr_sol.objective);
title("Objective vs Iteration")

nexttile
plot(ptr_sol.delta_xp)
title("Stopping Criteria vs Iteration")

nexttile
plot(0:ptr_ops.iter_max, vecnorm(ptr_sol.Delta))
title("Defect Norm vs Iteration")