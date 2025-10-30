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
I = [50000; 150000; 150000] * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 3e-3; % [km] Distance from CoM to nozzle
m_dry = 2000; % [kg]
m_wet = 600; % [kg]
m_0 = m_dry + m_wet;
gimbal_max = deg2rad(8); % [rad]
g = 3.7114e-3; % [km / s2]

vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I);

% Problem Parameters
tf = 60; % [s]
N = 25; % []
r_0 = [0.5; -0.2; 2.0]; % [km]
v_0 = [0.0; 0.0; -0.0647] * 0.1; % [km / s]
q_0 = dcm2quat(angle2dcm(0, deg2rad(-90), 0,"ZYX"))'; % [rad]
w_0 = deg2rad([0; 10; 5]); % [rad / s]
glideslope_angle_max = deg2rad(60); % [rad]

q_f = dcm2quat(angle2dcm(0, deg2rad(-90), 0,"ZYX"))'; % [rad]

x_0 = [r_0; v_0; q_0; w_0];
x_f = [zeros(3, 1); zeros(3, 1); q_f; zeros(3, 1)];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

%% Get Dynamics
f = @(t, x, u, p) SymDynamicsEuler6DoF_q(x, u, [g, vehicle.m, vehicle.L, vehicle.I]);

%% Specify Constraints
% Convex state path constraints
glideslope_constraint = @(t, x, u, p) norm(x(1:3)) - x(3) / cos(glideslope_angle_max);
state_convex_constraints = {glideslope_constraint};

% Convex control constraints
max_thrust_constraint = @(t, x, u, p) u(4) - T_max;
min_thrust_constraint = @(t, x, u, p) T_min - u(4);
max_gimbal_constraint = @(t, x, u, p) u(4) - u(1) / cos(gimbal_max);
lcvx_thrust_constraint = @(t, x, u, p) norm(u(1:3))- u(4); 
control_convex_constraints = {min_thrust_constraint,max_gimbal_constraint,max_thrust_constraint,lcvx_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

%% Specify Objective
min_fuel_angular_velocity_objective = @(x, u, p) sum(u(3, :) / T_max + x(6, 1:Nu) .^ 2) * delta_t;
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p) sum(u(4, :)) * delta_t;
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p) sum((u(4, 1:(end - 1)) + u(4, 2:end)) / 2) * delta_t;
end

%% Create Guess
sl_guess = guess_6DoF_quat(x_0, x_f, N, Nu, delta_t, vehicle);

CasADi_sol = CasADi_solve_6DoF_Quaternion(x_0, x_f, sl_guess.x, sl_guess.u, vehicle, N, delta_t, glideslope_angle_max);

guess = CasADi_sol;
if u_hold == "ZOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "previous","extrap")';
elseif u_hold == "FOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "linear","extrap")';
end
guess.p = sl_guess.p;

%%
% figure
% plot_6DoF_trajectory(t_k, sl_guess.x, sl_guess.u, glideslope_angle_max, gimbal_max, T_min, T_max)
% 
% figure
% plot_6DoF_time_histories(t_k, sl_guess.x, sl_guess.u)

figure
plot_6DoF_trajectory_quaternion(t_k, guess.x, guess.u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
plot_6DoF_time_histories_quaternion(t_k, guess.x, guess.u)
