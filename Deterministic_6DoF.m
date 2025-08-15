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
m_dry = 2000; % [kg]
m_wet = 600; % [kg]
m_0 = m_dry + m_wet;
T_max = 3 * m_0 * 9.81e-3; % [kg km / s2]
T_min = 0.55 * T_max; % [kg km / s2]
I = [50000; 150000; 150000] * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 3e-3; % [km] Distance from CoM to nozzle
gimbal_max = deg2rad(8); % [rad]

vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I);

% Problem Parameters
tf = 30; % [s]
N = 50; % []
r_0 = [0; 0; 4.6]; % [km]
theta_0 = deg2rad(120); % [rad]
v_0 = [1 0; 0 0; 0 1] * make_R2(deg2rad(-60)) * [0.306; 0]; % [km / s]
[theta1_0, theta2_0, theta3_0] = dcm2angle(angle2dcm(0, deg2rad(-theta_0), 0,"ZYX"), "XYX"); % [rad]
theta_0 = [theta1_0; theta2_0; theta3_0];
w_0 = deg2rad([0; 0; 0]); % [rad / s]
glideslope_angle_max = deg2rad(60); % [rad]

[theta1_f, theta2_f, theta3_f] = dcm2angle(angle2dcm(0, deg2rad(-90), 0,"ZYX"), "XYX"); % [rad]
theta_f = [theta1_f; theta2_f; theta3_f];

x_0 = [r_0; v_0; theta_0; w_0];
x_f = [zeros(3, 1); zeros(3, 1); theta_f; zeros(3, 1)];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "FOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

% PTR algorithm parameters
ptr_ops.iter_max = 20;
ptr_ops.Delta_min = 1e-2;
ptr_ops.w_vc = 1e3;
ptr_ops.w_tr = ones(1, Nu) * 1e2;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = true;
ptr_ops.delta_tol = 2e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = true;

scale_hint.x_max = [max(r_0) * ones([3, 1]); max(v_0) * ones([3, 1]); pi * ones([3, 1]); max(w_0) * ones([3, 1])];
scale_hint.x_min = [-max(r_0) * ones([3, 1]); -max(v_0) * ones([3, 1]); -pi * ones([3, 1]); -max(w_0) * ones([3, 1])];
scale_hint.u_max = [T_max; sin(gimbal_max) * ones([2, 1]); T_max; pi / 4];
scale_hint.u_min = [T_min; -sin(gimbal_max) * ones([2, 1]); T_min; -pi / 4];
scale_hint.p_max = 80;
scale_hint.p_min = 40;

%% Get Dynamics
f = @(t, x, u, p) SymDynamicsEuler6DoF(x, u, vehicle.m, vehicle.L, vehicle.I);

%% Specify Constraints
% Convex state path constraints
glideslope_constraint = {1:N, @(t, x, u, p) norm(x(1:3)) - x(3) / cos(glideslope_angle_max)};
state_convex_constraints = {glideslope_constraint};

% Convex control constraints
max_thrust_constraint = {1:N, @(t, x, u, p) u(4) - T_max};
min_thrust_constraint = {1:N, @(t, x, u, p) T_min - u(4)};
max_gimbal_constraint = {1:N, @(t, x, u, p) u(4) - u(1) / cos(gimbal_max)};
lcvx_thrust_constraint = {1:N, @(t, x, u, p) norm(u(1:3))- u(4)}; 
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
sl_guess = guess_6DoF(x_0, x_f, N, Nu, delta_t, vehicle);

%CasADi_sol = CasADi_solve_6DoF(x_0, x_f, sl_guess.x, sl_guess.u, vehicle, N, delta_t, glideslope_angle_max);

sl_guess.u = sl_guess.u * T_max;

guess = sl_guess;
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

%cfigure
%plot_6DoF_trajectory(t_k, guess.x, guess.u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
plot_6DoF_time_histories(t_k, guess.x, guess.u)

%% Construct Problem Object
prob_6DoF = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective, scale = scale, scale_hint = scale_hint);

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
Delta = calculate_defect(prob_6DoF, guess.x, guess.u, guess.p);
norm(Delta)

%% Test Discretization on Initial Guess

[prob_6DoF, Delta_disc] = prob_6DoF.discretize(guess.x, guess.u, guess.p);

x_disc = prob_6DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_6DoF.cont_prop(guess.u, guess.p);

figure
comparison_plot_6DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

figure
comparison_plot_6DoF_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess.u, u_cont, guess.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%% Solve Problem with PTR
ptr_ops.w_vc = 7e3;
ptr_ops.w_tr = ones(1, Nu) * 5e0;
ptr_ops.iter_min = 2;
ptr_sol = ptr(prob_6DoF, ptr_ops);

%%
scvxstar_ops.D_x = 1;
scvxstar_ops.D_u = 1;
scvxstar_ops.D_p = 1;
scvxstar_ops.opt_tol = ptr_ops.delta_tol;
scvxstar_ops.feas_tol = 5e-4;%ptr_ops.Delta_min;
scvxstar_ops.eta_0 = 1;
scvxstar_ops.eta_1 = 0.5;
scvxstar_ops.eta_2 = 0.1;
scvxstar_ops.alpha_1 = 2;
scvxstar_ops.alpha_2 = 4;
scvxstar_ops.beta = 2;
scvxstar_ops.gamma = 0.95;
scvxstar_ops.w_0 = 1e1;
scvxstar_ops.w_max = 1e6;
scvxstar_ops.r_0 = 0.1;
scvxstar_ops.r_min = 1e-8;
scvxstar_ops.r_max = 1;
scvxstar_ops.tau = 1.1;
scvxstar_ops.iter_max = ptr_ops.iter_max*5;
scvxstar_ops.iter_min = 2;

scvxstar_sol = SCvx_star(prob_6DoF, scvxstar_ops, "CVX");


%%
if ~ptr_sol.converged
    ptr_sol.converged_i = ptr_ops.iter_max;
end

%%
% ptr_ops.iter_max = 30;
% ptr_ops.w_vse = 1e6;
% ptr_ops.w_tr = 3e3;
% ptr_ops.w_prime = 1e2;
% ptr_sol = ptr_virtual_state(prob_6DoF, ptr_ops, "CVX");

%%
nx = 13;
nu = 5;

Rxtheta1 = [1 0 0; 0 cos(theta1_0) sin(theta1_0); 0 -sin(theta1_0) cos(theta1_0)];
Rxtheta2 = [cos(theta2_0) 0 -sin(theta2_0); 0 1 0; sin(theta2_0) 0 cos(theta2_0)];
Rxtheta3 = [1 0 0; 0 cos(theta3_0) sin(theta3_0); 0 -sin(theta3_0) cos(theta3_0)];
C_be = Rxtheta3 * Rxtheta2 * Rxtheta1;

v_0_b = angle2dcm(theta1_0, theta2_0, theta3_0, "XYX") * v_0% - cross(w_0, r_0);

opt_time = t_k;
control_inputs = guess.u(:, :);
control_inputs(5, :) = deg2rad(1);
input_vector = 1:nu;
x_opt = guess.x(:, :);
state_vector = 1:nx;
r_0_6DoF = r_0;
v_0_6DoF = v_0;
v_0_b_6DoF = v_0_b;
[r0, p0, y0] = dcm2angle(angle2dcm(theta1_0, theta2_0, theta3_0,"XYX"), "ZYX");
rpy_0_6DoF = -[r0; p0; y0];
w_0_6DoF = w_0;

x_0_6DoF = [r_0_6DoF; v_0_6DoF; rpy_0_6DoF; w_0_6DoF];
I_matrix = diag(I);

%%
nx = 13;
nu = 5;

Rxtheta1 = [1 0 0; 0 cos(theta1_0) sin(theta1_0); 0 -sin(theta1_0) cos(theta1_0)];
Rxtheta2 = [cos(theta2_0) 0 -sin(theta2_0); 0 1 0; sin(theta2_0) 0 cos(theta2_0)];
Rxtheta3 = [1 0 0; 0 cos(theta3_0) sin(theta3_0); 0 -sin(theta3_0) cos(theta3_0)];
C_be = Rxtheta3 * Rxtheta2 * Rxtheta1;

v_0_b = angle2dcm(theta1_0, theta2_0, theta3_0, "XYX") * v_0% - cross(w_0, r_0);

opt_time = t_k;
control_inputs = ptr_sol.u(:, :, i);
input_vector = 1:nu;
x_opt = ptr_sol.x(:, :, i);
state_vector = 1:nx;
r_0_6DoF = r_0;
v_0_6DoF = v_0;
v_0_b_6DoF = v_0_b;
rpy_0_6DoF = dcm2angle(angle2dcm(theta1_0, theta2_0, theta3_0,"XYX"), "ZYX");
w_0_6DoF = w_0;

x_0_6DoF = [r_0_6DoF; v_0_6DoF; rpy_0_6DoF; w_0_6DoF];
I_matrix = diag(I);
% 
% T_e = zeros([2, N]);
% v_B = zeros([2, N]);
% 
% for j = 1:N
%     T_e(:, j) = make_R2(x_opt(5, j)) * control_inputs(1:2, j) .* exp(x_opt(7, j));
%     v_B(:, j) = make_R2(x_opt(5, j))' * x_opt(3:4, j) - eye(2, 3) * cross([0; 0; x_opt(6, j)], [x_opt(1:2, j); 0]);
% end



%%
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [prob_6DoF.objective(prob_6DoF.guess.x, prob_6DoF.guess.u, prob_6DoF.guess.p), [ptr_sol.info.J]]); hold on
yline(CasADi_sol.objective); hold off
legend("PTR Iterations", "CasADi Solution")
title("Objective vs Iteration")
grid on

nexttile
plot(ptr_sol.delta_xp)
title("Stopping Criteria vs Iteration")
grid on

nexttile
plot(0:ptr_sol.converged_i, vecnorm(ptr_sol.Delta(:, 1:(ptr_sol.converged_i + 1)), 2, 1))
title("Defect Norm vs Iteration")
grid on

%%
i = ptr_sol.converged_i;

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_6DoF.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));

figure
plot_6DoF_trajectory(t_k, ptr_sol.x(:, :, i), ptr_sol.u(:, :, i), glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
comparison_plot_6DoF_trajectory({guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, ["Guess", "Continuous Propagation", "Solution Output"], glideslope_angle_max, linestyle = [":", "-", "--", "-"], title = "Continuous vs Discrete Propagation of Solution")

figure
comparison_plot_6DoF_time_histories({t_k, t_cont_sol, t_k}, {guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, {guess.u, u_cont_sol, ptr_sol.u(:, :, i)}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")

