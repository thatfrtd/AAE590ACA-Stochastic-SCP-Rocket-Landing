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

vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I, alpha = alpha);

% Problem Parameters
tf = 90; % [s]
N = 20; % []
r_0 = [1.5; 0; 3.0]; % [km]
v_0 = [0.0385; 0; -0.0647]; % [km / s]
[theta1_0, theta2_0, theta3_0] = dcm2angle(angle2dcm(0, deg2rad(-90), 0,"ZYX"), "XYX"); % [rad]
theta_0 = [theta1_0; theta2_0; theta3_0];
w_0 = deg2rad([0; 0; 0]); % [rad / s]
glideslope_angle_max = deg2rad(60); % [rad]

[theta1_f, theta2_f, theta3_f] = dcm2angle(angle2dcm(0, deg2rad(-90), 0,"ZYX"), "XYX"); % [rad]
theta_f = [theta1_f; theta2_f; theta3_f];

x_0 = [r_0; v_0; theta_0; w_0; log(m_0)];
x_f = [zeros(3, 1); zeros(3, 1); theta_f; zeros(3, 1)];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

% PTR algorithm parameters
ptr_ops.iter_max = 100;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 2e3;
ptr_ops.w_tr = ones(1, Nu) * 5e0;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = false;

%% Get Dynamics
f = @(t, x, u, p) SymDynamicsEuler6DoF_convex(x, u, vehicle.L, vehicle.I, vehicle.alpha);

%% Specify Constraints
% Convex state path constraints
glideslope_constraint = @(t, x, u, p) norm(x(1:3)) - x(3) / cos(glideslope_angle_max);
state_convex_constraints = {glideslope_constraint};

z_lb = @(t) log(m_0 - alpha * T_max * t);
z_lb_k = z_lb(t_k);

% Convex control constraints
max_thrust_constraint = @(t, x, u, p) u(4) - T_max * exp(-z_lb(t)) * (1 - (x(13) - z_lb(t)));
min_thrust_constraint = @(t, x, u, p) T_min * exp(-x(13)) - u(4);
max_gimbal_constraint = @(t, x, u, p) u(4) - u(1) / cos(gimbal_max);
lcvx_thrust_constraint = @(t, x, u, p) norm(u(1:3))- u(4); 
control_convex_constraints = {min_thrust_constraint,max_gimbal_constraint,max_thrust_constraint,lcvx_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Terminal boundary conditions
terminal_bc = @(x, u, p) [x(1:12, :) - x_f; 0];

%% Specify Objective
min_fuel_angular_velocity_objective = @(x, u, p) sum(u(4, :) / T_max + x(6, 1:Nu) .^ 2) * delta_t;
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p) sum(u(4, :)) * delta_t;
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p) sum((u(4, 1:(end - 1)) + u(4, 2:end)) / 2) * delta_t;
end

%% Create Guess
sl_guess = guess_6DoF(x_0(1:12), x_f, N, Nu, delta_t, vehicle);
if u_hold == "ZOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [cumsum(sl_guess.u(4, :) * delta_t), sum(sl_guess.u(4, :)) * delta_t]];
elseif u_hold == "FOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * cumsum(sl_guess.u(4, :) * delta_t)];
end
sl_guess.x(13, :) = log(sl_guess.x(13, :));
sl_guess.u = sl_guess.u * T_max .* exp(-sl_guess.x(13, 1:Nu));
%%
%\CasADi_sol = CasADi_solve_6DoF_mass_convexified(x_0, x_f, guess.x, guess.u, vehicle, N, delta_t, glideslope_angle_max);
%%
guess = sl_guess;
if u_hold == "ZOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "previous","extrap")';
elseif u_hold == "FOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "linear","extrap")';
end
guess.p = sl_guess.p;
%guess.x = x_disc;
%%
% figure
% plot_6DoF_trajectory(t_k, sl_guess.x, sl_guess.u .* [exp(sl_guess.x(13, 1:Nu)); exp(sl_guess.x(13, 1:Nu)); exp(sl_guess.x(13, 1:Nu)); exp(sl_guess.x(13, 1:Nu)); ones([1, Nu])], glideslope_angle_max, gimbal_max, T_min, T_max)
% 
% figure
% plot_6DoF_time_histories(t_k, sl_guess.x, sl_guess.u)

figure
plot_6DoFc_trajectory(t_k, guess.x, guess.u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
plot_6DoFc_time_histories(t_k, guess.x, guess.u)

%% Construct Problem Object
prob_6DoF = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc);

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

guess.u(1, :) = T_min .* exp(-sl_guess.x(13, 1:Nu)) * 1.5;
guess.u(2, :) = T_min .* exp(-sl_guess.x(13, 1:Nu)) * -0.001;
guess.u(3, :) = T_min .* exp(-sl_guess.x(13, 1:Nu)) * -0.004;
guess.u(4, :) = vecnorm(guess.u(1:3, :));

[prob_6DoF, Delta_disc] = prob_6DoF.discretize(guess.x, guess.u, guess.p);

x_disc = prob_6DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_6DoF.cont_prop(guess.u, guess.p);

figure
comparison_plot_6DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

figure
comparison_plot_6DoF_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess.u, u_cont, guess.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

prob_6DoF = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc);

%% Solve Problem with PTR
ptr_sol = ptr(prob_6DoF, ptr_ops);

%%
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [prob_6DoF.objective(prob_6DoF.guess.x, prob_6DoF.guess.u, prob_6DoF.guess.p), [ptr_sol.info.J]]); hold on
%yline(CasADi_sol.objective); hold off
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
plot_6DoFc_trajectory(t_k, ptr_sol.x(:, :, i), ptr_sol.u(:, :, i), glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
comparison_plot_6DoF_trajectory({guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, ["Guess", "Continuous Propagation", "Solution Output"], glideslope_angle_max, linestyle = [":", "-", "--", "-"], title = "Continuous vs Discrete Propagation of Solution")

figure
comparison_plot_6DoFc_time_histories({t_k, t_cont_sol, t_k}, {guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, {guess.u, u_cont_sol, ptr_sol.u(:, :, i)}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")

