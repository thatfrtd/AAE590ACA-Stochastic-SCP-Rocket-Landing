%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 14 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 
load("3DoF_deterministic_Workspace.mat");
clear tri;
%% Initialize

%% Stochastic Optimization Parameters
nu = 2;
nr = 2;
n_sigma_99 = sigma_mag_confidence(1e-2, nu);
n_sigma_99p9 = sigma_mag_confidence(1e-3, nu);

n_probit_99p9 = norminv(1 - 1e-3);

%Initial Covariance values
%% Define Stochastic Elements
% Initial estimated state
sigma_xhat0 = [10e-3; ... % r_x
            20e-3; ... % r_y
            1e-3; ... % v_x
            1e-3; ... % v_y
            10e-3; ... % theta
            1e-3; ... % w
            1e-4]; % mass
Phat0 = diag(sigma_xhat0 .^ 2);

% Disturbance
sigma_accelx = 0.5e-3;
sigma_accely = 0.1e-3;
sigma_theta = 0.3e-3;
sigma_ang_vel = 0.2e-3;
sigma_m = 1e-7;

delta_t = 1e-1;
G = @(t, x, u, p) sqrt(delta_t) * [zeros([2, 3]); ... % velocity
                   [sigma_accelx; sigma_accely] .* eye([2, 3]); ... % acceleration
                   zeros([1, 3]); ... % angular velocity
                   sigma_theta * [1, 1, 0];
                   [zeros([1, 2]), sigma_m]]; ... % mass flow 

% Initial state estimation error
sigma_xtilde0 = 2 * [1e-4; ... % r_x
            5e-4; ... % r_y
            1e-4; ... % v_x
            1e-4; ... % v_y
            1e-4; ... % theta
            1e-5; ... % w
            1e-5]; % mass
Ptilde0 = diag(sigma_xtilde0 .^ 2);

% Final state covariance
% Final state
sigma_xf = [1e-3; ... % r_x
            1e-4; ... % r_y
            1e-4; ... % v_x
            1e-4; ... % v_y
            1e-3; ... % theta
            1e-4; ... % w
            1e-5]; % mass
Pf = diag(sigma_xf .^ 2);

% PTR algorithm parameters are defined in Deterministic


%% Get Dynamics
f = @(t, x, u, p) SymDynamics3DoF_mass_noumag(t, x, u, vehicle.m, vehicle.L, vehicle.I(2), vehicle.alpha);

%% Measurement Model
% Measurement model (identity with noise)
g_0_stds = [10e-5; ... % r_x
            50e-5; ... % r_y
            5e-5; ... % v_x
            5e-5; ... % v_y
            10e-3; ... % theta
            5e-5; ... % w
            1e-5]; ... % mass

f_0 = @(t, x, u, p) x;
g_0 = @(t, x, u, p) diag(g_0_stds);

% Gain 
K_k = -3*repmat([0, 1, 0, 0, 0, 0, 0; 0, 0, 0, 0, -1, 0, 0], 1, 1, prob_3DoF.Nu);

%% Specify Constraints
% Convex state path constraints
glideslope_angle_max = deg2rad(80);
h_glideslope = calculate_glideslope_offset(sigma_xf(1:2) * norminv(1 - 1e-3 / 2), glideslope_angle_max);
glideslope_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) (norm(x(1:2) + [0; h_glideslope]) + 0*sigma_mag_confidence(1e-3, nr) * norm(X_k_ref(1:2, (tri(k - 1, nx) + 1):tri(k, nx)))) * cos(glideslope_angle_max) - (x(2) + h_glideslope - norminv(1 - 1e-3) * norm(X_k_ref(2, (tri(k - 1, nx) + 1):tri(k, nx))));
state_nonconvex_constraints = {glideslope_constraint};
state_convex_constraints = {};
control_convex_constraints ={};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];
% Nonconvex control constraints
lcvx_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) norm(u(1:2)) + sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k) - thrust_magnitude_bound(S_k_ref, u_ref, k, t_k, T_max, m_0, alpha, nu);
min_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) T_min / m_0 - (norm(u(1:2)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k));
control_nonconvex_constraints = {lcvx_thrust_constraint};
control_nonconvex_constraints = {};
% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];
%% Specify Objective

% Size Definitions
nu = prob_3DoF.n.u;
nx = prob_3DoF.n.x;

% Objective
min_fuel_angular_velocity_objective = @(x, u, p, X_K, S_K) sum(u(3, :) / T_max + x(6, 1:Nu) .^ 2) * delta_t;
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p) sum(u(3, :)) * delta_t;
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p) sum((u(3, 1:(end - 1)) + u(3, 2:end)) / 2) * delta_t;
end

stochastic_min_fuel_objective = @(x, u, p, X_k, S_k) einsum(@(k) norm(u(1:2, k), 2) + sigma_mag_confidence(1e-2, nu) * norm(S_k(:, (tri(k - 1, nx) + 1):tri(k, nx)), 2), 1:Nu) * (t_k(2) - t_k(1));

%% Create Guess
sl_guess = guess_3DoF(x_0(1:6), x_f + [0; 0; 0; 0; 0; 0], N, Nu, delta_t, vehicle);
if u_hold == "ZOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [cumsum(sl_guess.u(3, :) * delta_t), sum(sl_guess.u(3, :)) * delta_t]];
elseif u_hold == "FOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * cumsum(sl_guess.u(3, :) * delta_t)]
end

%CasADi_sol = CasADi_solve_mass(x_0, sl_guess.x, sl_guess.u, vehicle, N, delta_t, glideslope_angle_max);

guess = CasADi_sol;
guess.u = guess.u(1:2, :);
if u_hold == "ZOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "previous","extrap")';
elseif u_hold == "FOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "linear","extrap")';
end
guess.p = sl_guess.p;


% figure
% plot_3DoF_trajectory(t_k, sl_guess.x, sl_guess.u, glideslope_angle_max, gimbal_max, T_min, T_max)
% 
% figure
% plot_3DoF_time_histories(t_k, sl_guess.x, sl_guess.u)

% figure
% plot_3DoF_trajectory(t_k, guess.x, guess.u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)
% 
% figure
% plot_3DoF_time_histories(t_k, guess.x, guess.u)

%% Construct Problem Object
mod_ptr_sol = ptr_sol;
mod_ptr_sol.x = ptr_sol.x(:, :, ptr_sol.converged_i);
mod_ptr_sol.u = ptr_sol.u(1:2, :, ptr_sol.converged_i);
mod_ptr_sol.p = ptr_sol.p(:, ptr_sol.converged_i);
ptr_ops.iter_min = 2;
stoch_terminal_bc = @(x, p) [x(1:6) - prob_3DoF.xf; 0]; 

prob_3DoF = StochasticProblem(x_0, x_f, Phat0, Ptilde0, Pf*10, N, u_hold, tf, f, G, f_0, g_0, mod_ptr_sol, convex_constraints, stochastic_min_fuel_objective, scale = scale, nonconvex_constraints = nonconvex_constraints,  terminal_bc = stoch_terminal_bc);

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
Delta = calculate_defect(prob_3DoF, prob_3DoF.guess.x, prob_3DoF.guess.u, prob_3DoF.guess.p);
norm(Delta)

%% Test Discretization on Initial Guess

%guess.u(1, :) = T_min * 1.5;
%guess.u(2, :) = T_min / 2000;
%guess.u(3, :) = vecnorm(guess.u(1:2, :));

[prob_3DoF, Delta_disc] = prob_3DoF.discretize(prob_3DoF.guess.x, prob_3DoF.guess.u, prob_3DoF.guess.p);

x_disc = prob_3DoF.disc_prop(guess.x, guess.u, guess.p, K_k);

[t_cont, x_cont, u_cont] = prob_3DoF.cont_prop(guess.x, guess.u, guess.p, K_k);

% figure
% comparison_plot_3DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")
% 
% figure
% comparison_plot_3DoF_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess.u, u_cont, guess.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%% Solve Problem with PTR
ptr_sol = Stochastic_ptr(prob_3DoF, ptr_ops);

%%
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [prob_3DoF.objective(prob_3DoF.guess.x, prob_3DoF.guess.u, prob_3DoF.guess.p), [ptr_sol.info.J]]); hold on
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

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));

figure
plot_3DoF_trajectory(t_k, ptr_sol.x(:, :, i), ptr_sol.u(:, :, i), glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
comparison_plot_3DoF_trajectory({guess.x, x_cont_sol, ptr_sol.x(:, :, i), CasADi_sol.x}, ["Guess", "Continuous Propagation", "Solution Output", "CasADi"], glideslope_angle_max, linestyle = [":", "-", "--", "-"], title = "Continuous vs Discrete Propagation of Solution")

figure
comparison_plot_3DoF_time_histories({t_k, t_cont_sol, t_k}, {guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, {guess.u, u_cont_sol, ptr_sol.u(:, :, i)}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")

