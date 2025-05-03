%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter, Atharva Awasthi
% Created On: 27 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 29 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear; 
%load("3DoF_deterministic_Workspace.mat");
%clear tri;
%% Initialize
Deterministic_6DoF_convexified
%% Stochastic Optimization Parameters
nu = 4;
nr = 3;
nx = 13;

n_probit_99p9 = norminv(1 - 1e-3);

%Initial Covariance values
%% Define Stochastic Elements
% Initial estimated state
sigma_xhat0 = [10e-3; ... % r_x
            10e-3; ... % r_y
            20e-3; ... % r_z
            1e-3; ... % v_x
            1e-3; ... % v_y
            1e-3; ... % v_z
            10e-3; ... % theta_1
            10e-3; ... % theta_2
            10e-3; ... % theta_3
            1e-3; ... % w_x
            1e-3; ... % w_y
            1e-3; ... % w_z
            1e-4]; % mass
Phat0 = diag(sigma_xhat0 .^ 2);

% Disturbance
sigma_accel = [0.5e-4; 0.5e-4; 0.1e-4];
sigma_ang_accel = [0.2e-5; 0.2e-5; 0.2e-5];
sigma_m = 1e-7;

delta_t = 1e-1;
G = @(t, x, u, p) sqrt(delta_t) * [zeros([3, 7]); ... % velocity
                                   sigma_accel .* eye([3, 7]); ... % acceleration
                                   zeros([3, 7]); ... % angular velocity
                                   sigma_ang_accel .* [zeros(3, 3), eye(3), zeros(3, 1)]; % angular acceleration
                                   sigma_m * [zeros([1, 6]), 1]]; ... % mass flow 

% Initial state estimation error
sigma_xtilde0 = [1e-4; ... % r_x
            1e-4; ... % r_y
            5e-4; ... % r_z
            1e-4; ... % v_x
            1e-4; ... % v_y
            1e-4; ... % v_z
            1e-4; ... % theta_1
            1e-4; ... % theta_2
            1e-4; ... % theta_3
            1e-5; ... % w_x
            1e-5; ... % w_y
            1e-5; ... % w_z
            1e-5]; % mass
Ptilde0 = diag(sigma_xtilde0 .^ 2);

% Final state covariance
% Final state
sigma_xf = [1e-2; ... % r_x
            1e-2; ... % r_y
            1e-3; ... % r_z
            1e-3; ... % v_x
            1e-3; ... % v_y
            1e-3; ... % v_z
            2e-2; ... % theta_1
            2e-2; ... % theta_2
            2e-2; ... % theta_3
            1e-2; ... % w_x
            1e-2; ... % w_y
            1e-2; ... % w_z
            1e-4]; % mass
Pf = diag(sigma_xf .^ 2);

% PTR algorithm parameters are defined in Deterministic


%% Get Dynamics
f = @(t, x, u, p) SymDynamicsEuler6DoF_convex_noumag(x, u, vehicle.L, vehicle.I, vehicle.alpha);

%% Measurement Model
% Measurement model (identity with noise)
g_0_stds = [10e-5; ... % r_x
            10e-5; ... % r_y
            50e-5; ... % r_z
            5e-5; ... % v_x
            5e-5; ... % v_y
            5e-5; ... % v_z
            1e-4; ... % theta_1
            1e-4; ... % theta_2
            1e-4; ... % theta_3
            5e-5; ... % w_x
            5e-5; ... % w_y
            5e-5; ... % w_z
            1e-5]; ... % mass

f_0 = @(t, x, u, p) x;
g_0 = @(t, x, u, p) diag(g_0_stds);

%% Specify Constraints
% Convex state path constraints
glideslope_angle_max = deg2rad(60);
h_glideslope = calculate_glideslope_offset(sigma_xf(1:nr) * norminv(1 - 1e-3 / 2), glideslope_angle_max);
glideslope_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) (norm(x(1:nr) + [zeros([nr - 1, 1]); h_glideslope]) + 0*sigma_mag_confidence(1e-3, nr) * norm(X_k_ref(1:nr, (tri(k - 1, nx) + 1):tri(k, nx)))) * cos(glideslope_angle_max) - (x(nr) + h_glideslope - norminv(1 - 1e-3) * norm(X_k_ref(nr, (tri(k - 1, nx) + 1):tri(k, nx))));
state_nonconvex_constraints = {glideslope_constraint};

state_convex_constraints = {};
control_convex_constraints ={};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex control constraints
max_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) norm(u(1:nr)) + sigma_mag_confidence(1e-3 / 2, nr) * norm(S_k(1:nr, :)) - thrust_magnitude_bound(S_k_ref, u_ref, k, t_k, T_max, m_0, alpha, nr, nx);
min_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) T_min / m_0 - (norm(u_ref(1:nr, k)) + (u_ref(1:nr, k)' / norm(u_ref(1:nr, k))) * (u(1:nr) - u_ref(1:nr, k)) - sigma_mag_confidence(1e-3 / 2, nr) * norm(S_k(1:nr, :)));
control_nonconvex_constraints = {max_thrust_constraint};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];
%% Specify Objective
stochastic_min_fuel_objective = @(x, u, p, X_k, S_k) einsum(@(k) norm(u(1:nr, k), 2) + sigma_mag_confidence(1e-2, nr) * norm(S_k(1:3, (tri(k - 1, nx) + 1):tri(k, nx)), 2), 1:Nu) * (t_k(2) - t_k(1));

%% Construct Problem Object
mod_ptr_sol = ptr_sol;
mod_ptr_sol.u = ptr_sol.u([1:nr, 5], :, :);
ptr_ops.iter_min = 4;
ptr_ops.Delta_min = 1e-5;
ptr_ops.w_vc = 1e6;

prob_6DoF.cont.f = f;
prob_6DoF.xf(nr) = sigma_xf(nr) * norminv(1 - 1e-3);
stoch_terminal_bc = @(x, p) [x(1:(nx - 1)) - prob_6DoF.xf; 0];
stoch_prob_6DoF = StochasticProblem.stochastify_discrete_problem(prob_6DoF, G, f_0, g_0, Phat0, Ptilde0, Pf, sol = mod_ptr_sol, objective = stochastic_min_fuel_objective, f = f, convex_constraints = convex_constraints, nonconvex_constraints = nonconvex_constraints, terminal_bc = stoch_terminal_bc, delta_t = delta_t);
[stoch_prob_6DoF, Delta] = stoch_prob_6DoF.discretize(stoch_prob_6DoF.guess.x, stoch_prob_6DoF.guess.u, stoch_prob_6DoF.guess.p);

norm(Delta)

%% Test Discretization on Initial Guess

%guess.u(1, :) = T_min * 1.5;
%guess.u(2, :) = T_min / 2000;
%guess.u(3, :) = vecnorm(guess.u(1:2, :));

% [prob_3DoF, Delta_disc] = prob_3DoF.discretize(prob_3DoF.guess.x, prob_3DoF.guess.u, prob_3DoF.guess.p);
% 
% x_disc = prob_3DoF.disc_prop(guess.x, guess.u, guess.p, K_k);
% 
% [t_cont, x_cont, u_cont] = prob_3DoF.cont_prop(guess.x, guess.u, guess.p, K_k);

% figure
% comparison_plot_3DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")
% 
% figure
% comparison_plot_3DoF_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess.u, u_cont, guess.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%% Solve Problem with PTR
stoch_ptr_sol = Stochastic_ptr(stoch_prob_6DoF, ptr_ops);

%% MC Simulations with Optimized Feedback Gain
K_k_opt = recover_gain_matrix(stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i));

K_k_ck = zeros([stoch_prob_6DoF.N, 1]);
for km = 1:(stoch_prob_6DoF.N - 1)
    K_k_ck(km) = norm(K_k_opt(:, :, km) * stoch_ptr_sol.X(:, (tri(km - 1, nx) + 1):tri(km, nx), stoch_ptr_sol.converged_i) - stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), stoch_ptr_sol.converged_i));
end

%%
for km = 1:(stoch_prob_6DoF.N - 1)
    X_norms(km) = norm(stoch_ptr_sol.X(3, (tri(km - 1, nx) + 1):tri(km, nx), stoch_ptr_sol.converged_i))
end

plot(X_norms)


%% Check Convergence for Gamma_k
Gamma_k = zeros([stoch_prob_6DoF.N, stoch_ptr_sol.converged_i]);
thrust_mag_k = zeros([stoch_prob_6DoF.N, stoch_ptr_sol.converged_i]);
thrust_min_k = zeros([stoch_prob_6DoF.N, stoch_ptr_sol.converged_i]);
thrust_mag_nom_k = zeros([stoch_prob_6DoF.N, stoch_ptr_sol.converged_i]);
max_thrust_constraint_evals = zeros([stoch_prob_6DoF.N, stoch_ptr_sol.converged_i]);
min_thrust_constraint_evals = zeros([stoch_prob_6DoF.N, stoch_ptr_sol.converged_i]);
for ms = 1:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_6DoF.N - 1)
        Gamma_k(km, ms) = thrust_magnitude_bound(stoch_ptr_sol.S(:, :, ms), squeeze(stoch_ptr_sol.u(:, :, ms)), km, t_k, T_max, m_0, alpha, nr, nx);
        thrust_mag_k(km, ms) = norm(stoch_ptr_sol.u(1:nr, km, ms)) + sigma_mag_confidence(1e-3 / 2, nr) * norm(stoch_ptr_sol.S(1:nr, (tri(km - 1, nx) + 1):tri(km, nx), ms));
        thrust_min_k(km, ms) = norm(stoch_ptr_sol.u(1:nr, km, ms)) - sigma_mag_confidence(1e-3 / 2, nr) * norm(stoch_ptr_sol.S(1:nr, (tri(km - 1, nx) + 1):tri(km, nx), ms));
        thrust_mag_nom_k(km, ms) = norm(stoch_ptr_sol.u(1:nr, km, ms));
    end
    ms
end

for ms = 2:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_6DoF.N - 1)
        max_thrust_constraint_evals(km, ms) = max_thrust_constraint(stoch_ptr_sol.x(1:nr, km, ms), stoch_ptr_sol.u(1:nr, km, ms), 0, 0, stoch_ptr_sol.S(1:nr, (tri(km - 1, nx) + 1):tri(km, nx), ms), stoch_ptr_sol.x(1:nr, km, ms - 1), stoch_ptr_sol.u(1:nr, :, ms - 1), 0, 0, stoch_ptr_sol.S(:, :, ms - 1), km);
    end
end

for ms = 2:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_6DoF.N - 1)
        min_thrust_constraint_evals(km, ms) = T_min - (vecnorm(stoch_ptr_sol.u(1:nr, km, ms)) * exp(stoch_ptr_sol.x(13, km, ms)) - sigma_mag_confidence(1e-3 / 2, nr) * norm(stoch_ptr_sol.S(1:nr, (tri(km - 1, nx) + 1):tri(km, nx), ms)));
    end
end
%%
figure
tiledlayout(2, 2)
nexttile
stairs(t_k(1:(end - 1)), Gamma_k(1:(end - 1), :) .* squeeze(exp(stoch_ptr_sol.x(7, 1:(end - 1), 1:stoch_ptr_sol.converged_i))))
title("\Gamma_k vs Time for All Iterations")
xlabel("Time [s]")
ylabel("[km / s2]")
legend("Iter " + string(1:stoch_ptr_sol.converged_i), Location="southeast")
grid on

nexttile
for ms = 1:(stoch_ptr_sol.converged_i - 1)
    stairs(t_k, abs(Gamma_k(:, ms + 1) - Gamma_k(:, ms))); hold on
end
hold off
title("\Delta\Gamma_k vs Time for All Iterations")
yscale("log")
xlabel("Time [s]")
ylabel("[km / s2]")
legend("Iter " + string(2:stoch_ptr_sol.converged_i) + " minus " + string(1:(stoch_ptr_sol.converged_i - 1)), Location="southeast")
grid on

nexttile
stairs(t_k(1:(end - 1)), thrust_mag_k(1:(end - 1),:) .* squeeze(exp(stoch_ptr_sol.x(13, 1:(end - 1), 1:stoch_ptr_sol.converged_i))));
title("||u_k|| + 99.9% Uncertainty Bound vs Time for All Iterations")
xlabel("Time [s]")
ylabel("[km / s2]")
legend("Iter " + string(1:stoch_ptr_sol.converged_i), Location="southeast")
grid on

nexttile
stairs(t_k(1:(end - 1)), thrust_mag_nom_k(1:(end - 1),:) .* squeeze(exp(stoch_ptr_sol.x(13, 1:(end - 1), 1:stoch_ptr_sol.converged_i))));
title("||u_k|| vs Time for All Iterations")
xlabel("Time [s]")
ylabel("[km / s2]")
legend("Iter " + string(1:stoch_ptr_sol.converged_i), Location="southeast")
grid on


sgtitle("\Gamma_k Convergence Plots")
%%
stoch_prob_6DoF.tolerances.AbsTol = 1e-8;
stoch_prob_6DoF.tolerances.RelTol = 1e-8;

m = 10;

t_ofb = zeros([stoch_prob_6DoF.N, m]);
x_ofb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
xhat_ofb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
Phat_ofb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
u_ofb = zeros([stoch_prob_6DoF.n.u, stoch_prob_6DoF.Nu, m]);

for i = 1:m
    [t_ofb(:, i), x_ofb(:, :, i), xhat_ofb(:, :, i), Phat_ofb(:, :, :, i), u_ofb(:, :, i)] = stoch_prob_6DoF.disc_prop(stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.p(:, stoch_ptr_sol.converged_i), K_k_opt);
    i
end

%%
defect = calculate_defect(stoch_prob_6DoF, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i),  stoch_ptr_sol.p(:, stoch_ptr_sol.converged_i))

%plot_6DoFc_trajectory(t_cont_sol, x_cont_sol, u_cont_sol, glideslope_angle_max, gimbal_max, T_min, T_max, step=1)

%% MC Simulations with Non-Optimized Feedback Gain

K_k = zeros([nu, nx, stoch_prob_6DoF.Nu]);
K_k(1, 1, :) = -1e-2;
K_k(2, 2, :) = -1e-2;

m = 1;

t_fb = zeros([stoch_prob_6DoF.N, m]);
x_fb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
xhat_fb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
Phat_fb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
u_fb = zeros([stoch_prob_6DoF.n.u, stoch_prob_6DoF.Nu, m]);

for i = 1:m
    [t_fb(:, i), x_fb(:, :, i), xhat_fb(:, :, i), Phat_fb(:, :, :, i), u_fb(:, :, i)] = stoch_prob_6DoF.disc_prop(stoch_prob_6DoF.guess.x, stoch_prob_6DoF.guess.u, stoch_prob_6DoF.guess.p, K_k);
    i
end

%% MC Simulations with No Feedback Control
t_no_fb = zeros([stoch_prob_6DoF.N, m]);
x_no_fb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
xhat_no_fb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
Phat_no_fb = zeros([stoch_prob_6DoF.n.x, stoch_prob_6DoF.n.x, stoch_prob_6DoF.N, m]);
u_no_fb = zeros([stoch_prob_6DoF.n.u, stoch_prob_6DoF.Nu, m]);

for i = 1:m
    [t_no_fb(:, i), x_no_fb(:, :, i), xhat_no_fb(:, :, i), Phat_no_fb(:, :, :, i), u_no_fb(:, :, i)] = stoch_prob_6DoF.disc_prop(stoch_prob_6DoF.guess.x, stoch_prob_6DoF.guess.u, stoch_prob_6DoF.guess.p, K_k * 0);
    i
end

%%
[Phat_k_opt, Pu_k_opt] = recover_est_covariances(stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i));

P_k_opt = Phat_k_opt + stoch_prob_6DoF.disc.Ptilde_k;
%%

plot_6DoF_MC_trajectories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, x_ofb, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, P_k_opt, t_k, xhat_no_fb, Pf, glideslope_angle_max, h_glideslope);

%%
plot_6DoF_MC_time_histories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i), t_k, x_ofb, u_ofb, t_k, stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i), t_k, xhat_no_fb, T_max, T_min, gimbal_max, true)

%%
plot_6DoF_MC_trajectories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, x_fb, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, P_k_opt, t_k, xhat_no_fb, Pf, glideslope_angle_max, h_glideslope)

%%
figure
covariance_plot(stoch_ptr_sol.x(1:9, end, stoch_ptr_sol.converged_i), squeeze(x_ofb(:, end, :)), squeeze(P_k_opt(1:9, 1:9, end)), Phat0(1:9, 1:9, :) + Ptilde0(1:9, 1:9, :), ["x [km]", "y [km]", "z", "v_x [km / s]", "v_y [km / s]", "vz", "\theta [rad]", "\theta [rad]", "\theta [rad]"], "State Dispersion at Final Node")

%% 3D plots/ animations???