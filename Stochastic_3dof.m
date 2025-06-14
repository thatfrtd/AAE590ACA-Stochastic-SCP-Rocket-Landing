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
Deterministic_3DoF_with_mass_convexified
%load("Deterministic_convexified_3DoF_workspace.mat");
%% Stochastic Optimization Parameters
nu = 2;
nr = 2;
nx = 7;
n_sigma_99 = sigma_mag_confidence(1e-2, nu);
n_sigma_99p9 = sigma_mag_confidence(1e-3, nu);

n_probit_99p9 = norminv(1 - 1e-3);

%Initial Covariance values
%% Define Stochastic Elements
% Initial estimated state
sigma_xhat0 = [30e-3; ... % r_x
            10e-3; ... % r_y
            6e-3; ... % v_x
            3e-3; ... % v_y
            deg2rad(1); ... % theta
            deg2rad(0.5); ... % w
            1e-4]; % mass
Phat0 = diag(sigma_xhat0 .^ 2);
Phat0(1:2, 1:2) = make_R2(-deg2rad(60)) * Phat0(1:2, 1:2) * make_R2(-deg2rad(60))';
Phat0(3:4, 3:4) = make_R2(-deg2rad(60)) * Phat0(3:4, 3:4) * make_R2(-deg2rad(60))';

%%

% Disturbance
sigma_accelx = 0.5e-3;
sigma_accely = 0.2e-3;
sigma_theta = 0;
sigma_ang_vel = deg2rad(0.05);
sigma_m = 1e-7;

delta_t = 1e-1;
G = @(t, x, u, p) sqrt(delta_t) * [zeros([2, 3]); ... % velocity
                   [sigma_accelx; sigma_accely] .* eye([2, 3]); ... % acceleration
                   zeros([1, 3]); ... % angular velocity
                   sigma_theta * [1, 1, 0];
                   [zeros([1, 2]), sigma_m]]; ... % mass flow 

% Initial state estimation error
sigma_xtilde0 = [3e-4; ... % r_x
            3e-4; ... % r_y
            1e-4; ... % v_x
            1e-4; ... % v_y
            deg2rad(0.05); ... % theta
            deg2rad(0.05); ... % w
            1e-5]; % mass
Ptilde0 = diag(sigma_xtilde0 .^ 2);

% Final state covariance
% Final state
sigma_xf = [1e-3; ... % r_x
            0.3e-3; ... % r_y
            1e-3; ... % v_x
            1e-3; ... % v_y
            deg2rad(1); ... % theta
            deg2rad(0.5); ... % w
            1e-4]; % mass
Pf = diag(sigma_xf .^ 2);

% PTR algorithm parameters are defined in Deterministic


%% Get Dynamics
f = @(t, x, u, p) SymDynamics3DoF_mass_convexified_noumag(t, x, u, vehicle.L, vehicle.I(2), vehicle.alpha);

%% Measurement Model
% Measurement model (identity with noise)
g_0_stds = [3e-4; ... % r_x
            3e-4; ... % r_y
            1e-4; ... % v_x
            1e-4; ... % v_y
            deg2rad(0.05); ... % theta
            deg2rad(0.05); ... % w
            1e-5]; ... % mass

f_0 = @(t, x, u, p) x;
g_0 = @(t, x, u, p) diag(g_0_stds);

% Gain 
K_k = -3*repmat([0, 1, 0, 0, 0, 0, 0; 0, 0, 0, 0, -1, 0, 0], 1, 1, prob_3DoF.Nu);

%% Specify Constraints
% Convex state path constraints
glideslope_angle_max = deg2rad(50);
h_glideslope = calculate_glideslope_offset(sigma_xf(1:2) * norminv(1 - 1e-3 / 2), glideslope_angle_max);
glideslope_constraint = @(x, u, p, X_k, S_k) norm(x(1)) + sigma_mag_confidence(1e-3, nr - 1) * (norm(X_k(1, :))) - tan(glideslope_angle_max) * (x(2) + h_glideslope - norminv(1 - 1e-3) * norm(X_k(2, :)));

% sigma_max = [200e-3; ... % r_x
%             200000e-3; ... % r_y
%             20000e-3; ... % v_x
%             200000e-3; ... % v_y
%             deg2rad(25000); ... % theta
%             deg2rad(10000)]; ... % w
% max_uncertainty_constraint = @(x, u, p, X_k, S_k) norminv(1 - 1e-3) * norms(X_k(1:6, :), 2, 2) - sigma_max;
state_convex_constraints = {glideslope_constraint};

gimbal_angle_max = deg2rad(6);
gimbal_constraint = @(x, u, p, X_k, S_k) norm(u(2)) + sigma_mag_confidence(1e-3, nr - 1) * (norm(S_k(2, :))) - tan(gimbal_angle_max) * (u(1) - norminv(1-1e-3) * (norm(S_k(1, :))));
control_convex_constraints ={gimbal_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state path constraints
state_nonconvex_constraints = {};

% Nonconvex control constraints
max_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) norm(u(1:2)) + sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k) - thrust_magnitude_bound(S_k_ref, u_ref, k, t_k, T_max, m_0, alpha, nu, nx);
%min_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) T_min / m_0 - (norm(u_ref(:, k)) + (u_ref(:, k)' / norm(u_ref(:, k))) * (u - u_ref(:, k)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k));
min_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) T_min * exp(-x_ref(7)) - (norm(u_ref(:, k)) + (u_ref(:, k)' / norm(u_ref(:, k))) * (u - u_ref(:, k)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k));
control_nonconvex_constraints = {max_thrust_constraint, min_thrust_constraint};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];
%% Specify Objective
stochastic_min_fuel_objective = @(x, u, p, X_k, S_k) einsum(@(k) norm(u(1:2, k), 2) + sigma_mag_confidence(1e-2, nu) * norm(S_k(:, (tri(k - 1, nx) + 1):tri(k, nx)), 2), 1:Nu) * (t_k(2) - t_k(1));

%% Construct Problem Object
mod_ptr_sol = ptr_sol;
mod_ptr_sol.x = ptr_sol.x(:, :, ptr_sol.converged_i);
mod_ptr_sol.u = ptr_sol.u(1:2, :, ptr_sol.converged_i);
mod_ptr_sol.p = ptr_sol.p(:, ptr_sol.converged_i);
ptr_ops.iter_min = 4;
ptr_ops.w_vc = 1e5;
ptr_ops.w_tr(:) = 1e-1; 
ptr_ops.delta_tol = 2e-3;

prob_3DoF.cont.f = f;
prob_3DoF.xf = [0; sigma_xf(2) * norminv(1 - 1e-3); 0; 0; pi/2; 0];
stoch_terminal_bc = @(x, p) [x(1:6) - prob_3DoF.xf; 0];
stoch_prob_3DoF = StochasticProblem.stochastify_discrete_problem(prob_3DoF, G, f_0, g_0, Phat0, Ptilde0, Pf, sol = mod_ptr_sol, objective = stochastic_min_fuel_objective, f = f, convex_constraints = convex_constraints, nonconvex_constraints = nonconvex_constraints, terminal_bc = stoch_terminal_bc);
[stoch_prob_3DoF, Delta] = stoch_prob_3DoF.discretize(stoch_prob_3DoF.guess.x, stoch_prob_3DoF.guess.u, stoch_prob_3DoF.guess.p);

%stoch_prob_3DoF = StochasticProblem(x_0, x_f, Phat0, Ptilde0, Pf, N, u_hold, tf, f, G, f_0, g_0, mod_ptr_sol, convex_constraints, stochastic_min_fuel_objective, scale = scale, nonconvex_constraints = nonconvex_constraints,  terminal_bc = stoch_terminal_bc);

norm(Delta)
%% Give Python Solver Info
stoch_prob_3DoF.params = [T_min, T_max, alpha, glideslope_angle_max, gimbal_max, h_glideslope, m_0];

%% Solve Problem with PTR
%stoch_prob_3DoF.guess.x = stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i);
%stoch_prob_3DoF.guess.u = stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i);
%X_k_guess = stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i);
%S_k_guess = stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i);
stoch_ptr_sol = Stochastic_ptr(stoch_prob_3DoF, ptr_ops, "CVXPY");%, X_k_guess = X_k_guess, S_k_guess = S_k_guess);

%%
%stoch_ptr_sol.converged_i = 2;

%% MC Simulations with Optimized Feedback Gain
K_k_opt = recover_gain_matrix(stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i));

K_k_ck = zeros([stoch_prob_3DoF.N, 1]);
for km = 1:(stoch_prob_3DoF.N - 1)
    K_k_ck(km) = norm(K_k_opt(:, :, km) * stoch_ptr_sol.X(:, (tri(km - 1, nx) + 1):tri(km, nx), stoch_ptr_sol.converged_i) - stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), stoch_ptr_sol.converged_i));
end


%% Check Convergence for Gamma_k
Gamma_k = zeros([stoch_prob_3DoF.N, stoch_ptr_sol.converged_i]);
thrust_mag_k = zeros([stoch_prob_3DoF.N, stoch_ptr_sol.converged_i]);
thrust_min_k = zeros([stoch_prob_3DoF.N, stoch_ptr_sol.converged_i]);
thrust_mag_nom_k = zeros([stoch_prob_3DoF.N, stoch_ptr_sol.converged_i]);
max_thrust_constraint_evals = zeros([stoch_prob_3DoF.N, stoch_ptr_sol.converged_i]);
min_thrust_constraint_evals = zeros([stoch_prob_3DoF.N, stoch_ptr_sol.converged_i]);
for ms = 1:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_3DoF.N - 1)
        Gamma_k(km, ms) = thrust_magnitude_bound(stoch_ptr_sol.S(:, :, ms), squeeze(stoch_ptr_sol.u(:, :, ms)), km, t_k, T_max, m_0, alpha, nu, nx);
        thrust_mag_k(km, ms) = norm(stoch_ptr_sol.u(1:2, km, ms)) + sigma_mag_confidence(1e-3 / 2, nu) * norm(stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), ms));
        thrust_min_k(km, ms) = norm(stoch_ptr_sol.u(1:2, km, ms)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), ms));
        thrust_mag_nom_k(km, ms) = norm(stoch_ptr_sol.u(1:2, km, ms));
    end
    ms
end

for ms = 2:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_3DoF.N - 1)
        max_thrust_constraint_evals(km, ms) = max_thrust_constraint(stoch_ptr_sol.x(1:2, km, ms), stoch_ptr_sol.u(1:2, km, ms), 0, 0, stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), ms), stoch_ptr_sol.x(1:2, km, ms - 1), stoch_ptr_sol.u(1:2, :, ms - 1), 0, 0, stoch_ptr_sol.S(:, :, ms - 1), km);
    end
end

for ms = 2:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_3DoF.N - 1)
        min_thrust_constraint_evals(km, ms) = T_min - (vecnorm(stoch_ptr_sol.u(1:2, km, ms)) * exp(stoch_ptr_sol.x(7, km, ms)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), ms)));
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
stairs(t_k(1:(end - 1)), thrust_min_k(1:(end - 1),:) .* squeeze(exp(stoch_ptr_sol.x(7, 1:(end - 1), 1:stoch_ptr_sol.converged_i))));
title("||u_k|| - 99.9% Uncertainty Bound vs Time for All Iterations")
xlabel("Time [s]")
ylabel("[km / s2]")
legend("Iter " + string(1:stoch_ptr_sol.converged_i), Location="southeast")
grid on

nexttile
stairs(t_k(1:(end - 1)), thrust_mag_nom_k(1:(end - 1),:) .* squeeze(exp(stoch_ptr_sol.x(7, 1:(end - 1), 1:stoch_ptr_sol.converged_i))));
title("||u_k|| vs Time for All Iterations")
xlabel("Time [s]")
ylabel("[km / s2]")
legend("Iter " + string(1:stoch_ptr_sol.converged_i), Location="southeast")
grid on


sgtitle("\Gamma_k Convergence Plots")
%%
m = 100;

t_ofb = zeros([stoch_prob_3DoF.N, m]);
x_ofb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
xhat_ofb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
Phat_ofb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
u_ofb = zeros([stoch_prob_3DoF.n.u, stoch_prob_3DoF.Nu, m]);

parfor i = 1:m
    [t_ofb(:, i), x_ofb(:, :, i), xhat_ofb(:, :, i), Phat_ofb(:, :, :, i), u_ofb(:, :, i)] = stoch_prob_3DoF.disc_prop(stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.p(:, stoch_ptr_sol.converged_i), K_k_opt);
    i
end

%% MC Simulations with Non-Optimized Feedback Gain
K_k = -1e-2*repmat([1, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 0, 0, 0], 1, 1, stoch_prob_3DoF.Nu);

m = 100;

t_fb = zeros([stoch_prob_3DoF.N, m]);
x_fb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
xhat_fb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
Phat_fb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
u_fb = zeros([stoch_prob_3DoF.n.u, stoch_prob_3DoF.Nu, m]);

parfor i = 1:m
    [t_fb(:, i), x_fb(:, :, i), xhat_fb(:, :, i), Phat_fb(:, :, :, i), u_fb(:, :, i)] = stoch_prob_3DoF.disc_prop(stoch_prob_3DoF.guess.x, stoch_prob_3DoF.guess.u, stoch_prob_3DoF.guess.p, K_k);
    i
end

%% MC Simulations with No Feedback Control
t_no_fb = zeros([stoch_prob_3DoF.N, m]);
x_no_fb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
xhat_no_fb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
Phat_no_fb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.n.x, stoch_prob_3DoF.N, m]);
u_no_fb = zeros([stoch_prob_3DoF.n.u, stoch_prob_3DoF.Nu, m]);

parfor i = 1:m
    [t_no_fb(:, i), x_no_fb(:, :, i), xhat_no_fb(:, :, i), Phat_no_fb(:, :, :, i), u_no_fb(:, :, i)] = stoch_prob_3DoF.disc_prop(stoch_prob_3DoF.guess.x, stoch_prob_3DoF.guess.u, stoch_prob_3DoF.guess.p, K_k * 0);
    i
end

%%
[Phat_k_opt, Pu_k_opt] = recover_est_covariances(stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i));

P_k_opt = Phat_k_opt + stoch_prob_3DoF.disc.Ptilde_k;
%%

plot_3DoF_MC_trajectories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, x_ofb, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, P_k_opt, t_k, xhat_no_fb, Pf, glideslope_angle_max, h_glideslope, x_ref_solution = ptr_sol.x(:, :, ptr_sol.converged_i));

%%
plot_3DoF_MC_time_histories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i), t_k, x_ofb, u_ofb, t_k, stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i), t_k, xhat_no_fb, glideslope_angle_max, h_glideslope, T_max, T_min, gimbal_angle_max, true)

%%
plot_3DoF_MC_trajectories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, x_fb, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, P_k_opt, t_k, xhat_no_fb, Pf, glideslope_angle_max, h_glideslope)

%%
figure
covariance_plot(stoch_ptr_sol.x(:, end, stoch_ptr_sol.converged_i), squeeze(x_ofb(:, end, :)), squeeze(P_k_opt(:, :, end)), Pf, ["x [km]", "y [km]", "v_x [km / s]", "v_y [km / s]", "\theta [rad]", "\omega [rad / s]", "m [kg]"], "State Dispersion at Final Node")
%%
figure
mc_dv = squeeze(sum(vecnorm(u_ofb, 2, 1) * (t_k(2) - t_k(1)), 2));
mc_max = max(mc_dv);

ptr_dv = sum(vecnorm(stoch_ptr_sol.u(:,:,stoch_ptr_sol.converged_i), 2, 1) * (t_k(2) - t_k(1)));

ptr_bound = stoch_ptr_sol.info.J;

histogram(mc_dv); hold on
xline(ptr_dv, "--", "\Delta V nominal", LabelHorizontalAlignment="right",LineWidth=1); hold on
xline(ptr_bound, "--", "\Delta V_{99} (bound)", LabelHorizontalAlignment="right",LineWidth=1); hold on
xline(mc_max, "--", "\Delta V_{99} (MC)", LabelHorizontalAlignment="right",LineWidth=1);
title("Delta V Distribution")
xlabel("Delta V [km / s]")
grid on