%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter, Atharva Awasthi
% Created On: 19 April, 2025
% Description: Stochastic 2DoF (all translational) landing of rocket using 
% PTR SCP algorithm
% Most Recent Change: 27 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HOW DOES FEEDBACK WITH THE THRUST MAGNITUDE PART OF THE CONTROL VECTOR
% WORK?? DO THE DYNAMICS JUST NOT HAVE IT IN THEM?
% - changed it so the control is just [u1; u2] for stochastic case

%% Stochastic Optimization Parameters
nu = 2;
nr = 2;
nx = 5;

%% Get deterministic solution for an initial guess 
Deterministic_2DoF_linear

t_k = linspace(0, prob_2DoF.tf, prob_2DoF.N);

%% Define Stochastic Elements
% Initial estimated state
sigma_xhat0 = [30e-3; ... % r_x
            10e-3; ... % r_y
            6e-3; ... % v_x
            3e-3; ... % v_y
            1e-4]; % mass
Phat0 = diag(sigma_xhat0 .^ 2);
Phat0(1:2, 1:2) = make_R2(-deg2rad(60)) * Phat0(1:2, 1:2) * make_R2(-deg2rad(60))';
Phat0(3:4, 3:4) = make_R2(-deg2rad(60)) * Phat0(3:4, 3:4) * make_R2(-deg2rad(60))';


% Initial state estimation error
sigma_xtilde0 = [3e-4; ... % r_x
            3e-4; ... % r_y
            1e-4; ... % v_x
            1e-4; ... % v_y
            1e-5]; % mass
Ptilde0 = diag(sigma_xtilde0 .^ 2);

% Final state
sigma_xf = [1e-3; ... % r_x
            0.3e-3; ... % r_y
            1e-3; ... % v_x
            1e-3; ... % v_y
            1e-4]; % mass
Pf = diag(sigma_xf .^ 2);

% Disturbance
sigma_accelx = 0.5e-3;
sigma_accely = 0.2e-3;
sigma_m = 1e-7;

delta_t = 1e-1;
G = @(t, x, u, p) sqrt(delta_t) * [zeros([2, 3]); ... % velocity
                   [sigma_accelx; sigma_accely] .* eye([2, 3]); ... % acceleration
                   [zeros([1, 2]), sigma_m]]; ... % mass flow 
% Brownian noise approximation
w = @(n) randn([3, n]);

% Noise timespan
tspan = 0:delta_t:prob_2DoF.tf;

% Integration tolerances
tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

% Measurement functions
g_0_stds = [3e-4; ... % r_x
            3e-4; ... % r_y
            1e-4; ... % v_x
            1e-4; ... % v_y
            1e-5]; ... % mass

f_0 = @(t, x, u, p) x; % full state measurement
g_0 = @(t, x, u, p) diag(g_0_stds);

%% Stochastify Objective and Constraints
% Objective
stochastic_min_fuel_objective = @(x, u, p, X_k, S_k) einsum(@(k) norm(u(1:2, k), 2) + sigma_mag_confidence(1e-2, nu) * norm(S_k(:, (tri(k - 1, nx) + 1):tri(k, nx)), 2), 1:Nu) * (t_k(2) - t_k(1));

% Convex state path constraints
glideslope_angle_max = deg2rad(50);
h_glideslope = calculate_glideslope_offset(sigma_xf(1:2) * norminv(1 - 1e-3 / 2), glideslope_angle_max);
glideslope_constraint = @(x, u, p, X_k, S_k) norm(x(1)) + sigma_mag_confidence(1e-3, nr - 1) * (norm(X_k(1, :))) - tan(glideslope_angle_max) * (x(2) + h_glideslope - norminv(1 - 1e-3) * norm(X_k(2, :)));
state_convex_constraints = {glideslope_constraint};

% Nonconvex state path constraints
state_nonconvex_constraints = {};

% Convex control constraints
%lcvx_constraint = @(x, u, p, X_k, S_k, u_mag) norm(u(1:2)) - u_mag;
control_convex_constraints = {};%{lcvx_constraint, min_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex control constraints
max_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) norm(u) + sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k) - thrust_magnitude_bound(S_k_ref, u_ref, k, t_k, T_max, m_0, alpha, nu, nx);
min_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) T_min * exp(-x(5)) - (norm(u_ref(:, k)) + (u_ref(:, k)' / norm(u_ref(:, k))) * (u - u_ref(:, k)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k));

control_nonconvex_constraints = {max_thrust_constraint, min_thrust_constraint};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

%% Set Up StochasticProblem
f_stoch = @(t, x, u, p) SymDynamics2DoF_linear_noumag(t, x, u, vehicle.alpha);
ptr_sol_mod = ptr_sol;
ptr_sol_mod.u = ptr_sol_mod.u(1:2, :, :);
prob_2DoF.xf = [0; sigma_xf(2) * norminv(1 - 1e-3); 0; 0];
stoch_terminal_bc = @(x, p) [x(1:4) - prob_2DoF.xf; 0];
stoch_prob_2DoF = StochasticProblem.stochastify_discrete_problem(prob_2DoF, G, f_0, g_0, Phat0, Ptilde0, Pf, f = f_stoch, sol = ptr_sol_mod, objective = stochastic_min_fuel_objective, convex_constraints = convex_constraints, nonconvex_constraints = nonconvex_constraints, terminal_bc = stoch_terminal_bc);
[stoch_prob_2DoF, Delta] = stoch_prob_2DoF.discretize(stoch_prob_2DoF.guess.x, stoch_prob_2DoF.guess.u, stoch_prob_2DoF.guess.p);

%% Solve Stochastic Optimization Problem with PTR
ptr_ops.iter_max = 20;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 1e4;
ptr_ops.w_tr = ones(1, Nu) * 1e-1;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;
%%
tic
stoch_ptr_sol = Stochastic_ptr(stoch_prob_2DoF, ptr_ops, slack_control = false);
toc
stoch_ptr_sol.converged_i
%% MC Simulations with Optimized Feedback Gain
K_k_opt = recover_gain_matrix(stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i));

K_k_ck = zeros([stoch_prob_2DoF.N, 1]);
for km = 1:(stoch_prob_2DoF.N - 1)
    K_k_ck(km) = norm(K_k_opt(:, :, km) * stoch_ptr_sol.X(:, (tri(km - 1, nx) + 1):tri(km, nx), stoch_ptr_sol.converged_i) - stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), stoch_ptr_sol.converged_i));
end


%% Check Convergence for Gamma_k
Gamma_k = zeros([stoch_prob_2DoF.N, stoch_ptr_sol.converged_i]);
thrust_mag_k = zeros([stoch_prob_2DoF.N, stoch_ptr_sol.converged_i]);
thrust_min_k = zeros([stoch_prob_2DoF.N, stoch_ptr_sol.converged_i]);
thrust_mag_nom_k = zeros([stoch_prob_2DoF.N, stoch_ptr_sol.converged_i]);
max_thrust_constraint_evals = zeros([stoch_prob_2DoF.N, stoch_ptr_sol.converged_i]);
min_thrust_constraint_evals = zeros([stoch_prob_2DoF.N, stoch_ptr_sol.converged_i]);
for ms = 1:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_2DoF.N - 1)
        Gamma_k(km, ms) = thrust_magnitude_bound(stoch_ptr_sol.S(:, :, ms), squeeze(stoch_ptr_sol.u(:, :, ms)), km, t_k, T_max, m_0, alpha, nu, nx);
        thrust_mag_k(km, ms) = norm(stoch_ptr_sol.u(1:2, km, ms)) + sigma_mag_confidence(1e-3 / 2, nu) * norm(stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), ms));
        thrust_min_k(km, ms) = norm(stoch_ptr_sol.u(1:2, km, ms)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), ms));
        thrust_mag_nom_k(km, ms) = norm(stoch_ptr_sol.u(1:2, km, ms));
    end
    ms
end

for ms = 2:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_2DoF.N - 1)
        max_thrust_constraint_evals(km, ms) = max_thrust_constraint(stoch_ptr_sol.x(1:2, km, ms), stoch_ptr_sol.u(1:2, km, ms), 0, 0, stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), ms), stoch_ptr_sol.x(1:2, km, ms - 1), stoch_ptr_sol.u(1:2, :, ms - 1), 0, 0, stoch_ptr_sol.S(:, :, ms - 1), km);
    end
end

for ms = 2:stoch_ptr_sol.converged_i
    for km = 1:(stoch_prob_2DoF.N - 1)
        min_thrust_constraint_evals(km, ms) = T_min - (vecnorm(stoch_ptr_sol.u(1:2, km, ms)) * exp(stoch_ptr_sol.x(5, km, ms)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(stoch_ptr_sol.S(:, (tri(km - 1, nx) + 1):tri(km, nx), ms)));
    end
end
%%
figure
tiledlayout(2, 2)
nexttile
stairs(t_k(1:(end - 1)), Gamma_k(1:(end - 1), :) .* squeeze(exp(stoch_ptr_sol.x(5, 1:(end - 1), 1:stoch_ptr_sol.converged_i))))
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
stairs(t_k(1:(end - 1)), thrust_min_k(1:(end - 1),:) .* squeeze(exp(stoch_ptr_sol.x(5, 1:(end - 1), 1:stoch_ptr_sol.converged_i))));
title("||u_k|| - 99.9% Uncertainty Bound vs Time for All Iterations")
xlabel("Time [s]")
ylabel("[km / s2]")
legend("Iter " + string(1:stoch_ptr_sol.converged_i), Location="southeast")
grid on

nexttile
stairs(t_k(1:(end - 1)), thrust_mag_nom_k(1:(end - 1),:) .* squeeze(exp(stoch_ptr_sol.x(5, 1:(end - 1), 1:stoch_ptr_sol.converged_i))));
title("||u_k|| vs Time for All Iterations")
xlabel("Time [s]")
ylabel("[km / s2]")
legend("Iter " + string(1:stoch_ptr_sol.converged_i), Location="southeast")
grid on


sgtitle("\Gamma_k Convergence Plots")
%%
m = 100;

t_ofb = zeros([stoch_prob_2DoF.N, m]);
x_ofb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
xhat_ofb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
Phat_ofb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
u_ofb = zeros([stoch_prob_2DoF.n.u, stoch_prob_2DoF.Nu, m]);

parfor i = 1:m
    [t_ofb(:, i), x_ofb(:, :, i), xhat_ofb(:, :, i), Phat_ofb(:, :, :, i), u_ofb(:, :, i)] = stoch_prob_2DoF.disc_prop(stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.p(:, stoch_ptr_sol.converged_i), K_k_opt);
    i
end

%% MC Simulations with Non-Optimized Feedback Gain
K_k = -1e-2*repmat([1, 0, 0, 0, 0; 0, 1, 0, 0, 0], 1, 1, prob_2DoF.Nu);

m = 100;

t_fb = zeros([stoch_prob_2DoF.N, m]);
x_fb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
xhat_fb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
Phat_fb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
u_fb = zeros([stoch_prob_2DoF.n.u, stoch_prob_2DoF.Nu, m]);

parfor i = 1:m
    [t_fb(:, i), x_fb(:, :, i), xhat_fb(:, :, i), Phat_fb(:, :, :, i), u_fb(:, :, i)] = stoch_prob_2DoF.disc_prop(stoch_prob_2DoF.guess.x, stoch_prob_2DoF.guess.u, stoch_prob_2DoF.guess.p, K_k);
    i
end

%% MC Simulations with No Feedback Control
t_no_fb = zeros([stoch_prob_2DoF.N, m]);
x_no_fb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
xhat_no_fb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
Phat_no_fb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
u_no_fb = zeros([stoch_prob_2DoF.n.u, stoch_prob_2DoF.Nu, m]);

parfor i = 1:m
    [t_no_fb(:, i), x_no_fb(:, :, i), xhat_no_fb(:, :, i), Phat_no_fb(:, :, :, i), u_no_fb(:, :, i)] = stoch_prob_2DoF.disc_prop(stoch_prob_2DoF.guess.x, stoch_prob_2DoF.guess.u, stoch_prob_2DoF.guess.p, K_k * 0);
    i
end
%% Calculate X_k and S_k from P_k and K_k
% [A_k, B_k, E_k, c_k, G_k, Delta] = discretize_stochastic_dynamics_ZOH(stoch_prob_2DoF.cont.f, prob_2DoF.cont.A, prob_2DoF.cont.B, prob_2DoF.cont.E, prob_2DoF.cont.c, G, prob_2DoF.N, t_k, ptr_sol.x(:, :, ptr_sol.converged_i), stoch_ptr_sol.u(:, :, ptr_sol.converged_i), ptr_sol.p(:, ptr_sol.converged_i), prob_2DoF.tolerances);
% 
% Ptilde_k = zeros([prob_2DoF.n.x, prob_2DoF.n.x, numel(t_k)]);
% Ptilde_k(:, :, 1) = diag(sigma_xtilde0 .^ 2);
% for k = 1:(numel(t_k) - 1)
%     Ptilde_k(:, :, k + 1) = (A_k(:, :, k) + B_k(:, :, k) * K_k_opt(:, :, k)) * Ptilde_k(:, :, k) * (A_k(:, :, k) + B_k(:, :, k) * K_k_opt(:, :, k))' + G_k(:, :, k) * G_k(:, :, k)';
% end
% 
% X_k = zeros(size(Ptilde_k));
% for k = 1:numel(t_k)
%     X_k(:, :, k) = chol(Ptilde_k(:, :, k), "lower");
% end
% S_k = pagemtimes(K_k, X_k(:, :, 1:prob_2DoF.Nu));

%%
[Phat_k_opt, Pu_k_opt] = recover_est_covariances(stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i));

P_k_opt = Phat_k_opt + stoch_prob_2DoF.disc.Ptilde_k;
%%

plot_2DoF_MC_trajectories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, x_ofb, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, P_k_opt, t_k, xhat_no_fb, Pf, glideslope_angle_max, h_glideslope, x_ref_solution = ptr_sol.x(:, :, ptr_sol.converged_i), title = ""); hold on

%%
plot_2DoF_MC_time_histories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i), t_k, x_ofb, u_ofb, t_k, stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i), t_k, xhat_no_fb, T_max, T_min, true)

%%
%plot_2DoF_MC_trajectories(t_k, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, x_fb, stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), t_k, P_k_opt, t_k, xhat_no_fb, Pf, glideslope_angle_max, h_glideslope)

%%
figure
covariance_plot(stoch_ptr_sol.x(1:4, end, stoch_ptr_sol.converged_i), squeeze(x_ofb(1:4, end, :)), squeeze(P_k_opt(1:4, 1:4, end)), Pf, ["x [km]", "y [km]", "v_x [km / s]", "v_y [km / s]"], "")
