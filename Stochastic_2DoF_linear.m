%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 19 April, 2025
% Description: Stochastic 2DoF (all translational) landing of rocket using 
% PTR SCP algorithm
% Most Recent Change: 20 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% HOW DOES FEEDBACK WITH THE THRUST MAGNITUDE PART OF THE CONTROL VECTOR
% WORK?? DO THE DYNAMICS JUST NOT HAVE IT IN THEM?

%% Stochastic Optimization Parameters
nu = 2;
n_sigma_99 = sigma_mag_confidence(1e-2, nu);
n_sigma_99p9 = sigma_mag_confidence(1e-3, nu);

%% Get deterministic solution for an initial guess 
Deterministic_2DoF_linear

t_k = linspace(0, prob_2DoF.tf, prob_2DoF.N);
if prob_2DoF.u_hold == "ZOH"
    u_func = @(t, x) interp1(t_k(1:prob_2DoF.Nu), ptr_sol.u(:, :, ptr_sol.converged_i)', t, "previous", "extrap")';
elseif prob_2DoF.u_hold == "FOH"
    u_func = @(t, x) interp1(t_k(1:prob_2DoF.Nu), ptr_sol.u(:, :, ptr_sol.converged_i)', t)';
end

%% Define Stochastic Elements
% Initial state
sigma_x0 = [10e-3; ... % r_x
            50e-3; ... % r_y
            5e-3; ... % v_x
            5e-3; ... % v_y
            1e-4]; % mass
P0 = diag(sigma_x0 .^ 2);

% Final state
sigma_xf = [10e-3; ... % r_x
            50e-3; ... % r_y
            5e-3; ... % v_x
            5e-3; ... % v_y
            1e-4]; % mass
Pf = diag(sigma_xf .^ 2);

% Disturbance
sigma_accelx = 0.5e-4;
sigma_accely = 0.1e-4;
sigma_m = 1e-4;
G = @(t, x, u, p) [zeros([2, 3]); ... % velocity
                   [sigma_accelx; sigma_accely] .* eye([2, 3]); ... % acceleration
                   [zeros([1, 2]), sigma_m]]; ... % mass flow 
% Brownian noise approximation
w = @(n) randn([3, n]);
delta_t = 1e-0;

% Noise timespan
tspan = 0:delta_t:prob_2DoF.tf;

% Integration tolerances
tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

% Measurement functions
g_0_stds = [10e-3; ... % r_x
            50e-3; ... % r_y
            5e-3; ... % v_x
            5e-3; ... % v_y
            1e-2]; ... % mass

f_0 = @(t, x, u, p) x; % full state measurement
g_0 = @(t, x, u, p) diag(g_0_stds);

%% Stochastify Objective and Constraints
% Objective
stochastic_min_fuel_objective = @(x, u, p, X_k, S_k) einsum(@(k) norm(u(1:2, k), 2) + n_sigma_99 * norm(S_k(:, (k - 1) * prob_2DoF.n.x + (1:prob_2DoF.n.x)), 2), 1:Nu) * delta_t;

% Define bounds
z_lb = @(t) log(m_0 - alpha * T_max * t);
z_lb_k = z_lb(t_k);

% Convex state path constraints
glideslope_constraint = @(t, x, u, p) norm(x(1:2)) - x(2) / cos(pi/2 - gamma_min);
min_mass_constraint = @(t, x, u, p) z_lb(t) - x(5);
state_convex_constraints = {}; % Ignoring state constriants for now !!!!

% Convex control constraints
control_convex_constraints = {};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex control constraints
lcvx_thrust_constraint = @(t, x, u, p, x_ref, u_ref, p_ref, S_K_ref, k) norm(u(1:2)) + n_sigma_99p9 * norm(S_k) - T_max / m_0 * exp(einsum(@(i) alpha * (t_k(2) - t_k(1)) * max(0, norm(u_ref(1:2, i)) - sigma_mag_confidence(1e-3 / (2 * k), nu) * norm(S_k_ref)), 1:(k - 1)));
control_nonconvex_constraints = {lcvx_thrust_constraint};

% Combine nonconvex constraints
nonconvex_constraints = [control_nonconvex_constraints];

%% Set Up StochasticProblem
stoch_prob_2DoF = StochasticProblem.stochastify_discrete_problem(prob_2DoF, G, f_0, g_0, P0, Pf, sol = ptr_sol, objective = stochastic_min_fuel_objective, convex_constraints = convex_constraints, nonconvex_constraints = nonconvex_constraints);
[stoch_prob_2DoF, Delta] = stoch_prob_2DoF.discretize(stoch_prob_2DoF.guess.x, stoch_prob_2DoF.guess.u, stoch_prob_2DoF.guess.p);

%% Solve Stochastic Optimization Problem with PTR
stoch_ptr_sol = Stochastic_ptr(stoch_prob_2DoF, ptr_ops);

%% MC Simulations with Non-Optimized Feedback Gain
K_k = -1e-2*repmat([1, 0, 0, 0, 0; 0, 1, 0, 0, 0; 1, 1, 0, 0, 0], 1, 1, prob_2DoF.Nu);

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
[A_k, B_k, E_k, c_k, G_k, Delta] = discretize_stochastic_dynamics_ZOH(prob_2DoF.cont.f, prob_2DoF.cont.A, prob_2DoF.cont.B, prob_2DoF.cont.E, prob_2DoF.cont.c, G, prob_2DoF.N, t_k, ptr_sol.x(:, :, ptr_sol.converged_i), ptr_sol.u(:, :, ptr_sol.converged_i), ptr_sol.p(:, ptr_sol.converged_i), prob_2DoF.tolerances);

P_k = zeros([prob_2DoF.n.x, prob_2DoF.n.x, numel(t_k)]);
P_k(:, :, 1) = diag(sigma_x0 .^ 2);
for k = 1:(numel(t_k) - 1)
    P_k(:, :, k + 1) = (A_k(:, :, k) + B_k(:, :, k) * K_k(:, :, k)) * P_k(:, :, k) * (A_k(:, :, k) + B_k(:, :, k) * K_k(:, :, k))' + G_k(:, :, k) * G_k(:, :, k)';
end

X_k = zeros(size(P_k));
for k = 1:numel(t_k)
    X_k(:, :, k) = chol(P_k(:, :, k), "lower");
end
S_k = pagemtimes(K_k, X_k(:, :, 1:prob_2DoF.Nu));

%%
plot_2DoF_MC_trajectories(t_k, stoch_prob_2DoF.guess.x, t_k, xhat_fb, stoch_prob_2DoF.guess.x, t_k, X_k, t_k, xhat_no_fb, Pf, pi / 2 - gamma_min)
plot_2DoF_MC_time_histories(t_k, stoch_prob_2DoF.guess.x, stoch_prob_2DoF.guess.u, t_k, xhat_fb, u_fb, t_k, X_k, S_k, t_k, xhat_no_fb, T_max, T_min, true)

%%
