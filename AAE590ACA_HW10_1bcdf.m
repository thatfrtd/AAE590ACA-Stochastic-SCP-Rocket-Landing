%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 14 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get HW4c Solution
AAE590ACA_HW4_Q1c_solution

%%

nu = 2;
n_sigma_99 = sigma_mag_confidence(1e-2, nu);
n_sigma_99p9 = sigma_mag_confidence(1e-3, nu);

tri = @(k) (k + 1) * k / 2 * 4;

dt = 0.1;
tspan_control = 0:dt:tf;
N = numel(tspan_control);

%% Define Stochastic Elements
% Covariance
sigma_pos = 3e-4;
sigma_vel = 3e-4;

Phat0 = diag([sigma_pos, sigma_pos, sigma_vel, sigma_vel] .^ 2);
Ptilde0 = zeros([4,4]);
Pf = diag([3e-3, 3e-3, 3e-3, 3e-3] .^ 2);

% Disturbance
sigma_dist = 3e-2;
G = @(t, x, u, p)  [zeros([2, 2]); ... % velocity
                   sigma_dist * eye(2)]; ... % mass flow 
% Brownian noise approximation
delta_t = 1e-3;
w = @(n) randn([3, n]);

% Noise timespan
tspan = 0:delta_t:tf;

% Integration tolerances
tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

%% Stochastify Objective and Constraints
% Objective
stochastic_min_fuel_objective = @(x, u, p, X_k, S_k) einsum(@(k) norm(u(1:2, k), 2) + n_sigma_99 * norm(S_k(:, (tri(k - 1) + 1):tri(k)), 2), 1:(N - 1)) * dt;

% Convex state path constraints
state_convex_constraints = {};

% Nonconvex state path constraints
state_nonconvex_constraints = {};

% Convex control constraints
lcvx_thrust_constraint = @(x, u, p, X_k, S_k) norm(u(1:2)) + n_sigma_99p9 * norm(S_k) - u_max;
control_convex_constraints = {lcvx_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex control constraints
control_nonconvex_constraints = {};

% Combine nonconvex constraints
nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

%% Set Up StochasticProblem
u_hold = "ZOH";

f = @(t, x, u, p) state_equation(x, u, mu, B);

f_0 = @(t, x, u, p) x; % Perfect measurment
g_0 = @(t, x, u, p) zeros([4, 4]);

guess.x = interp1(t, x', tspan_control)';
guess.u = interp1(t, u', tspan_control(1:(N - 1)))';
guess.p = [];

stoch_prob = StochasticProblem(x0, xf, Phat0, Ptilde0, Pf, N, u_hold, tf, f, G, f_0, g_0, guess, convex_constraints, stochastic_min_fuel_objective, nonconvex_constraints = nonconvex_constraints, scale = false, delta_t = delta_t);
[stoch_prob, Delta] = stoch_prob.discretize(stoch_prob.guess.x, stoch_prob.guess.u, stoch_prob.guess.p);

% Estimate X and S
K_k_guess = -1e1 * repmat([0, 0, 1, 0; 0, 0, 0, 1], 1, 1, stoch_prob.Nu);
Phat_k = zeros([stoch_prob.n.x, stoch_prob.n.x, stoch_prob.N]);
Phat_k(:, :, 1) = Phat0;
guess.X = zeros([4, stoch_prob.n.x * stoch_prob.N * (stoch_prob.N + 1) / 2]);
guess.X(1:4, 1:4, 1) = diag([sigma_pos, sigma_pos, sigma_vel, sigma_vel]);
guess.S = zeros([2, stoch_prob.n.x * stoch_prob.N * (stoch_prob.N + 1) / 2]);
for k = 1:stoch_prob.Nu
    Phat_k(:, :, k + 1) = (stoch_prob.disc.A_k(:, :, k) + stoch_prob.disc.B_k(:, :, k) * K_k_guess(:, :, k)) * Phat_k(:, :, k) * (stoch_prob.disc.A_k(:, :, k) + stoch_prob.disc.B_k(:, :, k) * K_k_guess(:, :, k))' + stoch_prob.disc.G_k(:, :, k) * stoch_prob.disc.G_k(:, :, k)';
    guess.X(1:4, 1:4, k + 1) = sqrtm(Phat_k(:, :, k + 1));
    guess.S(:, (tri(k - 1) + 1):tri(k)) = K_k_guess(:, :, k) * guess.X(:, (tri(k - 1) + 1):tri(k));
end

stoch_prob.guess = guess;

%% Solve Stochastic Optimization Problem with PTR
ptr_ops.iter_max = 20;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 1e2;
ptr_ops.w_tr = ones(1, stoch_prob.Nu) * 1e-2;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;
%stoch_ptr_sol = Stochastic_ptr_memoryefficient(stoch_prob, ptr_ops);

%stoch_ptr_sol.converged_i

%% MC Simulations with no Feedback
m = 20;

N_sub = 1;
numel_N_sub = N_sub * (N - 1) + 1;

t_no_fb = zeros([numel_N_sub, m]);
x_no_fb = zeros([stoch_prob.n.x, numel_N_sub, m]);
xhat_no_fb = zeros([stoch_prob.n.x, numel_N_sub, m]);
Phat_no_fb = zeros([stoch_prob.n.x, stoch_prob.n.x, numel_N_sub, m]);
u_no_fb = zeros([stoch_prob.n.u, numel_N_sub - 1, m]);

t_ref = t;
u_ref = @(t, x) interp1(t_ref, u', t)';

parfor i = 1:m
    [t_no_fb(:, i), x_no_fb(:, :, i), u_no_fb(:, :, i)] = stoch_prob.cont_prop_without_feedback_control(u_ref, stoch_prob.guess.p, N_sub = N_sub);
    %[t_no_fb(:, i), x_no_fb(:, :, i), xhat_no_fb(:, :, i), Phat_no_fb(:, :, :, i), u_no_fb(:, :, i)] = stoch_prob.disc_prop(stoch_prob.guess.x, u, stoch_prob.guess.p, K_k_guess, u_type = "cont", t_cont = t);
    i
end

%% Plot MC with No Correction Results
K_k_opt = -1e1 * repmat([0, 0, 1, 0; 0, 0, 0, 1], 1, 1, stoch_prob.Nu) * 0;
Phat_k = zeros([stoch_prob.n.x, stoch_prob.n.x, stoch_prob.N]);
Phat_k(:, :, 1) = Phat0;
for k = 1:stoch_prob.Nu
    Phat_k(:, :, k + 1) = (stoch_prob.disc.A_k(:, :, k) + stoch_prob.disc.B_k(:, :, k) * K_k_opt(:, :, k)) * Phat_k(:, :, k) * (stoch_prob.disc.A_k(:, :, k) + stoch_prob.disc.B_k(:, :, k) * K_k_opt(:, :, k))' + stoch_prob.disc.G_k(:, :, k) * stoch_prob.disc.G_k(:, :, k)';
end
%%

[P_eigvecs, P_eigvals] = pageeig(Phat_k(1:2, 1:2, :)); % Looks more correct then projecting ellipsoid...

thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
ellipse_3sigma = zeros([2, 100, stoch_prob.N]);
for k = 1:stoch_prob.N
    ellipse_3sigma(:, :, k) = stoch_prob.guess.x(1:2, k) + P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
end
figure
plot(squeeze(x_no_fb(1, :, :)), squeeze(x_no_fb(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
plot(stoch_prob.guess.x(1, :), stoch_prob.guess.x(2, :), Color = [30, 144, 255] / 256, LineWidth=1, DisplayName="Nominal"); hold on
x_lim = xlim;
plot(squeeze(ellipse_3sigma(1, :, 1)), squeeze(ellipse_3sigma(2, :, 1)), Color = "k", DisplayName="Covariance"); hold on
plot(squeeze(ellipse_3sigma(1, :, 2:end)), squeeze(ellipse_3sigma(2, :, 2:end)), Color = "k", HandleVisibility='off'); hold off
title("Without Optimized Feedback Policy")
xlabel("X [km]")
ylabel("Y [km]")
legend(location = "best")
axis equal
grid on

%%
earth_trajectory = earth_pos(t);
mars_trajectory = mars_pos(t);

figure
plot(squeeze(x_no_fb(1, :, :)), squeeze(x_no_fb(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
plot(x(1, :), x(2, :), DisplayName="Minimum Fuel Transfer", Color = [30, 144, 255] / 256, LineWidth=1); hold on
plot(earth_trajectory(1, :), earth_trajectory(2, :), DisplayName="Earth Orbit"); hold on;
plot(mars_trajectory(1, :), mars_trajectory(2, :), DisplayName="Mars Orbit"); hold on
plot(squeeze(ellipse_3sigma(1, :, 1)), squeeze(ellipse_3sigma(2, :, 1)), Color = "k", DisplayName="Covariance"); hold on
plot(squeeze(ellipse_3sigma(1, :, 2:end)), squeeze(ellipse_3sigma(2, :, 2:end)), Color = "k", HandleVisibility='off'); hold on
scatter(x0(1), x0(2), DisplayName="Earth at t = 0"); hold on
scatter(xf(1), xf(2), DisplayName=sprintf("Mars at t = %.g", tf)); hold off;
title("Minimum Fuel Transfer Between Earth and Mars with \sigma Increased by 30x")
subtitle(sprintf("For final time of %.g", tf))
xlabel("X [AU]")
ylabel("Y [AU]")
legend(location = "best")
grid on
axis equal
%%

tiledlayout(2,2)

nexttile
plot(t_no_fb(:, 1), squeeze(x_no_fb(1, :, :)), HandleVisibility='off'); hold off
title("r_{x_{MC}}")
xlabel("Time []")
ylabel("X Position [AU]")
grid on

nexttile
plot(t_no_fb(:, 1), squeeze(x_no_fb(2, :, :)), HandleVisibility='off'); hold on
title("r_{y_{MC}}")
xlabel("Time []")
ylabel("Y Position [AU]")
grid on

nexttile
plot(t_no_fb(:, 1), squeeze(x_no_fb(3, :, :)), HandleVisibility='off'); hold off
title("v_{x_{MC}}")
xlabel("Time []")
ylabel("X Velocity []")
grid on

nexttile
plot(t_no_fb(:, 1), squeeze(x_no_fb(4, :, :)), HandleVisibility='off'); hold off
title("v_{y_{MC}}")
xlabel("Time []")
ylabel("Y Velocity []")
grid on

sgtitle("Monte Carlo Simulations State History using \sigma Increased by 10x")

%%
x_ref = interp1(t_ref, x', t_no_fb(:, 1))';

tiledlayout(2,2)

nexttile
plot(t_no_fb(:, 1), x_ref(1, :)' - squeeze(x_no_fb(1, :, :)), HandleVisibility='off'); hold on
plot(tspan_control, squeeze(sqrt(Phat_k(1, 1, :))) * 3, Color = "k", LineStyle="--", DisplayName="3\sigma Bound"); hold on
plot(tspan_control, -squeeze(sqrt(Phat_k(1, 1, :))) * 3, Color = "k", LineStyle="--", HandleVisibility='off'); hold off
title("r_{x_{ref}} - r_{x_{MC}}")
xlabel("Time []")
ylabel("X Position Difference [AU]")
legend(Location="best")
grid on

nexttile
plot(t_no_fb(:, 1), x_ref(2, :)' - squeeze(x_no_fb(2, :, :)), HandleVisibility='off'); hold on
plot(tspan_control, squeeze(sqrt(Phat_k(2, 2, :))) * 3, Color = "k", LineStyle="--", DisplayName="3\sigma Bound"); hold on
plot(tspan_control, -squeeze(sqrt(Phat_k(2, 2, :))) * 3, Color = "k", LineStyle="--", HandleVisibility='off'); hold off
title("r_{y_{ref}} - r_{y_{MC}}")
xlabel("Time []")
ylabel("Y Position Difference [AU]")
legend(Location="best")
grid on

nexttile
plot(t_no_fb(:, 1), x_ref(3, :)' - squeeze(x_no_fb(3, :, :)), HandleVisibility='off'); hold on
plot(tspan_control, squeeze(sqrt(Phat_k(3, 3, :))) * 3, Color = "k", LineStyle="--", DisplayName="3\sigma Bound"); hold on
plot(tspan_control, -squeeze(sqrt(Phat_k(3, 3, :))) * 3, Color = "k", LineStyle="--", HandleVisibility='off'); hold off
title("v_{x_{ref}} - v_{x_{MC}}")
xlabel("Time []")
ylabel("X Velocity Difference []")
legend(Location="best")
grid on

nexttile
plot(t_no_fb(:, 1), x_ref(4, :)' - squeeze(x_no_fb(4, :, :)), HandleVisibility='off'); hold on
plot(tspan_control, squeeze(sqrt(Phat_k(4, 4, :))) * 3, Color = "k", LineStyle="--", DisplayName="3\sigma Bound"); hold on
plot(tspan_control, -squeeze(sqrt(Phat_k(4, 4, :))) * 3, Color = "k", LineStyle="--", HandleVisibility='off'); hold off
title("v_{y_{ref}} - v_{y_{MC}}")
xlabel("Time []")
ylabel("Y Velocity Difference []")
legend(Location="best")
grid on

sgtitle("MC State Difference from Reference using \sigma Increased by 30x")

%%
covariance_plot_no_target(x(:, end), squeeze(x_no_fb(:, end, :)), Phat_k(:,:,end), ["x [km]", "y [km]", "v_x [km / s]", "v_y [km / s]"], "State Dispersion at Final Node for \sigma Increased by 30x")


%% Run MC Correction Results
stoch_ptr_sol.converged_i = 1;
K_k_opt = recover_gain_matrix(stoch_ptr_sol.X, stoch_ptr_sol.S);

K_k_ck = zeros([stoch_prob.N, 1]);
for km = 1:(stoch_prob.N - 1)
    K_k_ck(km) = norm(K_k_opt(:, :, km) * stoch_ptr_sol.X(:, (tri(km - 1) + 1):tri(km)) - stoch_ptr_sol.S(:, (tri(km - 1) + 1):tri(km)));
end

m = 20;

N_sub = 2;
numel_N_sub = N_sub * (N - 1) + 1;

t_ofb = zeros([numel_N_sub, m]);
x_ofb = zeros([stoch_prob.n.x, numel_N_sub, m]);
u_ofb = zeros([stoch_prob.n.u, numel_N_sub - 1, m]);

parfor i = 1:m
    %[t_no_fb(:, i), x_no_fb(:, :, i), u_no_fb(:, :, i)] = stoch_prob.cont_prop_without_feedback_control(stoch_prob.guess.u, stoch_prob.guess.p);
    %[t_ofb(:, i), x_ofb(:, :, i), xhat_ofb(:, :, i), Phat_ofb(:, :, :, i), u_ofb(:, :, i)] = stoch_prob.disc_prop(stoch_ptr_sol.x, stoch_ptr_sol.u, stoch_ptr_sol.p(:), K_k_opt, N_sub = N_sub);
    [t_ofb(:, i), x_ofb(:, :, i), u_ofb(:, :, i)] = cont_prop_feedback_no_kalman_filter(stoch_prob, stoch_ptr_sol.x, stoch_ptr_sol.u, stoch_ptr_sol.p, K_k_opt, N_sub = N_sub);
    i
end

%%
AU = 149597898; % [km]
mu_sun = 132712440017.99; % [km3 / s2]
mean_motion = sqrt(mu_sun / AU ^ 3); % [1 / s]

cost_fb = max(squeeze(sum(vecnorm(u_ofb, 2, 1), 2) * (t_ofb(2, 1) - t_ofb(1, 1)))) * AU * mean_motion

cost_no_fb = mean(squeeze(sum(vecnorm(u_no_fb, 2, 1), 2) * dt)) * AU * mean_motion

%% Plot MC with Correction Results

[t_mean, x_mean, u_mean] = cont_prop_feedback_no_kalman_filter(stoch_prob, stoch_ptr_sol.x, stoch_ptr_sol.u, stoch_ptr_sol.p, K_k_opt * 0, N_sub = N_sub, noise_multiplier = 0, x_0 = stoch_prob.x0);
%%
X_k = stoch_ptr_sol.X;
S_k = stoch_ptr_sol.S;
[P_k, Pu_k] = recover_est_covariances(X_k, S_k);

[P_eigvecs, P_eigvals] = pageeig(P_k(1:2, 1:2, :)); % Looks more correct then projecting ellipsoid...

thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
ellipse_3sigma_opt = zeros([2, 100, stoch_prob.N]);
ellipse_3sigma_normal = zeros([2, 100, stoch_prob.N]);
for k = 1:stoch_prob.N
    ellipse_3sigma_opt(:, :, k) = stoch_ptr_sol.x(1:2, k) + 10 * P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
    ellipse_3sigma_normal(:, :, k) = stoch_ptr_sol.x(1:2, k) + P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
end

tiledlayout(1, 2,"TileSpacing","tight")

nexttile
plot(x_mean(1, :)' + 10 * (squeeze(x_ofb(1, :, :)) - x_mean(1, :)'), x_mean(2, :)' + 10 * (squeeze(x_ofb(2, :, :)) - x_mean(2, :)'), Color = [192, 192, 192] / 256,LineWidth = 1,  HandleVisibility='off'); hold on
plot(squeeze(x_mean(1, :)), squeeze(x_mean(2, :)), Color = [30, 144, 255] / 256, LineWidth=1, DisplayName="Nominal"); hold on
plot(earth_trajectory(1, :), earth_trajectory(2, :), DisplayName="Earth Orbit"); hold on;
plot(mars_trajectory(1, :), mars_trajectory(2, :), DisplayName="Mars Orbit"); hold on
x_lim = xlim;
plot(squeeze(ellipse_3sigma_opt(1, :, 1)), squeeze(ellipse_3sigma_opt(2, :, 1)), Color = "k", DisplayName="Covariance"); hold on
plot(squeeze(ellipse_3sigma_opt(1, :, 2:2:end)), squeeze(ellipse_3sigma_opt(2, :, 2:2:end)), Color = "k", HandleVisibility='off'); hold on
scatter(x0(1), x0(2), DisplayName="Earth at t = 0"); hold on
scatter(xf(1), xf(2), DisplayName=sprintf("Mars at t = %.g", tf)); hold off;
title("With Optimized Feedback Policy")
text(-0.27, 2.2, "Position Deviation Scaled by 10")
xlabel("X [AU]")
ylabel("Y [AU]")
legend(location = "best")
axis equal
grid on

x_mean_no = interp1(t, x', t_no_fb(:, 1))';

nexttile
plot(x_mean_no(1, :)' + 2 * (squeeze(x_no_fb(1, :, :)) - x_mean_no(1, :)'), x_mean_no(2, :)' + 2 * (squeeze(x_no_fb(2, :, :)) - x_mean_no(2, :)'), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
plot(x(1, :), x(2, :), DisplayName="Minimum Fuel Transfer", Color = [30, 144, 255] / 256, LineWidth=1); hold on
plot(earth_trajectory(1, :), earth_trajectory(2, :), DisplayName="Earth Orbit"); hold on;
plot(mars_trajectory(1, :), mars_trajectory(2, :), DisplayName="Mars Orbit"); hold on
plot(squeeze(ellipse_3sigma(1, :, 1)), squeeze(ellipse_3sigma(2, :, 1)), Color = "k", DisplayName="Covariance"); hold on
plot(squeeze(ellipse_3sigma(1, :, 2:2:end)), squeeze(ellipse_3sigma(2, :, 2:2:end)), Color = "k", HandleVisibility='off'); hold on
scatter(x0(1), x0(2), DisplayName="Earth at t = 0"); hold on
scatter(xf(1), xf(2), DisplayName=sprintf("Mars at t = %.g", tf)); hold off;
text(-1, 2.48, "Position Deviation Scaled by 2")
title("No Feedback Policy")
xlabel("X [AU]")
ylabel("Y [AU]")
%legend(location = "best")
grid on
axis equal

% Plot zoomed in terminal section
proj_Pf_r = project_ellipsoid(Pf, [1,2]);
[Pf_eigvecs, Pf_eigvals] = eigs(proj_Pf_r);
terminal_ellipse = xf(1:2) + 10 * Pf_eigvecs * [3 * sqrt(Pf_eigvals(1, 1)) * cos(thetas); 3 * sqrt(Pf_eigvals(2, 2)) * sin(thetas)];

axes('Position',[.3 .16 .09 .14])
box on
plot(squeeze(x_ofb(1, :, :)), squeeze(x_ofb(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
plot(x_mean(1, :), x_mean(2, :), Color = [30, 144, 255] / 256, LineWidth=1, DisplayName="Nominal"); hold on
plot(terminal_ellipse(1, :), terminal_ellipse(2, :), Color = "k", HandleVisibility='off', LineWidth=1, LineStyle="-."); hold on
plot(squeeze(ellipse_3sigma_normal(1, :, end)), squeeze(ellipse_3sigma_normal(2, :, end)), Color = "k", DisplayName="Solution", LineWidth=0.7); hold off
xlim(xf(1) + 2 * [-max(3 * sqrt(Pf_eigvals(1, 1)), 3 * sqrt(Pf_eigvals(2, 2))), max(3 * sqrt(Pf_eigvals(1, 1)), 3 * sqrt(Pf_eigvals(2, 2)))])
ylim(xf(2) + 2 * [-max(3 * sqrt(Pf_eigvals(1, 1)), 3 * sqrt(Pf_eigvals(2, 2))), max(3 * sqrt(Pf_eigvals(1, 1)), 3 * sqrt(Pf_eigvals(2, 2)))])
x_lim = xlim;
grid on

earth_trajectory = earth_pos(t);
mars_trajectory = mars_pos(t);

sgtitle(sprintf("Stochastic Optimized Minimum Fuel Transfer Between Earth and Mars\n"))

%% States
titles = ["X Position", "Y Position", "X Velocity", "Y Velocity"];
ylabels = ["r_x [AU]", "r_y [AU]", "v_x []", "v_y []"];

ops = {@(x) x, @(x) x, @(x) x, @(x) x};

figure 
tiledlayout(2,2)
for xi = 1:4
    %x_3sigbound = 3 * sqrt(squeeze(project_ellipsoid(P_k, x)))';
    x_3sigbound = 3 * sqrt(squeeze(P_k(xi, xi, :)))';
    x_3sigbound_cont = x_3sigbound;
        
    nexttile
    for i = 1:m
        plot(t_ofb(:, 1), ops{xi}(x_mean(xi, :)) - ops{xi}(x_ofb(xi, :, i)), HandleVisibility='off'); hold on
    end
    plot(tspan_control, x_3sigbound_cont, LineStyle="--", Color = "k", LineWidth=1, DisplayName="3\sigma Bound"); hold on
    plot(tspan_control, -x_3sigbound_cont, LineStyle="--", Color = "k", LineWidth=1, HandleVisibility='off'); hold off
    title(titles(xi) + " Difference vs Time")
    xlabel("Time []")
    ylabel(ylabels(xi))
    legend("Location","southeast")
    grid on
end
sgtitle("MC State Difference from Reference using \sigma_{dist} with Optimized Feedback Policies")

%% Thrust Magnitude
u_mean = stoch_ptr_sol.u;


figure
tiledlayout(1, 3)

%Pu_k = pagemtimes(S_k, pagetranspose(S_k));
%thrust_3sigbound = 3 * sqrt(squeeze(project_ellipsoid(Pu_k, 3)))';
nu = 2;
epsilon = 1e-3; % 99.9% 
S_k_norm = zeros([N - 1, 1]);
thrustx_3sigbound = zeros([N - 1, 1]);
thrusty_3sigbound = zeros([N - 1, 1]);
for j = 1:(N - 1)
    S_k_norm(j) = norm(S_k(:, (tri(j - 1) + 1):tri(j)));
    thrustx_3sigbound(j) = norminv(1 - 1e-3) * norm(S_k(1, (tri(j - 1) + 1):tri(j)));
    thrusty_3sigbound(j) = norminv(1 - 1e-3) * norm(S_k(2, (tri(j - 1) + 1):tri(j)));
end
thrust_3sigbound = squeeze(sigma_mag_confidence(epsilon / 2, nu) * S_k_norm);
u_mean_full = interp1(tspan_control(1:size(u_mean, 2)), u_mean', tspan_control, "previous", "extrap")';
thrust_3sigbound_full = interp1(tspan_control(1:size(u_mean, 2)), thrust_3sigbound, tspan_control, "previous", "extrap");
thrustx_3sigbound_full = interp1(tspan_control(1:size(u_mean, 2)), thrustx_3sigbound, tspan_control, "previous", "extrap");
thrusty_3sigbound_full = interp1(tspan_control(1:size(u_mean, 2)), thrusty_3sigbound, tspan_control, "previous", "extrap");

% Thrust x
nexttile
for i = 1:m
    stairs(t_ofb(1:(end-1), 1), u_ofb(1, :, i), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
end

stairs(t, u(1, :), Color = "m",LineWidth=1, DisplayName="Reference"); hold on
stairs(tspan_control, u_mean_full(1, :), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
stairs(tspan_control, u_mean_full(1, :) + thrustx_3sigbound_full, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
stairs(tspan_control, u_mean_full(1, :) - thrustx_3sigbound_full, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1,HandleVisibility='off'); hold off
title("Thrust in X")
%legend(Location="southeast")
xlabel("Time [s]")
ylabel("Thrust [N]")
grid on

% Thrust y
nexttile
for i = 1:m
    stairs(t_ofb(1:(end-1), 1), u_ofb(2, :, i), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
end

stairs(t, u(2, :), Color = "m",LineWidth=1, DisplayName="Reference"); hold on
stairs(tspan_control, u_mean_full(2, :), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
stairs(tspan_control, u_mean_full(2, :) + thrusty_3sigbound_full, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
stairs(tspan_control, u_mean_full(2, :) - thrusty_3sigbound_full, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1,HandleVisibility='off'); hold off
title("Thrust in Y")
legend(Location="southeast")
xlabel("Time [s]")
ylabel("Thrust [N]")
grid on

% Thrust mag
nexttile
for i = 1:m
    stairs(t_ofb(1:(end-1), 1), vecnorm(u_ofb(1:2, :, i), 2, 1), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
end

stairs(t, vecnorm(u, 2, 1)', Color = "m",LineWidth=1, DisplayName="Reference"); hold on
stairs(tspan_control, vecnorm(u_mean_full(1:2, :),2,1), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
stairs(tspan_control, (vecnorm(u_mean_full(1:2, :),2,1) + thrust_3sigbound_full), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
stairs(tspan_control, max((vecnorm(u_mean_full(1:2, :),2,1) - thrust_3sigbound_full), 0), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1,HandleVisibility='off'); hold on
yline(u_max, LineWidth = 1, LineStyle="--", Color="k", DisplayName = "Constraint"); hold on
title("Thrust Magnitude")
%legend(Location="southeast")
xlabel("Time [s]")
ylabel("Thrust [N]")
grid on

sgtitle("Control vs Time with Optimized Feedback Policies")

%%
figure
covariance_plot(x_mean(:, end), squeeze(x_ofb(:, end, :)), squeeze(P_k(:, :, end)), Phat0 + Ptilde0, ["x [km]", "y [km]", "v_x [km / s]", "v_y [km / s]"], "State Dispersion at Final Node with Optimal Feedback")


%% Helper Functions
function [xdot] = state_equation(x, u, mu, B)
    rvec = x(1:2);
    vvec = x(3:4);
    r = sqrt(rvec(1) ^ 2 + rvec(2) ^ 2);
    
    f0 = [vvec; -mu / r ^ 3 * rvec];

    xdot = f0 + B * u;
end