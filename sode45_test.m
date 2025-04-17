%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 15 April, 2025
% Description: Test of stochastic propagation without optimized feedback 
% control
% Most Recent Change: 15 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do HW 9 Q1d
f = @(t, x, u, p) [0; 0];
u = @(t, x) [0; 0];
p = 0;

sigma_1 = 2;
sigma_2 = 3;
G = @(t, x, u, p) [sigma_1, 0; 0, sigma_2];

x0 = [0; 0];

w = @(n) randn([2, n]);
delta_t = 1e-3;

tspan = 0:delta_t:1;

m = 200;

x = zeros([numel(tspan), 2, m]);
t = zeros([numel(tspan), m]);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

parfor i = 1:m
    [t(:, i), x(:, :, i)] = sode45(f, G, u, p, w, tspan, delta_t, x0, tolerances);
    i
end
%% Analyze distribution of trajectories to check if it matches expectations
tiledlayout(2, 2);

nexttile
plot(t, squeeze(x(:, 1, :))); hold on
plot(t, 3 * sigma_1 * sqrt(t), Color = "k"); hold on
plot(t, -3 * sigma_1 * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_1")
legend("x_1", "", "3 \sigma Bound")
title("x_1 Trajectories")

nexttile
plot(t, squeeze(x(:, 2, :))); hold on
plot(t, 3 * sigma_2 * sqrt(t), Color = "k"); hold on
plot(t, -3 * sigma_2 * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_2")
legend("x_2", "", "3 \sigma Bound")
title("x_2 Trajectories")

% Estimate standard deviation over time normalized by sqrt(t) to check if
% it matches sigma_1, sigma_2
stds = std(x,0,3) ./ sqrt(t(:, 1));

nexttile
histogram(stds(:, 1))
xlabel("Normalized x_1")
title("x_1 Standard Deviation Distribution Normalized by sqrt(t)")

nexttile
histogram(stds(:, 2))
xlabel("Normalized x_2")
title("x_2 Standard Deviation Distribution Normalized by sqrt(t)")

%%
xs = squeeze(x(:, 2, :));
hist3([t(:), xs(:)], [30, 30],'CdataMode','auto');

%% Do deterministic 3DoF
Deterministic_3DoF

%% Monte Carlo
t_k = linspace(0, prob_3DoF.tf, prob_3DoF.N);
if prob_3DoF.u_hold == "ZOH"
    u_func = @(t, x) interp1(t_k(1:prob_3DoF.Nu), ptr_sol.u(:, :, ptr_sol.converged_i)', t, "previous", "extrap")';
elseif prob_3DoF.u_hold == "FOH"
    u_func = @(t, x) interp1(t_k(1:prob_3DoF.Nu), ptr_sol.u(:, :, ptr_sol.converged_i)', t)';
end

%%
sigma_x0 = [10e-3; ... % r_x
            50e-3; ... % r_y
            5e-3; ... % v_x
            5e-3; ... % v_y
            deg2rad(0.2); ... % theta
            deg2rad(0.02)];  % w

sigma_accelx = 0.5e-3 * 0;
sigma_accely = 0.1e-3 * 0;
sigma_alpha = deg2rad(2e-2) * 0;
G = @(t, x, u, p) [[zeros([2, 3]); ... % velocity
                   [sigma_accelx; sigma_accely] .* eye([2, 3]); ... % acceleration
                   zeros([1, 3]); ... % angular velocity
                   sigma_alpha * [0, 0, 1]], ... % angular acceleration
                   zeros(6, 3)]; 

w = @(n) [randn([3, n]); zeros([3, n])];
delta_t = 1e-0;

tspan = 0:delta_t:prob_3DoF.tf;
tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(ptr_sol.u(:, :, ptr_sol.converged_i), ptr_sol.p(:, ptr_sol.converged_i), tspan = tspan);

m = 20;

x_0_m = prob_3DoF.x0 + sigma_x0 .* randn([prob_3DoF.n.x, m]);

t = zeros([numel(tspan), m]);
x = zeros([prob_3DoF.n.x, numel(tspan), m]);
u = zeros([prob_3DoF.n.u, numel(tspan), m]);

parfor i = 1:m
    [t(:, i), x(:, :, i), u_applied(:, :, i)] = sode45(prob_3DoF.cont.f, G, u_func, guess.p, w, tspan, delta_t, x_0_m(:, i), tolerances);
    i
end
%%
N_sub = 15;
K_k = -3*repmat([0, 1, 0, 0, 0, 0; 0, 0, 0, 0, -1, 0; 0, 1, 0, 0, 1, 0], 1, 1, prob_3DoF.Nu);

t_fb = zeros([N_sub * (numel(t_k) - 1) + 1, m]);
x_fb = zeros([prob_3DoF.n.x, N_sub * (numel(t_k) - 1) + 1, m]);
u_fb = zeros([prob_3DoF.n.u, N_sub * (numel(t_k) - 1) + 1, m]);

parfor i = 1:m
    [t_fb(:, i), x_fb(:, :, i), u_fb(:, :, i)] = propagate_cont_feedback_no_kalman_filter(x_0_m(:, i), ptr_sol.x(:, :, ptr_sol.converged_i), ptr_sol.u(:, :, ptr_sol.converged_i), K_k, prob_3DoF.cont.f, G, t_k, N_sub, w, delta_t, tolerances);
    i
end
%%
[A_k, B_k, E_k, c_k, G_k, Delta] = discretize_stochastic_dynamics_ZOH(prob_3DoF.cont.f, prob_3DoF.cont.A, prob_3DoF.cont.B, prob_3DoF.cont.E, prob_3DoF.cont.c, G, prob_3DoF.N, t_k, ptr_sol.x(:, :, ptr_sol.converged_i), ptr_sol.u(:, :, ptr_sol.converged_i), ptr_sol.p(:, ptr_sol.converged_i), prob_3DoF.tolerances);
%%
P_k = zeros([prob_3DoF.n.x, prob_3DoF.n.x, numel(t_k)]);
P_k(:, :, 1) = diag(sigma_x0 .^ 2);
for k = 1:(numel(t_k) - 1)
    P_k(:, :, k + 1) = (A_k(:, :, k) + B_k(:, :, k) * K_k(:, :, k)) * P_k(:, :, k) * (A_k(:, :, k) + B_k(:, :, k) * K_k(:, :, k))' + G_k(:, :, k) * G_k(:, :, k)';
end

%%

proj_P_r = project_ellipsoid(P_k, [1,2]);

%%
proj_P_rx = project_ellipsoid(P_k, 4)

%%

[P_eigvecs, P_eigvals] = pageeig(proj_P_r);

X_k = zeros(size(P_k));
thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
ellipse_3sigma = zeros([2, 100, numel(t_k)]);
for k = 1:numel(t_k)
    ellipse_3sigma(:, :, k) = ptr_sol.x(1:2, k, ptr_sol.converged_i) + P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
    X_k(:, :, k) = chol(P_k(:, :, k), "lower");
end

figure
plot(squeeze(x(1, :, :)), squeeze(x(2, :, :)), Color = [192, 192, 192] / 256); hold on
plot(x_cont_sol(1, :), x_cont_sol(2, :), Color = [30, 144, 255] / 256, LineWidth=1); hold on
plot(squeeze(ellipse_3sigma(1, :, :)), squeeze(ellipse_3sigma(2, :, :)), Color = "k"); hold off
title("Monte Carlo Simulation of 3DoF Rocket Landing with No Feedback Control")
xlabel("X [km]")
ylabel("Y [km]")
grid on

%%
[P_eigvecs, P_eigvals] = pageeig(proj_P_r);

X_k = zeros(size(P_k));
thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
ellipse_3sigma = zeros([2, 100, numel(t_k)]);
for k = 1:numel(t_k)
    ellipse_3sigma(:, :, k) = ptr_sol.x(1:2, k, ptr_sol.converged_i) + P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
    X_k(:, :, k) = chol(P_k(:, :, k), "lower");
end

figure
plot(squeeze(x_fb(1, :, :)), squeeze(x_fb(2, :, :)), Color = [192, 192, 192] / 256); hold on
plot(x_cont_sol(1, :), x_cont_sol(2, :), Color = [30, 144, 255] / 256, LineWidth=1); hold on
plot(squeeze(ellipse_3sigma(1, :, :)), squeeze(ellipse_3sigma(2, :, :)), Color = "k"); hold off
title("Monte Carlo Simulation of 3DoF Rocket Landing with No Feedback Control")
xlabel("X [km]")
ylabel("Y [km]")
grid on

%%
S_k = pagemtimes(K_k, X_k(:, :, 1:prob_3DoF.Nu));
plot_3DoF_MC_time_histories(tspan, x_cont_sol, ptr_sol.u(:, :, ptr_sol.converged_i), t_fb(:, 1), x_fb, u_fb, t_k, X_k, S_k, t(:, 1), x, vehicle.max_thrust, vehicle.min_thrust, vehicle.max_gimbal, true)

%%
plot_3DoF_MC_trajectories(tspan, x_cont_sol, t_fb(:, 1), x_fb, ptr_sol.x(:, :, ptr_sol.converged_i), t_k, X_k, t(:, 1), x, P_k(:, :, 1)/10, glideslope_angle_max)

%% Constraint Histogram

glideslope_func = @(x, i) max(squeeze(acosd(x(2, 1:(end-i), :) ./ vecnorm(x(1:2, 1:(end-i), :), 2, 1))), [], 1);
plot_MC_constraint_hist(glideslope_func, [100, 10], x_fb, x, rad2deg(glideslope_angle_max), "Max Glideslope Angle", "Max Glideslope Angle [deg]", 7)

%% Analyze distribution of trajectories to check if it matches expectations
tiledlayout(2, 2);

nexttile
plot(t(:, 1), squeeze(x(1, :, :))); hold on
plot(t(:, 1), sigma_accelx * sqrt(t(:, 1)), Color = "k"); hold on
plot(t(:, 1), -sigma_accelx * sqrt(t(:, 1)), Color = "k"); hold off
xlabel("Time")
ylabel("x_1")
legend("x_1", "", "3 \sigma Bound")
title("x_1 Trajectories")

nexttile
plot(t, squeeze(x(2, :, :))); hold on
plot(t, 3 * sigma_alpha * sqrt(t), Color = "k"); hold on
plot(t, -3 * sigma_alpha * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_2")
legend("x_2", "", "3 \sigma Bound")
title("x_2 Trajectories")

% Estimate standard deviation over time normalized by sqrt(t) to check if
% it matches sigma_1, sigma_2
stds = std(x,0,3) ./ sqrt(t(:, 1))';

nexttile
histogram(stds(1, :), 10)
xlabel("Normalized x_1")
title("x_1 Standard Deviation Distribution Normalized by sqrt(t)")

nexttile
histogram(stds(2, :), 10)
xlabel("Normalized x_2")
title("x_2 Standard Deviation Distribution Normalized by sqrt(t)")

%%
xs = squeeze(x(1, :, :));
hist3([t(:), xs(:)], [30, 30],'CdataMode','auto');
view(2)