%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 15 April, 2025
% Description: Test of stochastic propagation without feedback control
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

m = 20;

x = zeros([numel(tspan), 2, m]);
t = zeros([numel(tspan), m]);

parfor i = 1:m
    [t(:, i), x(:, :, i)] = sode45(f, G, u, p, w, tspan, delta_t, x0, tolerances);
    i
end
%% Analyze distribution of trajectories to check if it matches expectations
tiledlayout(2, 2);

nexttile
plot(t, squeeze(x(:, 1, :))); hold on
plot(t, sigma_1 * sqrt(t), Color = "k"); hold on
plot(t, -sigma_1 * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_1")
legend("x_1", "", string(sigma_1) + " \sigma Bound")
title("x_1 Trajectories")

nexttile
plot(t, squeeze(x(:, 2, :))); hold on
plot(t, sigma_2 * sqrt(t), Color = "k"); hold on
plot(t, -sigma_2 * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_2")
legend("x_2", "", string(sigma_2) + " \sigma Bound")
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
sigma_accelx = 0.2e-3;
sigma_accely = 0.1e-3;
sigma_alpha = 3e-4;
G = @(t, x, u, p) [zeros([2, 3]); ... % velocity
                   [sigma_accelx; sigma_accely] .* eye([2, 3]); ... % acceleration
                   zeros([1, 3]); ... % angular velocity
                   sigma_alpha * [0, 0, 1]]; % angular acceleration

w = @(n) randn([3, n]);
delta_t = 1e-0;

tspan = 0:delta_t:tf;

m = 100;

x = zeros([numel(tspan), prob_3DoF.n.x, m]);
t = zeros([numel(tspan), m]);

parfor i = 1:m
    [t(:, i), x(:, :, i)] = sode45(prob_3DoF.cont.f, G, u_func, guess.p, w, tspan, delta_t, x_0, tolerances);
    i
end
%%
figure
plot(squeeze(x(:, 1, :)), squeeze(x(:, 2, :)))
title("Monte Carlo Simulation of 3DoF Rocket Landing with No Feedback Control")
xlabel("X [km]")
ylabel("Y [km]")
grid on

%% Analyze distribution of trajectories to check if it matches expectations
tiledlayout(2, 2);

nexttile
plot(t, squeeze(x(:, 1, :))); hold on
plot(t, sigma_accelx * sqrt(t), Color = "k"); hold on
plot(t, -sigma_accelx * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_1")
legend("x_1", "", string(sigma_accelx) + " \sigma Bound")
title("x_1 Trajectories")

nexttile
plot(t, squeeze(x(:, 2, :))); hold on
plot(t, sigma_alpha * sqrt(t), Color = "k"); hold on
plot(t, -sigma_alpha * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_2")
legend("x_2", "", string(sigma_alpha) + " \sigma Bound")
title("x_2 Trajectories")

% Estimate standard deviation over time normalized by sqrt(t) to check if
% it matches sigma_1, sigma_2
stds = std(x,0,3) ./ sqrt(t(:, 1));

nexttile
histogram(stds(:, 1), 10)
xlabel("Normalized x_1")
title("x_1 Standard Deviation Distribution Normalized by sqrt(t)")

nexttile
histogram(stds(:, 2), 10)
xlabel("Normalized x_2")
title("x_2 Standard Deviation Distribution Normalized by sqrt(t)")

%%
xs = squeeze(x(:, 1, :));
hist3([t(:), xs(:)], [30, 30],'CdataMode','auto');
view(2)