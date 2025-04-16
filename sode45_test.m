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
sigma_x0 = [50e-3; ... % r_x
            50e-3; ... % r_y
            5e-3; ... % v_x
            5e-3; ... % v_y
            deg2rad(2); ... % theta
            deg2rad(0.2)];  % w

sigma_accelx = 0.2e-3;
sigma_accely = 0.1e-3;
sigma_alpha = deg2rad(2e-2);
G = @(t, x, u, p) [[zeros([2, 3]); ... % velocity
                   [sigma_accelx; sigma_accely] .* eye([2, 3]); ... % acceleration
                   zeros([1, 3]); ... % angular velocity
                   sigma_alpha * [0, 0, 1]], ... % angular acceleration
                   zeros(6, 3)]; 

w = @(n) [randn([3, n]); zeros([3, n])];
delta_t = 1e-1;

tspan = 0:delta_t:prob_3DoF.tf;
tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(ptr_sol.u(:, :, ptr_sol.converged_i), ptr_sol.p(:, ptr_sol.converged_i), tspan = tspan);

m = 20;

x_0_m = prob_3DoF.x0 + sigma_x0 .* randn([prob_3DoF.n.x, m]);

t = zeros([numel(tspan), m]);
x = zeros([numel(tspan), prob_3DoF.n.x, m]);
u = zeros([numel(tspan), prob_3DoF.n.u, m]);

parfor i = 1:m
    [t(:, i), x(:, :, i), u_applied(:, :, i)] = sode45(prob_3DoF.cont.f, G, u_func, guess.p, w, tspan, delta_t, x_0_m(:, i), tolerances);
    i
end
%%
[A_k, B_k, E_k, c_k, G_k, Delta] = discretize_stochastic_dynamics_ZOH(prob_3DoF.cont.f, prob_3DoF.cont.A, prob_3DoF.cont.B, prob_3DoF.cont.E, prob_3DoF.cont.c, G, prob_3DoF.N, t_k, ptr_sol.x(:, :, ptr_sol.converged_i), ptr_sol.u(:, :, ptr_sol.converged_i), ptr_sol.p(:, ptr_sol.converged_i), prob_3DoF.tolerances);
%%
K_k = zeros([prob_3DoF.n.u, prob_3DoF.n.x, numel(t_k)]);

P_k = zeros([prob_3DoF.n.x, prob_3DoF.n.x, numel(t_k)]);
P_k(:, :, 1) = diag(sigma_x0 .^ 2);
for k = 1:(numel(t_k) - 1)
    P_k(:, :, k + 1) = (A_k(:, :, k) + B_k(:, :, k) * K_k(:, :, k)) * P_k(:, :, k) * (A_k(:, :, k) + B_k(:, :, k) * K_k(:, :, k))' + G_k(:, :, k) * G_k(:, :, k)';
end

%%
[P_eigvals, P_eigvecs] = pageeig(P_k);


figure
plot(squeeze(x(:, 1, :)), squeeze(x(:, 2, :)))
title("Monte Carlo Simulation of 3DoF Rocket Landing with No Feedback Control")
xlabel("X [km]")
ylabel("Y [km]")
grid on

%%
plot_3DoF_MC_time_histories(tspan, x_cont_sol, u_cont_sol, x, u_applied, 0, x)

%% Analyze distribution of trajectories to check if it matches expectations
tiledlayout(2, 2);

nexttile
plot(t, squeeze(x(:, 1, :))); hold on
plot(t, sigma_accelx * sqrt(t), Color = "k"); hold on
plot(t, -sigma_accelx * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_1")
legend("x_1", "", "3 \sigma Bound")
title("x_1 Trajectories")

nexttile
plot(t, squeeze(x(:, 2, :))); hold on
plot(t, 3 * sigma_alpha * sqrt(t), Color = "k"); hold on
plot(t, -3 * sigma_alpha * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_2")
legend("x_2", "", "3 \sigma Bound")
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