%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath(pwd));

%% Initialize
% Vehicle Parameters
u_max = 0.1;
mu = 1;
a_Earth = 1; a_Mars = 1.524;
theta0_Earth = 0; theta0_Mars = pi;
omega_Earth = sqrt(mu / a_Earth^3);
omega_Mars = sqrt(mu / a_Mars^3);

% Problem Parameters
tf = 10; % [s]
N = 30; % []
r_0 = [a_Earth; 0];
v_0 = [0; sqrt(mu/a_Earth)];

P_earth_over_P_mars = (a_earth / a_mars) ^ (3/2);
mars_pos = @(t) a_mars .* [cos(t * P_earth_over_P_mars + nu0_mars); ...
                           sin(t * P_earth_over_P_mars + nu0_mars)];
thetaf_Mars = @(tf) tf * P_earth_over_P_mars + theta0_Mars;

x_0 = [r_0; v_0];
x_f = [];% [r_f; v_f];

tspan = [0, 1];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

parser = "CVX";

% PTR algorithm parameters
ptr_ops.iter_max = 100;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 1e-3;
ptr_ops.w_vc = 1e2;
ptr_ops.w_tr = ones(1, Nu) * 1e-4;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = false;

nx = 4;
nu = 2;
np = 1;

%% Get Dynamics
f = @(t, x, u, p) dynamics_asteroidmining(t, x, u, p);

%% Specify Constraints
% Convex state path constraints

% Convex control constraints
max_thrust_constraint = @(t, x, u, p) norm(u, 2)-u_max;
control_convex_constraints = {max_thrust_constraint};

% Combine convex constraints
convex_constraints = control_convex_constraints;

% Terminal boundary conditions
mars_xf = @(tf) [mars_pos(tf); v_circ(mars_pos(tf), thetaf_Mars(tf), mu)];
terminal_bc = @(t, x, u, p) x - mars_xf(p(1));
terminal_bc_linearized = linearize_constraint(terminal_bc, nx, nu, np, "p", 1);
terminal_bc_linearized = @(x, p, x_ref, p_ref) terminal_bc_linearized(0, x, 0, p, x_ref, 0, p_ref);

%% Specify Objective
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p, x_ref, u_ref, p_ref) sum(norms(u, 2, 1)) * p_ref/(N + 1) + sum(norms(u_ref, 2, 1)) * (p-p_ref)/(N + 1);
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p) sum((u(3, 1:(end - 1)) + u(3, 2:end)) / 2) * delta_t;
end

%% Create Guess
AU_guess = interp1(tspan, [a_Earth, a_Mars]', t_k);
theta_guess = interp1(tspan, [theta0_Earth, thetaf_Mars(tf)]', t_k);
r_guess = [AU_guess.*cos(theta_guess); AU_guess.*sin(theta_guess)];
v_guess = v_circ_array(r_guess, theta_guess, mu);

function [v_guess] = v_circ_array(r_guess, theta_guess, mu)
    r = vecnorm(r_guess, 2, 1);
    v = sqrt(mu ./ r);
    v_guess = v .* [-sin(theta_guess); cos(theta_guess)];
end

guess.x = [r_guess; v_guess];
guess.u = interp1(tspan, zeros(2, 2)', t_k(1:Nu))';
guess.p = tf;

%% Construct Problem Object
problem = DeterministicProblem(x_0, x_f, N, u_hold, 1, f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc_linearized, integration_tolerance = 1e-12);

%%
Delta = calculate_defect(problem, guess.x, guess.u, guess.p);
norm(Delta)

%% Test Discretization on Initial Guess
[problem, Delta_disc] = problem.discretize(guess.x, guess.u, guess.p);

x_disc = problem.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = problem.cont_prop(guess.u, guess.p);

%% Solve Problem with PTR
ptr_ops.w_tr = ones(1, Nu) * 5e-2;
ptr_sol_vc = ptr(problem, ptr_ops);

%%
% ptr_ops.w_vse = 1e4;
% ptr_ops.w_tr = 5e-2;
% ptr_ops.w_prime = 1e2;
% ptr_sol_vs = ptr_virtual_state(problem, ptr_ops, "CVX");

%%
ptr_sol = ptr_sol_vc;

if ~ptr_sol.converged
    ptr_sol.converged_i = ptr_ops.iter_max;
end

x = ptr_sol.x(:, :, ptr_sol.converged_i + 1);
u = ptr_sol.u(:, :, ptr_sol.converged_i + 1);
p = ptr_sol.p(:, ptr_sol.converged_i + 1);

t_k_scaled = t_k * ptr_sol.p(1, ptr_sol.converged_i + 1);

i = ptr_sol.converged_i + 1;
[t_cont_sol, x_cont_sol, u_cont_sol] = problem.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));

earth_trajectory = earth_pos(t_cont_sol' * p(1));
mars_trajectory = mars_pos(t_cont_sol' * p(1));

figure
plot(x_cont_sol(1, :), x_cont_sol(2, :), DisplayName="Minimum Fuel Transfer", LineWidth=1); hold on
quiver(x(1, 1:Nu)',x(2, 1:Nu)',u(1, :)',u(2, :)', 1, DisplayName="Thrust"); hold on
plot(earth_trajectory(1, :), earth_trajectory(2, :), DisplayName="Earth Orbit"); hold on;
plot(mars_trajectory(1, :), mars_trajectory(2, :), DisplayName="Mars Orbit"); hold on
scatter(x_0(1), x_0(2), DisplayName="Earth at t = 0"); hold on
scatter(mars_trajectory(1, end), mars_trajectory(2, end), DisplayName=sprintf("Mars at t = %.4g", p(1))); hold off;
title("Minimum Fuel Transfer Between Earth and Mars")
subtitle(sprintf("For final time of %.4g", p(1)))
xlabel("X [AU]")
ylabel("Y [AU]")
legend(Location="eastoutside")
grid on
axis equal

%% Plot optimal Control
fuel = sum(squeeze(vecnorm(u(:, 1:end))) .* delta_t, 2);

figure
tiledlayout(1, 1, "TileSpacing","compact")

nexttile
if u_hold == "ZOH"
    stairs(t_k_scaled(1:Nu), u');
elseif u_hold == "FOH"
    plot(t_k_scaled(1:Nu), u);
end
title("Optimal Control History", Interpreter="latex")
subtitle(sprintf("Total Fuel Used of %.3f units", fuel(end)))
xlabel("Time []")
ylabel("Control Value")
xlim([0, tf])
legend("u_1", "u_2", location = "south", Orientation="horizontal")
grid on

%%

function [vvec] = v_circ(rvec, nu, mu)
    r = sqrt(rvec(1) ^ 2 + rvec(2) ^ 2);
    v = sqrt(mu / r);
    vvec = v * [-sin(nu); cos(nu)];
end