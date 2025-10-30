%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% HW4 Q1c
% Author: Travis Hastreiter 
% Created On: 26 February, 2025
% Description: Single shooting for minimum fuel transfer from Earth to
% Mars with fixed final time and max thrust constraint
% Most Recent Change: 26 February, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assumptions for minimum fuel transfer
% - Earth and Mars have circular orbits
% - No SRP or other perturbation except thrust
% - No gravity from Earth, Mars
% - Constant mass
% - Fixed final time

%% Initialize
a_earth = 1; % [AU]
nu0_earth = 0; % [rad]

a_mars = 1.524; % [AU]
nu0_mars = pi; % [rad]

% All times are nondimensionalized by Earth's mean motion
tf = 8; % [] arrival time (nondimensionalized)

mu = 1; % [] Sun's gravitaional parameter (nondimensionalized)

lambda0_guess = [0.8; 0.1; 0.2; 1.1]; % initial costate guess for rho = 1

u_max = 0.1;

B = [zeros(2); eye(2)];

default_tolerance = 1e-12;

rho_array = 10 .^ [0; -1; -2; -3];


%% c.3)
% Solution
P_earth_over_P_mars = (a_earth / a_mars) ^ (3/2);

earth_pos = @(t) a_earth .* [cos(t + nu0_earth); ...
                             sin(t + nu0_earth)];
mars_pos = @(t) a_mars .* [cos(t * P_earth_over_P_mars + nu0_mars); ...
                           sin(t * P_earth_over_P_mars + nu0_mars)];

nuf_mars = tf * P_earth_over_P_mars + nu0_mars;

x0 = [earth_pos(0); v_circ(earth_pos(0), nu0_earth, mu)];
xf = [mars_pos(tf); v_circ(mars_pos(tf), nuf_mars, mu)];

initial_violation = test_lambda0(lambda0_guess, x0, tf, xf, mu, B, u_max, rho_array(1), default_tolerance);

%% Solve for initial costate
options = optimoptions("fsolve", 'OptimalityTolerance', 1e-8);

lambda0_array = zeros(4, numel(rho_array));
n = 1000;
t = linspace(0, tf, n);
tf_array = t(2:end);
x_array = zeros(4, n, numel(rho_array));
lambda_array = zeros(4, n, numel(rho_array));
u_array = zeros(2, n, numel(rho_array));
ls = zeros(4, numel(rho_array));

for r = 1:numel(rho_array)
    rho = rho_array(r);
    lambda0 = fsolve(@(lambda0) test_lambda0(lambda0, x0, tf_array, xf, mu, B, u_max, rho, default_tolerance), lambda0_guess, options);
    lambda0_guess = lambda0;
    lambda0_array(:, r) = lambda0;
    ls(:,r) = lambda0;

    [x_sol, lambda, ~] = propagate_initial_conditions(x0, lambda0, tf_array, mu, B, u_max, rho, default_tolerance);
    u_sol = optimal_control(lambda', B, u_max, rho);
    x_array(:, :, r) = x_sol';
    lambda_array(:, :, r) = lambda';
    u_array(:, :, r) = u_sol;
end

x = x_array(:, :, end);
u = u_array(:, :, end);

%% Test solution
earth_trajectory = earth_pos(t);
mars_trajectory = mars_pos(t);

figure
plot(x(1, :)', x(2, :)', DisplayName="Minimum Fuel Transfer"); hold on
quiver(x(1, 1:30:end)',x(2, 1:30:end)',u(1, 1:30:end)',u(2, 1:30:end)', 1, DisplayName="Thrust"); hold on
plot(earth_trajectory(1, :), earth_trajectory(2, :), DisplayName="Earth Orbit"); hold on;
plot(mars_trajectory(1, :), mars_trajectory(2, :), DisplayName="Mars Orbit"); hold on
scatter(x0(1), x0(2), DisplayName="Earth at t = 0"); hold on
scatter(xf(1), xf(2), DisplayName=sprintf("Mars at t = %.g", tf)); hold off;
title("Minimum Fuel Transfer Between Earth and Mars")
subtitle(sprintf("For final time of %.g", tf))
xlabel("X [AU]")
ylabel("Y [AU]")
legend(Location="eastoutside")
grid on
axis equal

%% Plot optimal Control
fuel = sum(squeeze(vecnorm(u_array(:, 2:end, :))) .* diff(t)', 1);

figure
tiledlayout(1, 2, "TileSpacing","compact")

nexttile
plot(t, u)
title("Optimal Control History for $\rho = 10^{-3}$", Interpreter="latex")
subtitle(sprintf("Total Fuel Used of %.3f units", fuel(end)))
xlabel("Time []")
ylabel("Control Value")
xlim([0, tf])
legend("u_1", "u_2", location = "south", Orientation="horizontal")
grid on

nexttile
plot(t, lambda)
title("Optimal Costate History for $\rho = 10^{-3}$", Interpreter="latex")
xlabel("Time []")
ylabel("Costate Value")
xlim([0, tf])
legend("\lambda" + string(1:4), location = "south", Orientation="horizontal")
grid on

sgtitle("Optimal Control and Costate Histories for Minimum Fuel Transfer")

figure
tiledlayout(1, 2, "TileSpacing","compact")

nexttile
plot(t, squeeze(u_array(1, :, :)))
title("u_1 vs Time")
xlabel("Time []")
ylabel("Control Magnitude")
legend("$\rho = " + string(rho_array) + "$", location = "southoutside", Orientation="horizontal", Interpreter="latex")
grid on

nexttile
plot(t, squeeze(u_array(2, :, :)))
title("u_2 vs Time")
xlabel("Time []")
ylabel("Control Magnitude")
legend("$\rho = " + string(rho_array) + "$", location = "southoutside", Orientation="horizontal", Interpreter="latex")
grid on

sgtitle("Optimal Control vs Time for Different Smoothing Parameters")

figure
bar("\rho = " + string(rho_array), fuel, Interpreter="latex")
xlabel("Smoothing Parameter")
ylabel("Total Fuel")
title("Total Fuel Used vs Smoothing Parameter")
grid on

%% Plot Control Hamiltonian
H = zeros(size(t));

for i = 1:numel(t)
    H(i) = control_hamiltonian(x_array(:, i, end), lambda_array(:, i, end), mu, B, u_max, rho_array(end));
end

figure
plot(t, abs(H - H(1)))
title("Control Hamiltonian Change vs Time")
subtitle("Minimum Fuel Transfer")
xlabel("Time []")
ylabel("|H(t) - H_0|")
yscale("log")
grid on;

%% Helper Functions
function [xdot] = state_equation(x, lambda, mu, B, u_max, rho)
    rvec = x(1:2);
    vvec = x(3:4);
    r = norm(rvec);
    
    f0 = [vvec; -mu / r ^ 3 * rvec];

    u = optimal_control(lambda, B, u_max, rho);
    xdot = f0 + B * u;
end

function [lambdadot] = costate_equation(x, lambda, mu)
    rvec = x(1:2);
    r = norm(rvec);
    rhat = rvec / r;

    lambdadot = -[zeros(2), mu / r ^ 3 * (3 * (rhat * rhat.') - eye(2)); ...
                  eye(2), zeros(2)] * lambda;
end

function [vvec] = v_circ(rvec, nu, mu)
    r = norm(rvec);
    v = sqrt(mu / r);
    vvec = v * [-sin(nu); cos(nu)];
end

function [x, lambda, t] = propagate_initial_conditions(x0, lambda0, tf, mu, B, u_max, rho, error_tolerance)
    tolerances = odeset(RelTol=error_tolerance, AbsTol=error_tolerance);
    
    y0 = [x0; lambda0];
    
    ydot = @(x, lambda) [state_equation(x, lambda, mu, B, u_max, rho); ... 
                         costate_equation(x, lambda, mu)];

    [t, y] = ode45(@(t, y) ydot(y(1:4), y(5:8)), [0, tf], y0, tolerances);

    x = y(:, 1:4);
    lambda = y(:, 5:8);
end

function [u] = optimal_control(lambda, B, u_max, rho)
    p = -pagemtimes(B.', lambda);
    S = vecnorm(p) - 1;
    Gamma_tilde = (u_max / 2) .* (1 + tanh(S / rho));
    u = p ./ vecnorm(p) .* Gamma_tilde;
    u(isnan(u)) = 0;
end

function [H] = control_hamiltonian(x, lambda, mu, B, u_max, rho)
    u = optimal_control(lambda, B, u_max, rho);
    H = norm(u) + lambda.' * state_equation(x, lambda, mu, B, u_max, rho);
end

function [bc_violation] = test_lambda0(lambda0, x0, tf, xf_tar, mu, B, u_max, rho, error_tolerance)
    [x, lambda, ~] = propagate_initial_conditions(x0, lambda0, tf, mu, B, u_max, rho, error_tolerance);

    xf = x(end, :)';
    lambdaf = lambda(end, :)';

    bc_violation = xf - xf_tar;
end