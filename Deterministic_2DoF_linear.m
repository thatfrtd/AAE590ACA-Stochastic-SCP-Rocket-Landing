%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 19 April, 2025
% Description: 2DoF (all translational) landing of rocket using PTR SCP 
% algorithm
% Most Recent Change: 19 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Physical Parameters
g = [0; -3.7114e-3]; % [km / s2]
alpha = 0.5086; % [s / km]
T_min = 4.97; % [kg km / s2]
T_max = 13.26; % [kg km / s2]

% Problem Parameters
delta_t = 1.5; % [s]
N = 40; % []
r_0 = [1; 1.5]; % [km]
v_0 = [0; 0]; % [km / s]
m_0 = 2000; % [kg]
gamma_min = deg2rad(3); % [rad]

L = 3e-3;
vehicle = Vehicle(m_0 - 600, L, L * 3, 0, T_min, T_max, alpha = alpha);

x_0 = [r_0; v_0; log(m_0)];
x_f = zeros(4, 1);

tspan = [0, N * delta_t];
t_k = linspace(tspan(1), tspan(2), N);

u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

nx = 5;
nu = 3;

% Algorithm Parameters
default_tolerance = 1e-12;
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

% PTR algorithm parameters
ptr_ops.iter_max = 20;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 1e5;
ptr_ops.w_tr = ones(1, Nu) * 5e0;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

%% Get Dynamics
f = @(t, x, u, p) SymDynamics2DoF_linear(t, x, u, alpha);

%% Specify Constraints
z_lb = @(t) log(m_0 - alpha * T_max * t);
z_lb_k = z_lb(t_k);

% Convex state path constraints
glideslope_constraint = @(t, x, u, p) norm(x(1:2)) - x(2) / cos(pi/2 - gamma_min);
min_mass_constraint = @(t, x, u, p) z_lb(t) - x(5);
state_convex_constraints = {glideslope_constraint};

% Convex control constraints
max_thrust_constraint = @(t, x, u, p) u(3) - T_max * exp(-z_lb(t)) * (1 - (x(5) - z_lb(t)));
min_thrust_constraint = @(t, x, u, p) T_min * exp(-x(5)) - u(3);
lcvx_thrust_constraint = @(t, x, u, p) norm(u(1:2))- u(3); 
control_convex_constraints = {min_thrust_constraint,max_thrust_constraint,lcvx_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Terminal boundary condition
terminal_bc = @(x, p) [x(1:4) - x_f; 0];

%% Specify Objective
min_fuel_angular_velocity_objective = @(x, u, p) sum(u(3, :) / T_max + x(6, 1:Nu) .^ 2) * delta_t;
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p) sum(u(3, :)) * delta_t;
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p) sum((u(3, 1:(end - 1)) + u(3, 2:end)) / 2) * delta_t;
end

%% Create Guess
% Guess doesn't matter since problem is totally convex and the dynamics are
% linear
sl_guess = guess_3DoF([x_0(1:4); 0; 0], [x_f; 0; 0] + [0; 0; 0; 0; 0; 0], N, Nu, delta_t, vehicle);
sl_guess.x = sl_guess.x(1:4, :);
if u_hold == "ZOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [cumsum(sl_guess.u(3, :) * delta_t), sum(sl_guess.u(3, :)) * delta_t]];
elseif u_hold == "FOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * cumsum(sl_guess.u(3, :) * delta_t)];
end

guess = sl_guess;

%% Construct Problem Object
prob_2DoF = DeterministicProblem(x_0, x_f, N, u_hold, tspan(end), f, guess, convex_constraints, min_fuel_objective, scale = false, terminal_bc = terminal_bc);

%% Test Discretization
[prob_2DoF, Delta_disc] = prob_2DoF.discretize(guess.x, guess.u, guess.p);

%% Check with Matrix Exponential
A_k_exp = expm((t_k(2) - t_k(1)) * prob_2DoF.cont.A(0, x_0, guess.u(:, 1), 0));
A_k_ck = sum(pagenorm(prob_2DoF.disc.A_k(:, :, 1:Nu) - A_k_exp), "all") < default_tolerance; % Checks out

%% Solve Problem with PTR
ptr_sol = ptr(prob_2DoF, ptr_ops);

X = ptr_sol.x(:, :, ptr_sol.converged_i);
U = ptr_sol.u(:, :, ptr_sol.converged_i);

%% Plot Solution
tiledlayout(1, 4)
nexttile
plot(t_k, X(1:2, :)) % - also include continuous solution and look at error?
title("Position History")
xlabel("Time [s]")
ylabel("Position [km]")
legend("r_x", "r_y", Location="southoutside", Orientation="horizontal")
grid on

nexttile
plot(t_k, X(3:4, :) * 1000) % - also include continuous solution and look at error?
title("Velocity History")
xlabel("Time [s]")
ylabel("Velocity [m / s]")
legend("v_x", "v_y", Location="southoutside", Orientation="horizontal")
grid on

nexttile
plot(t_k, exp(X(5, :))) % - also include continuous solution and look at error?
title("Mass History")
xlabel("Time [s]")
ylabel("Mass [kg]")
grid on

nexttile
stairs(t_k(2:end), (U(:, :) .* exp(X(end, 2:end)))')
title("Control History")
xlabel("Time [s]")
ylabel("Thrust [kN]")
legend("T_x", "T_y", "\sigma", Location="southoutside", Orientation="horizontal")
grid on

sgtitle("State and Control Histories for Mars Optimal Fuel Rocket Landing")

%% Plot Solution 2D
figure
plot(X(1, :), X(2, :), DisplayName="Trajectory"); hold on
quiver(X(1, 2:end), X(2, 2:end), U(1, :), U(2, :), DisplayName = "Thrust")
grid on
title("2D Plot of Mars Optimal Fuel Rocket Landing")
xlabel("r_1 [km]")
ylabel("r_2 [km]")
legend(Location="southoutside", Orientation="horizontal")
axis equal

%% Quantify Lossiness from Convexification
figure
lcvx_err = abs(vecnorm(U(1:2, :)) - U(3, :));
plot(1:Nu, lcvx_err);
yscale("log")
grid on
xlabel("Time [s]")
ylabel("Error")
title("Lossless Convexification Error")