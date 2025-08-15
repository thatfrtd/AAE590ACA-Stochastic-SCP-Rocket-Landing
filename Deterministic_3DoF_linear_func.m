function [X_6DoF, U_6DoF] = Deterministic_3DoF_linear_func(x_0, tf, N, T_max, T_min, alpha, L, glideslope_angle_max, u_hold, g)

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
% Problem Parameters
delta_t = tf / (N - 1); % [s]
r_0 = x_0(1:3); % [km]
v_0 = x_0(4:6); % [km / s]
m_0 = x_0(end); % [kg]

vehicle = Vehicle(m_0 - 600, L, L * 3, 0, T_min, T_max, alpha = alpha);

x_0 = [r_0; v_0; log(m_0)];
x_f = zeros(6, 1);

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);

Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

nx = 7;
nu = 4;

% Algorithm Parameters
default_tolerance = 1e-12;
tolerances = odeset(RelTol=default_tolerance, AbsTol=default_tolerance);

% PTR algorithm parameters
ptr_ops.iter_max = 20;
ptr_ops.iter_min = 1;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 1e2;
ptr_ops.w_tr = ones(1, Nu) * 5e-4;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 3e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = false;

%% Get Dynamics
f = @(t, x, u, p) SymDynamics3DoF_linear(t, x, u, 1, alpha, g);

%% Specify Constraints
z_lb = @(t) log(m_0 - alpha * T_max * t);
z_lb_k = z_lb(t_k);

% Convex state path constraints
glideslope_constraint = {1:N, @(t, x, u, p) norm(x(1:3)) - x(3) / cos(glideslope_angle_max)};
min_mass_constraint = @(t, x, u, p) z_lb(t) - x(7);
state_convex_constraints = {glideslope_constraint};

% Convex control constraints
max_thrust_constraint = {1:N, @(t, x, u, p) u(4) - T_max * exp(-z_lb(t)) * (1 - (x(7) - z_lb(t)))};
min_thrust_constraint = {1:N, @(t, x, u, p) T_min * exp(-z_lb(t)) * (1 - (x(7) - z_lb(t)) + 0.5 * (x(7) - z_lb(t)) ^ 2) - u(4)};
%min_thrust_constraint = {1:N, @(t, x, u, p) T_min * exp(-x(7)) - u(4)};
lcvx_thrust_constraint = {1:N, @(t, x, u, p) norm(u(1:3))- u(4)}; 
final_thrust_constraint = {N, @(t, x, u, p) norm(u(1:3)) - u(1) / cosd(5)};
control_convex_constraints = {min_thrust_constraint,max_thrust_constraint,lcvx_thrust_constraint,final_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Terminal boundary condition
terminal_bc = @(x, p, x_ref, u_ref) [x(1:6) - x_f; 0];

%% Specify Objective
min_fuel_angular_velocity_objective = @(x, u, p) sum(u(3, :) / T_max + x(6, 1:Nu) .^ 2) * delta_t;
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p) sum(u(4, :)) * delta_t;
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p) sum((u(4, 1:(end - 1)) + u(4, 2:end)) / 2) * delta_t;
end

%% Create Guess
% Guess doesn't matter since problem is totally convex and the dynamics are
% linear
sl_guess = guess_3DoF([x_0(1:6); zeros([6, 1])], [x_f; zeros([6, 1])], N, Nu, delta_t, vehicle);
sl_guess.x = sl_guess.x(1:6, :);
if u_hold == "ZOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [cumsum(sl_guess.u(3, :) * delta_t), sum(sl_guess.u(3, :)) * delta_t]];
elseif u_hold == "FOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * cumsum(sl_guess.u(3, :) * delta_t)];
end
sl_guess.u = sl_guess.u ./ sl_guess.x(7, 1:Nu);
sl_guess.u(3:4, :) = [zeros([1, Nu]); sl_guess.u(3, :)];
sl_guess.x(7, :) = log(sl_guess.x(7, :));

guess = sl_guess;

%% Construct Problem Object
prob_3DoF = DeterministicProblem(x_0, x_f, N, u_hold, tspan(end), f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc);

%% Test Discretization
[prob_3DoF, Delta_disc] = prob_3DoF.discretize(guess.x, guess.u, guess.p);

%% Check with Matrix Exponential
A_k_exp = expm((t_k(2) - t_k(1)) * prob_3DoF.cont.A(0, x_0, guess.u(:, 1), 0));
A_k_ck = sum(pagenorm(prob_3DoF.disc.A_k(:, :, 1:(end - 1)) - A_k_exp), "all") < default_tolerance; % Checks out

%% Solve Problem with PTR
ptr_sol = ptr(prob_3DoF, ptr_ops);

X = ptr_sol.x(:, :, ptr_sol.converged_i);
U = ptr_sol.u(:, :, ptr_sol.converged_i);

%% Transform Solution to 6DoF
% 3DoF state: [r; v; m]
% 6DoF state: [r; v; q; w; m]

% Assume gimbal is always 0 so thrust direction is where vehicle points
x_B = U(1:3, :) ./ vecnorm(U(1:3, :));
cross_prod = cross(x_B, [ones([1, Nu]); zeros([2, Nu])]);
q_angle = asind(vecnorm(cross_prod));
q_dir = cross_prod ./ vecnorm(cross_prod);
q = [q_dir .* sind(q_angle / 2); cosd(q_angle / 2)];

w = zeros([3, N]);
X_6DoF = [X(1:6, :); q; w; exp(X(7, :))];

% 3DoF control: [a; |a|]
% 6DoF control: [T; gamma]
thrust = [vecnorm(U(1:3, :)) .* exp(X(7, 1:Nu)); zeros([1, Nu]); zeros([1, Nu])];
thrust_ck = quat_rot_array(q_conj(q), thrust);

gamma = zeros([1, Nu]);
U_6DoF = [thrust; gamma];

%% Plot Solution
[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(ptr_sol.u(:, :, ptr_sol.converged_i), ptr_sol.p(:, ptr_sol.converged_i));

t_scaled = t_cont_sol;

tiledlayout(1, 4)
nexttile
plot(t_scaled, x_cont_sol([1,3], :)) % - also include continuous solution and look at error?
title("Position History")
xlabel("Time [s]")
ylabel("Position [km]")
legend("r_x", "r_y", Location="southoutside", Orientation="horizontal")
grid on

nexttile
plot(t_scaled, x_cont_sol([4, 6], :) * 1000) % - also include continuous solution and look at error?
title("Velocity History")
xlabel("Time [s]")
ylabel("Velocity [m / s]")
legend("v_x", "v_y", Location="southoutside", Orientation="horizontal")
grid on

nexttile
plot(t_scaled, exp(x_cont_sol(7, :))) % - also include continuous solution and look at error?
title("Mass History")
xlabel("Time [s]")
ylabel("Mass [kg]")
grid on

nexttile
if u_hold == "ZOH"
    stairs(t_scaled(1:size(u_cont_sol, 2)), (u_cont_sol(:, :) .* exp(x_cont_sol(end, 1:size(u_cont_sol, 2))))')
elseif u_hold == "FOH"
    plot(t_scaled(1:size(u_cont_sol, 2)), (u_cont_sol(:, :) .* exp(x_cont_sol(end, 1:size(u_cont_sol, 2))))')
end
title("Control History")
xlabel("Time [s]")
ylabel("Thrust [kN]")
legend("T_x", "T_y", "\sigma", Location="southoutside", Orientation="horizontal")
grid on

sgtitle("State and Control Histories for Mars Optimal Fuel Rocket Landing")

%% Plot Solution 2D
figure
plot3(X(1, :), X(2, :), X(3, :), DisplayName="Trajectory"); hold on
quiver3(X(1, 1:Nu), X(2, 1:Nu), X(3, 1:Nu), U(1, :), U(2, :), U(3, :), DisplayName = "Thrust")
grid on
title("3D Plot of Mars Optimal Fuel Rocket Landing")
xlabel("r_1 [km]")
ylabel("r_2 [km]")
ylabel("r_3 [km]")
legend(Location="southoutside", Orientation="horizontal")
axis equal

%% Quantify Lossiness from Convexification
figure
lcvx_err = abs(vecnorm(U(1:3, :)) - U(4, :));
plot(1:Nu, lcvx_err);
yscale("log")
grid on
xlabel("Time [s]")
ylabel("Error")
title("Lossless Convexification Error")

end