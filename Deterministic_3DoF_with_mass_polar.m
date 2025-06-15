%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Vehicle Parameters
alpha = 0.5086; % [s / km]
m_dry = 1500; % [kg]
m_wet = 600; % [kg]
m_0 = m_dry + m_wet;
T_max = 3 * m_0 * 9.81e-3; % [kg km / s2]
T_min = 0.55 * T_max; % [kg km / s2]
I = 150000 * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 3e-3; % [km] Distance from CoM to nozzle
gimbal_max = deg2rad(6); % [rad]
 
vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I, alpha = alpha);
%vehicle_big_gimbs = Vehicle(m_dry, L, L * 3, deg2rad(20), T_min, T_max, I = I, alpha = alpha);

% Problem Parameters
tf = 35; % [s]
N = 15; % []
r_0 = [0; 4.6]; % [km]
theta_0 = deg2rad(120); % [rad]
v_0 = make_R2(-deg2rad(60)) * [0.306; 0]; % [km / s]
w_0 = deg2rad(0); % [rad / s]
glideslope_angle_max = deg2rad(70); % [rad]

x_0 = [r_0; v_0; theta_0; w_0; m_0];
x_f = [zeros(2, 1); zeros(2, 1); pi / 2; 0];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "FOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

initial_guess = "straight line"; % "CasADi" or "straight line"

% PTR algorithm parameters
ptr_ops.iter_max = 50;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 1e5;
ptr_ops.w_tr = ones(1, Nu) * 5e-2;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 2e-2;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = false;

%% Get Dynamics
f = @(t, x, u, p) SymDynamics3DoF_mass_polar(t, x, u, vehicle.L, vehicle.I(2), vehicle.alpha);

%% Specify Constraints
% Convex state path constraints
glideslope_constraint = @(t, x, u, p) norm(x(1:2)) - x(2) / cos(glideslope_angle_max);
state_convex_constraints = {glideslope_constraint};

% Convex control constraints
max_thrust_constraint = @(t, x, u, p) u(1) - T_max;
min_thrust_constraint = @(t, x, u, p) T_min - u(1);
max_gimbal_constraint = @(t, x, u, p) abs(u(2)) - gimbal_max;
control_convex_constraints = {min_thrust_constraint,max_gimbal_constraint,max_thrust_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Terminal boundary conditions
terminal_bc = @(x, u, p) [x(1:6, :) - x_f; 0];

%% Specify Objective
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p) -x(7, end);% sum(u(1, :)) * delta_t;
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p) sum((u(1, 1:(end - 1)) + u(1, 2:end)) / 2) * delta_t; %-x(7, end);
end

%% Create Guess
sl_guess = guess_3DoF(x_0(1:6), x_f + [0; 0; 0; 0; 0; 0], N, Nu, delta_t, vehicle);
if u_hold == "ZOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [cumsum(sl_guess.u(3, :) * delta_t), sum(sl_guess.u(3, :)) * delta_t]];
elseif u_hold == "FOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * cumsum(sl_guess.u(3, :) * delta_t)]
end

CasADi_sol = CasADi_solve_mass(x_0, sl_guess.x, sl_guess.u, vehicle, N, delta_t, glideslope_angle_max);

if initial_guess == "straight line"
    guess = sl_guess;
elseif initial_guess == "CasADi"
    guess = CasADi_sol;
end
if u_hold == "ZOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "previous","extrap")';
elseif u_hold == "FOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "linear","extrap")';
end
guess.p = sl_guess.p;

figure
plot_3DoF_trajectory(t_k, sl_guess.x, sl_guess.u, glideslope_angle_max, gimbal_max, T_min, T_max)

figure
plot_3DoF_time_histories(t_k, sl_guess.x, sl_guess.u)

figure
plot_3DoF_trajectory(t_k, guess.x, guess.u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
plot_3DoF_time_histories(t_k, guess.x, guess.u)

guess.u = [guess.u(3, :); atan2(guess.u(2, :), guess.u(1, :))];

%% Construct Problem Object
prob_3DoF = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc, discretization_method = "error", N_sub = 1);
prob_3DoF.vehicle = vehicle;

%% Test Scaling
% guess_scaled.x = prob_3DoF.scale_x(guess.x);
% guess_scaled.u = prob_3DoF.scale_u(guess.u);
% guess_scaled.p = prob_3DoF.scale_p(guess.p);
% 
% figure
% plot_3DoF_trajectory(t_k, guess_scaled.x, guess_scaled.u, glideslope_angle_max, gimbal_max, 0, 0)
% 
% figure
% plot_3DoF_time_histories(t_k, guess_scaled.x, guess_scaled.u)

%%
Delta = calculate_defect(prob_3DoF, guess.x, guess.u, guess.p);
norm(Delta)

%% Test Discretization on Initial Guess

[prob_3DoF, Delta_disc] = prob_3DoF.discretize(guess.x, guess.u, guess.p);

x_disc = prob_3DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_3DoF.cont_prop(guess.u, guess.p);
u_cont = [u_cont; vecnorm(u_cont)];

figure
comparison_plot_3DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%%
% x_sol_p = x_sol;
% x_sol_p(7, :) = exp(x_sol_p(7, :));
% u_sol_p = u_sol(1:2, :);
% u_sol_p(1, :) = u_sol(3, :) .* x_sol_p(7, :);
% u_sol_p(2, :) = atan2(u_sol(2, :), u_sol(1, :));
% 
% guess.x = x_sol_p;
% guess.u = u_sol_p;
% 
% [prob_3DoF, Delta_disc] = prob_3DoF.discretize(x_sol_p, u_sol_p, guess.p);
% 
% x_disc = prob_3DoF.disc_prop(u_sol_p, guess.p);
% 
% [t_cont, x_cont, u_cont] = prob_3DoF.cont_prop(u_sol_p, guess.p);
% 
% figure
% comparison_plot_3DoF_trajectory({x_sol_p, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%%
guess3 = guess;
guess3.u = [guess3.u; vecnorm(guess3.u)];

figure
comparison_plot_3DoF_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess3.u, u_cont, guess3.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%%
nx = 7;
nu = 2;
np = 0;

% Linearized matrices
t_sym = sym("t");
x_sym = sym("x", [nx, 1]);
u_sym = sym("u", [nu, 1]);
p_sym = sym("p", [np, 1]);

A = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
B = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), u_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
S = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), p_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);

A_k = prob_3DoF.disc.A_k;
B_k_plus = prob_3DoF.disc.B_plus_k;
B_k_minus = prob_3DoF.disc.B_minus_k;
S_k = prob_3DoF.disc.E_k;
d_k = prob_3DoF.disc.c_k;

%% Try RK4 Discretization
N_sub = 100;
[A_k_rk4, B_k_plus_rk4, B_k_minus_rk4, S_k_rk4, d_k_rk4, Delta_rk4] = discretize_error_dynamics_FOH_RK4(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk4 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk4 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk4 - B_k_plus), "all");
%S_err = sum(pagenorm(S_k_rk4 - S_k), "all");
d_err = sum(pagenorm(d_k_rk4 - d_k), "all");
Delta_err = norm(Delta_rk4 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.5f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);

%% Try RKV6(5) Discretization
N_sub = 1;
[A_k_rk65, B_k_plus_rk65, B_k_minus_rk65, S_k_rk65, d_k_rk65, Delta_rk65] = discretize_error_dynamics_FOH_RKV65(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk65 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk65 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk65 - B_k_plus), "all");
%S_err = sum(pagenorm(S_k_rk65 - S_k), "all");
d_err = sum(pagenorm(d_k_rk65 - d_k), "all");
Delta_err = norm(Delta_rk65 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.10f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);

%% Try RKV8(7) Discretization
N_sub = 1;
[A_k_rk87, B_k_plus_rk87, B_k_minus_rk87, S_k_rk87, d_k_rk87, Delta_rk87] = discretize_error_dynamics_FOH_RKV87(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk87 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk87 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk87 - B_k_plus), "all");
%S_err = sum(pagenorm(S_k_rk87 - S_k), "all");
d_err = sum(pagenorm(d_k_rk87 - d_k), "all");
Delta_err = norm(Delta_rk87 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.10f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);



%% Solve Problem with PTR
ptr_ops.w_tr = ones(1, Nu) * 5e-2;
ptr_sol_vc = ptr(prob_3DoF, ptr_ops);

%%
ptr_ops.w_vse = 1e6;
ptr_ops.w_tr = 5e-2;
ptr_ops.w_prime = 1e2;
ptr_sol_vs = ptr_virtual_state(prob_3DoF, ptr_ops, "CVX");

%%
ptr_sol = ptr_sol_vc;

%%
figure
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [prob_3DoF.objective(prob_3DoF.guess.x, prob_3DoF.guess.u, prob_3DoF.guess.p), [ptr_sol.info.J]]); hold on
yline(CasADi_sol.objective); hold off
legend("PTR Iterations", "CasADi Solution")
title("Objective vs Iteration")
grid on

nexttile
plot(ptr_sol.delta_xp)
yscale("log")
title("Stopping Criteria vs Iteration")
grid on

nexttile
plot(0:ptr_sol.converged_i, vecnorm(ptr_sol.Delta(:, 1:(ptr_sol.converged_i + 1)), 2, 1))
yscale("log")
title("Defect Norm vs Iteration")
grid on

%%
i = 12;% ptr_sol.converged_i;

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));

u_sol_p = ptr_sol.u(:, :, i);
u_sol_p = [u_sol_p(1, :) .* cos(u_sol_p(2, :)); u_sol_p(1, :) .* sin(u_sol_p(2, :)); u_sol_p(1, :)];

u_cont_sol_p = ptr_sol.u(:, :, i);
u_cont_sol_p = [u_cont_sol_p(1, :) .* cos(u_cont_sol_p(2, :)); u_cont_sol_p(1, :) .* sin(u_cont_sol_p(2, :)); u_cont_sol_p(1, :)];

figure
plot_3DoF_trajectory(t_k, ptr_sol.x(:, :, i), u_sol_p, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1)

figure
comparison_plot_3DoF_trajectory({guess.x, x_cont_sol, ptr_sol.x(:, :, i), CasADi_sol.x}, ["Guess", "Continuous Propagation", "Solution Output", "CasADi"], glideslope_angle_max, linestyle = [":", "-", "--", "-"], title = "Continuous vs Discrete Propagation of Solution")

figure
comparison_plot_3DoF_time_histories({t_k, t_cont_sol, t_k}, {guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, {guess3.u, u_cont_sol_p, u_sol_p}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")
