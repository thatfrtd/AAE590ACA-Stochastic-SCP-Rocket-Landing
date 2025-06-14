%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA
% AC Optimal Control
% Author: Travis Hastreiter 
% Created On: 14 May, 2025
% Description: Test of parts of customized PIPG conic convex solver
% Most Recent Change: 17 May, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO:
% - Work
% - See affect of constraining xi(0) = x(0) - probably better for mpc

%% Test with Planar 3DoF Rocket Landing
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
glideslope_angle_max = deg2rad(50); % [rad]

x_0 = [r_0; v_0; theta_0; w_0; m_0];
x_f = [zeros(2, 1); zeros(2, 1); pi / 2; 0];

u_0 = [T_min; 0];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "FOH";
Nu = N;

initial_guess = "straight line"; % "CasADi" or "straight line"

integration_tolerance = 1e-12;
tolerances = odeset(RelTol=integration_tolerance, AbsTol=integration_tolerance);

% PTR algorithm parameters
ptr_ops.iter_max = 20;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 1e5;
ptr_ops.w_tr = ones(1, Nu) * 5e-2;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = true;

nx = 7;
nu = 2;
np = 1;

%% Get Dynamics
f = @(t, x, u, p) SymDynamics3DoF_mass_polar(t, x, u, vehicle.L, vehicle.I(2), vehicle.alpha);

% Linearized matrices
t_sym = sym("t");
x_sym = sym("x", [nx, 1]);
u_sym = sym("u", [nu, 1]);
p_sym = sym("p", [np, 1]);

A = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
B = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), u_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
S = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), p_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);

% Initial guess
sl_guess = guess_3DoF(x_0(1:6), x_f, N, Nu, delta_t, vehicle);
if u_hold == "ZOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * [cumsum(sl_guess.u(3, :) * delta_t), sum(sl_guess.u(3, :)) * delta_t]];
elseif u_hold == "FOH"
    sl_guess.x = [sl_guess.x; m_0 - alpha * cumsum(sl_guess.u(3, :) * delta_t)];
end

guess = sl_guess;
guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "linear","extrap")';
guess.u = guess.u(1:2, :);
guess.u(:, 1) = u_0;
guess.p = 25;

% guess.x = x;
% guess.u = u_p;
% guess.p = s;

% Constraints
% State constraints
D_xi = SOC(1:2, tan(ones([1, N - 2]) * glideslope_angle_max)); % Glideslope
% Control constraints
D_u = [Box(1, ones([1, N - 1]) * T_min, ones([1, N - 1]) * T_max), ... % Min and max thrust
       Box(2, ones([1, N - 1]) * -gimbal_max, ones([1, N - 1]) * gimbal_max)]; % Gimbal 
% Parameter constraints
D_s = Box(1, 5, 100);

% Initial and final conditions
D_initial_x = Singleton(1:7, x_0);
D_final_x = Singleton(1:6, x_f);

D_initial_u = Singleton(1:2, u_0);

D_comb = [D_xi, D_u, D_s];
D_struct.xi = D_xi;
D_struct.u = D_u;
D_struct.s = D_s;
D_struct.ic_xi = D_initial_x;
D_struct.ic_u = D_initial_u;
D_struct.tc_xi = D_final_x;

convex_constraints = D_comb.get_constraint_functions;

% Discretize
terminal_bc = @(x, u, p) [x(1:6, :) - x_f; 0];
min_fuel_objective = @(x, u, p) sum((u(3, 1:(end - 1)) + u(3, 2:end)) / 2) * delta_t;
prob_3DoF = DeterministicProblem(x_0, x_f, N, u_hold, tf, f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc, discretization_method = "error", N_sub = 1);
[prob_3DoF, Delta_disc] = prob_3DoF.discretize(guess.x, guess.u, guess.p);

A_k = prob_3DoF.disc.A_k;
B_k_plus = prob_3DoF.disc.B_plus_k;
B_k_minus = prob_3DoF.disc.B_minus_k;
S_k = prob_3DoF.disc.E_k;
d_k = prob_3DoF.disc.c_k;

% Penalty weights
w_m = 1; % total fuel used penalty
w_vse = 1e6; % virtual state error
w_tr = 1e-2; % trust region
w_tr_s = 1; % parameter trust region

% qs
q_oc = [];
q_oc.x = [zeros([nx, N - 1]), [zeros([nx - 1, 1]); -w_m]];
q_oc.xi = zeros([nx, N]);
q_oc.u = zeros([nu, N]);
q_oc.s = zeros([np, 1]);

% Algorithm parameters
rho = 1.9;
tol_abs = 1e-8;
tol_rel = 1e-8;
tol_infeas = 1e-3;

j_check = 10;
j_max = 10000;

%% Try RK4 Discretization
N_sub = 1;
[A_k_rk4, B_k_plus_rk4, B_k_minus_rk4, S_k_rk4, d_k_rk4, Delta_rk4] = discretize_error_dynamics_FOH_RK4(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk4 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk4 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk4 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk4 - S_k), "all");
d_err = sum(pagenorm(d_k_rk4 - d_k), "all");
Delta_err = norm(Delta_rk4 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.5f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);

%% Try RKV6(5) Discretization
N_sub = 1;
[A_k_rk65, B_k_plus_rk65, B_k_minus_rk65, S_k_rk65, d_k_rk65, Delta_rk65] = discretize_error_dynamics_FOH_RKV65(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk65 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk65 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk65 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk65 - S_k), "all");
d_err = sum(pagenorm(d_k_rk65 - d_k), "all");
Delta_err = norm(Delta_rk65 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.10f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);

%% Try MEX File RKV6(5) Discretization
N_sub = 1;
[A_k_rk65, B_k_plus_rk65, B_k_minus_rk65, S_k_rk65, d_k_rk65, Delta_rk65] = discretize_error_dynamics_FOH_RKV65_3DoF_mex(N, tspan, guess.x, guess.u, guess.p, N_sub, L, I, alpha);

A_err = sum(pagenorm(A_k_rk65 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk65 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk65 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk65 - S_k), "all");
d_err = sum(pagenorm(d_k_rk65 - d_k), "all");
Delta_err = norm(Delta_rk65 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.10f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);


%% Try RKV8(7) Discretization
N_sub = 1;
[A_k_rk87, B_k_plus_rk87, B_k_minus_rk87, S_k_rk87, d_k_rk87, Delta_rk87] = discretize_error_dynamics_FOH_RKV87(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk87 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk87 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk87 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk87 - S_k), "all");
d_err = sum(pagenorm(d_k_rk87 - d_k), "all");
Delta_err = norm(Delta_rk87 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.10f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);

%% Test Discretization on Initial Guess
x_disc = prob_3DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_3DoF.cont_prop(guess.u, guess.p);
u_cont = [u_cont; vecnorm(u_cont)];

figure
comparison_plot_3DoF_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

guess3 = guess;
guess3.u = [guess3.u; vecnorm(guess3.u)];

figure
plot_3DoF_trajectory(t_k, guess.x, guess3.u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);

%% Construct General Problem
z_guess = [guess.x(:); guess.x(:); guess.u(:); guess.p(:)];

W_state = [w_tr + w_vse, -w_vse; -w_vse, w_vse];
Q_state = kron(W_state, eye(nx * N));
Q_u = w_tr * eye(nu * N);
Q_s = w_tr_s * eye(np);
Q = blkdiag(Q_state, Q_u, Q_s);

z_q_ref = z_guess .* [ones([nx * N, 1]); zeros([nx * N, 1]); ones([nu * N, 1]); 1];
z_q_ref_scaled = z_guess .* [ones([nx * N, 1]) * w_tr; zeros([nx * N, 1]); ones([nu * N, 1]) * w_tr; w_tr_s];
q = [q_oc.x(:); q_oc.xi(:); q_oc.u(:); q_oc.s(:)] - z_q_ref_scaled;
c = z_q_ref' * z_q_ref_scaled / 2;

H_x = zeros([(N - 1) * nx, N * nx]);
H_xi = zeros([(N - 1) * nx, N * nx]);
H_u = zeros([(N - 1) * nx, N * nu]);
H_s = zeros([(N - 1) * nx, np]);
for k = 1 : (N - 1)
    H_x((1 : nx) + nx * (k - 1), (1 : nx) + nx * (k - 1)) = A_k(:, :, k);
    H_x((1 : nx) + nx * (k - 1), (1 : nx) + nx * k) = -1 * eye(nx);
    H_u((1 : nx) + nx * (k - 1), (1 : nu) + nu * (k - 1)) = B_k_minus(:, :, k);
    H_u((1 : nx) + nx * (k - 1), (1 : nu) + nu * k) = B_k_plus(:, :, k);
    H_s((1 : nx) + nx * (k - 1), 1 : np) = S_k(:, :, k);
end
H = [H_x, H_xi, H_u, H_s];

h = -d_k(:);

D = [D_xi.shift_all_indices(reshape((N * nx + nx + 1) : (2 * N * nx - nx), nx, N - 2)), ...
     D_u.shift_all_indices(reshape((2 * N * nx + nu + 1) : (2 * N * nx + N * nu), nu, N - 1)), ...
     D_s.shift_all_indices(2 * N * nx + N * nu + 1)];

D_eq = [...%D_initial_x.shift_all_indices((1 : nx)'), ...
     D_initial_x.shift_all_indices(((N * nx + 1) : (N * nx + nx))'), ...
     D_final_x.shift_all_indices(reshape(((N - 1) * nx + N * nx + 1) : (2 * N * nx), nx, 1)), ...
     D_initial_u.shift_all_indices(reshape((2 * N * nx + 1) : (2 * N * nx + nu), nu, 1))];

convex_constraints_gn = D.get_constraint_functions;
convex_constraints_eq_gn = D_eq.get_constraint_functions;

lambda = 1.08 * norm(Q);

%%
S_z = blkdiag(kron(eye(N), prob_3DoF.scaling.S_x), kron(eye(N), prob_3DoF.scaling.S_x), kron(eye(N), prob_3DoF.scaling.S_u), prob_3DoF.scaling.S_p);
c_z = [kron(ones([N, 1]), prob_3DoF.scaling.c_x); kron(ones([N, 1]), prob_3DoF.scaling.c_x); kron(ones([N, 1]), prob_3DoF.scaling.c_u); prob_3DoF.scaling.c_p];
z_guess_hat = S_z \ (z_guess - c_z);

Qtilde = S_z' * Q * S_z;
qtilde = S_z' * (Q' * c_z + q);
Htilde = H * S_z;
htilde = -H * c_z + h;
ctilde = (1/2 * c_z' * Q + q') * c_z + c;
lambdatilde = 1.08 * norm(Qtilde);

%%
Dxi = D_xi.shift_all_indices(reshape((N * nx + 1) : (2 * N * nx), nx, N));
Du = D_u.shift_all_indices(reshape((2 * N * nx + 1) : (2 * N * nx + N * nu), nu, N));
Ds = D_s.shift_all_indices(2 * N * nx + N * nu + 1);
Dfinal = D_final_x.shift_all_indices(reshape(((N - 1) * nx + N * nx + 1) : (2 * N * nx), nx, 1));
Dinitialu = D_initial_u.shift_all_indices(reshape((2 * N * nx + 1) : (2 * N * nx + nu), nu, 1));

%%
% z_f = z_guess((2 * N * nx + 1) : (2 * N * nx + N * nu))
% z_p = Du(1:3).project_onto_cones(1, z_guess);
% z_pf = z_p((2 * N * nx + 1) : (2 * N * nx + N * nu))

%% Test Customized Preconditioning
% General
t1 = tic;
[qhat_gn, Hhat_gn, L_gn, L_inv_gn] = L2_hypersphere_preconditioning(Q, q, H, lambda);
[qhattilde_gn, Hhattilde_gn, Ltilde_gn, Ltilde_inv_gn] = L2_hypersphere_preconditioning(Qtilde, qtilde, Htilde, lambdatilde);
t_precond_gn = toc(t1);

% Custom 
t1 = tic;
[qhat, Ahat_1, Ahat_2, Bhat_minus, Bhat_plus, Shat, l, l_inv] = custom_L2_hypersphere_preconditioning(w_vse, w_tr, w_tr_s, q_oc, A_k, B_k_minus, B_k_plus, S_k, lambdatilde);
[qhat, Ahat_1, Ahat_2, Bhat_minus, Bhat_plus, Shat, l, l_inv] = custom_L2_hypersphere_preconditioning(w_vse, w_tr, w_tr_s, q_oc, A_k, B_k_minus, B_k_plus, S_k, lambdatilde);
t_precond = toc(t1);


% Compare general vs custom
fprintf("Custom L2 Hypersphere Preconditioning %.3f%s Faster than General\n", (t_precond_gn - t_precond) / t_precond_gn * 100, "%")

%% Test Customized Power Iteration
% General
t1 = tic;
sigma_gn = power_iteration(Hhat_gn, z_guess);
sigmatilde_gn = power_iteration(Hhattilde_gn, z_guess_hat);
t_power_gn = toc(t1);


% Custom
t1 = tic;
sigma = custom_power_iteration(Ahat_1, Ahat_2, Bhat_minus, Bhat_plus, Shat, l_inv.x1, l_inv.x2, guess.x, guess.x, guess.u, guess.p);
sigma = custom_power_iteration(Ahat_1, Ahat_2, Bhat_minus, Bhat_plus, Shat, l_inv.x1, l_inv.x2, guess.x, guess.x, guess.u, guess.p);
t_power = toc(t1);

% Exact
t1 = tic;
sigma_gn_exact = norm(Hhat_gn) ^ 2 * 1.1;
sigma_gn_exact = norm(Hhat_gn) ^ 2 * 1.1;
t_power_exact = toc(t1);

% Compare general vs custom
%assert(abs(sigma_gn - sigma) / sigma < 1e-8 && abs(sigma_gn - sigma_gn_exact) / sigma_gn_exact < 1e-2, "sigma for Hhat not matching")

fprintf("Custom Power Iteration %.3f%s Faster than General \n", (t_power_gn - t_power) / t_power_gn * 100, "%")
fprintf("Custom Power Iteration %.3f%s Faster than Exact \n", (t_power_exact - t_power) / t_power_exact * 100, "%")

%% Solve with CVX First
[z_star_cvx, sol_info_cvx] = cvx_PIPG(Q, q, c, H, h, convex_constraints_gn, convex_constraints_eq_gn);

%%
[z_star_cvx, sol_info_cvx] = cvx_PIPG_scaled(Q, q, c, H, h, convex_constraints_gn, convex_constraints_eq_gn, S_z, c_z);

%%
[z_star_cvx, sol_info_cvx] = cvx_PIPG_scaled_constraints_only(Qtilde, qtilde, ctilde, Htilde, htilde, convex_constraints_gn, convex_constraints_eq_gn, S_z, c_z);

%%
[z_star_cvx_g, sol_info_cvx_g] = general_cvx_PIPG(z_guess, H, h, convex_constraints_gn, convex_constraints_eq_gn, nx, nu, np, N, w_tr, w_tr_s, w_vse, w_m);

%%
x_i_0 = 1;
x_i_f = N * nx;
x_i = x_i_0 : x_i_f;

xi_i_0 = x_i_f + 1;
xi_i_f = x_i_f + N * nx;
xi_i = xi_i_0 : xi_i_f;

u_i_0 = xi_i_f + 1;
u_i_f = xi_i_f + N * nu;
u_i = u_i_0 : u_i_f;

s_i_0 = u_i_f + 1;
s_i_f = u_i_f + np;
s_i = s_i_0 : s_i_f;

val_m = obj_function(z_star_cvx, z_guess, x_i, xi_i, u_i, s_i, w_tr, w_tr_s, w_vse, w_m, Q_state, Q_u, Q_s, Q, q, c)
val_cvx = obj_function_cvx(z_star_cvx, z_guess, x_i, xi_i, u_i, s_i, w_tr, w_tr_s, w_vse, w_m)
val_Q = 1/2 * z_star_cvx' * Q * z_star_cvx + q' * z_star_cvx + c

%%
cval_cvx = D.eval_constraint_functions(z_star_cvx);

sparse(H * z_star_cvx - h)

%%
x = reshape(z_star_cvx(1 : (N * nx)), nx, N);
xi = reshape(z_star_cvx((N * nx + 1) : (N * nx * 2)), nx, N); 
u_p = reshape(z_star_cvx((N * nx * 2 + 1) : (N * nx * 2 + N * nu)), nu, N); 
s = z_star_cvx(N * nx * 2 + N * nu + 1); 
%%
u = u_p;
for k = 1 : N
    u(:, k) = [cos(u_p(2, k)) -sin(u_p(2, k)); sin(u_p(2, k)) cos(u_p(2, k))] * [u_p(1, k); 0];
end
u = [u; u_p(1, :)];

figure
plot_3DoF_trajectory(t_k, x, u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);

figure
plot_3DoF_time_histories(t_k, x, u);

%%
J_tr = w_tr * (sum(vecnorm(x - guess.x, 1, 2)))
J_tr_s = w_tr_s * norm(s - guess.p)
J_vc = w_vse * sum(vecnorm(x - xi, 1, 2))

J_state = 1/2 * [x(:); xi(:)]' * Q_state * [x(:); xi(:)]
J_u = 1/2 * [u_p(:)]' * Q_u * [u_p(:)]
J_s = 1/2 * s' * Q_s * s;

J = 1/2 * z_star_cvx' * Q * z_star_cvx + q' * z_star_cvx + c

%%
[x_disc_sol] = prob_3DoF.disc_prop(u_p, s);
[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(u_p, s);

figure
comparison_plot_3DoF_trajectory({guess.x, x_disc_sol, x_cont_sol, x}, ["Guess", "Discrete Propagation", "Continuous Propagation", "Solution Output"], glideslope_angle_max, linestyle = [":", ":", "-", "--"], title = "3DoF Solution Comparison")

%%
x_ref = guess.x;
u_ref = guess.u;
s_ref = guess.p;

%%
[prob_3DoF, Delta_disc] = prob_3DoF.discretize(x_ref, u_ref, s_ref);

sum(Delta_disc)

A_k = prob_3DoF.disc.A_k;
B_k_plus = prob_3DoF.disc.B_plus_k;
B_k_minus = prob_3DoF.disc.B_minus_k;
S_k = prob_3DoF.disc.E_k;
d_k = prob_3DoF.disc.c_k;

c_constraints = [];
c_constraints.x = D_struct.xi.get_constraint_functions;
c_constraints.u = D_struct.u.get_constraint_functions;
c_constraints.s = D_struct.s.get_constraint_functions;
c_x_i = D_struct.ic_xi.get_constraint_function;
c_x_f = D_struct.tc_xi.get_constraint_function;
c_u_i = D_struct.ic_u.get_constraint_function;
[x, xi, u_p, s, sol_info] = most_general_cvx_PIPG(x_ref, u_ref, s_ref, A_k, B_k_minus, B_k_plus, S_k, d_k, c_x_i, c_x_f, c_u_i, c_constraints, nx, nu, np, N, w_tr, w_tr_s, w_vse, w_m);

u = u_p;
for k = 1 : N
    u(:, k) = [cos(u_p(2, k)) -sin(u_p(2, k)); sin(u_p(2, k)) cos(u_p(2, k))] * [u_p(1, k); 0];
end
u = [u; u_p(1, :)];

figure
plot_3DoF_trajectory(t_k, x, u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);
figure
plot_3DoF_trajectory(t_k, xi, u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);

figure
plot_3DoF_time_histories(t_k, x, u);

x_ref = x;
u_ref = u_p;
s_ref = s;

[x_disc_sol] = prob_3DoF.disc_prop(u_p, s);
[t_cont_sol, x_cont_sol, u_cont_sol] = prob_3DoF.cont_prop(u_p, s);

figure
comparison_plot_3DoF_trajectory({x_ref, x_disc_sol, x_cont_sol, x}, ["Guess", "Discrete Propagation", "Continuous Propagation", "Solution Output"], glideslope_angle_max, linestyle = [":", ":", "-", "--"], title = "3DoF Solution Comparison")

%%
obj_function_x(x, xi, u_p, s, x_ref, u_ref, s_ref, w_tr, w_tr_s, w_vse, w_m)

%%
[z_star_cvx, sol_info_cvx] = cvx_PIPG(Q, q, c, H, h, convex_constraints_gn, convex_constraints_eq_gn);


%%
[z_star_cvx_pc, sol_info_cvx_pc] = cvx_PIPG_pc(lambda, qhat_gn, Hhat_gn, h, L_inv_gn, convex_constraints_gn, convex_constraints_eq_gn, eye(size(c_z, 1)), c_z * 0);

fprintf("CVX_pc Optimal value: %.5f, solution %.3f%s off from MOSEK\n", 1/2 * z_star_cvx_pc' * Q * z_star_cvx_pc + q' * z_star_cvx_pc + c, norm(z_star_cvx - z_star_cvx_pc) / norm(z_star_cvx) * 100, "%")

x = reshape(z_star_cvx_pc(1 : (N * nx)), nx, N);
xi = reshape(z_star_cvx_pc((N * nx + 1) : (N * nx * 2)), nx, N); 
u_p = reshape(z_star_cvx_pc((N * nx * 2 + 1) : (N * nx * 2 + N * nu)), nu, N); 
s = z_star_cvx_pc(N * nx * 2 + N * nu + 1); 

u = u_p;
for k = 1 : N
    u(:, k) = [cos(u_p(2, k)) -sin(u_p(2, k)); sin(u_p(2, k)) cos(u_p(2, k))] * [u_p(1, k); 0];
end
u = [u; u_p(1, :)];

figure
plot_3DoF_trajectory(t_k, x, u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);

figure
plot_3DoF_time_histories(t_k, x, u);

%%
[z_star_cvx_pc_sc, sol_info_cvx_pc_sc] = cvx_PIPG_pc(lambdatilde, qhattilde_gn, Hhattilde_gn, htilde, Ltilde_inv_gn, convex_constraints_gn, convex_constraints_eq_gn, S_z, c_z);

fprintf("CVX_pc Optimal value: %.5f, solution %.3f%s off from MOSEK\n", 1/2 * z_star_cvx_pc_sc' * Q * z_star_cvx_pc_sc + q' * z_star_cvx_pc_sc + c, norm(z_star_cvx - z_star_cvx_pc_sc) / norm(z_star_cvx) * 100, "%")

x = reshape(z_star_cvx_pc_sc(1 : (N * nx)), nx, N);
xi = reshape(z_star_cvx_pc_sc((N * nx + 1) : (N * nx * 2)), nx, N); 
u_p = reshape(z_star_cvx_pc_sc((N * nx * 2 + 1) : (N * nx * 2 + N * nu)), nu, N); 
s = z_star_cvx_pc_sc(N * nx * 2 + N * nu + 1); 

u = u_p;
for k = 1 : N
    u(:, k) = [cos(u_p(2, k)) -sin(u_p(2, k)); sin(u_p(2, k)) cos(u_p(2, k))] * [u_p(1, k); 0];
end
u = [u; u_p(1, :)];

figure
plot_3DoF_trajectory(t_k, x, u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);

figure
plot_3DoF_time_histories(t_k, x, u);

figure
comparison_plot_3DoF_time_histories({t_k, t_k}, {x, xi}, {u, u}, ["x", "\xi"], linestyle = ["-", "--"]);

%% Test Customized PIPG

j_max = 1000;
rho = 1;
j_check = 1;
omega = 1;
%[z_star, w_star, sol_info] = PIPG_orig_precond(qhat_gn, Hhat_gn, h, D, L_gn, L_inv_gn, lambda, sigma_gn, rho, tol_abs, tol_rel, tol_infeas, j_check, 900000000, j_max, z_star_cvx_pc_sc, Hhat_gn * (L_gn * z_star_cvx_pc_sc) - h, eye(size(S_z)), c_z * 0);
% General
%j_max = 100;
[z_star, what_star, sol_info] = PIPG_scaled_new(qhattilde_gn, Hhattilde_gn, htilde, D, Ltilde_gn, Ltilde_inv_gn, lambdatilde, sigmatilde_gn, omega, rho, tol_abs, tol_rel, tol_infeas, 1, j_max, S_z \ (z_star_cvx_pc_sc - c_z), Hhat_gn * (L_gn * z_star_cvx_pc_sc) - h, S_z, c_z, t_k, glideslope_angle_max, gimbal_max, T_min, T_max)
%PIPG Proportional Integral Projected Gradient convex conic solver);

x = reshape(z_star(1 : (N * nx)), nx, N);
xi = reshape(z_star((N * nx + 1) : (N * nx * 2)), nx, N); 
u_p = reshape(z_star((N * nx * 2 + 1) : (N * nx * 2 + N * nu)), nu, N); 
u = u_p;
for k = 1 : N
    u(:, k) = [cos(u_p(2, k)) -sin(u_p(2, k)); sin(u_p(2, k)) cos(u_p(2, k))] * [u_p(1, k); 0];
end
u = [u; u_p(1, :)];

s = z_star(N * nx * 2 + N * nu + 1); 

figure
plot_3DoF_trajectory(t_k, x, u, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);

figure
plot_3DoF_time_histories(t_k, x, u);

%%
figure
plot((1:numel(sol_info.what_inf_dj)) * j_check, abs(sol_info.what_inf_dj)); hold on
plot((1:numel(sol_info.zhat_inf_dj)) * j_check, abs(sol_info.zhat_inf_dj)); hold off
title("PIPG Iterations")
subtitle(sprintf("Status: %s, Primal Infeasible: %s, Dual Infeasible: %s", sol_info.solution_status, string(sol_info.primal_infeasible), string(sol_info.dual_infeasible)))
xlabel("Iteration")
ylabel("Convergence Metric")
yscale("log")
legend("w_{\infty\Delta j}", "z_{\infty\Delta j}")
grid on

%% Custom
what_ref = zeros([nx, N - 1]);

[x_star, xi_star, u_star, s_star, what_star, sol_info] = custom_PIPG(qhat, Ahat_1, Ahat_2, Bhat_minus, Bhat_plus, Shat, d_k, D_struct, l, l_inv, lambda, sigma, rho, tol_abs, tol_rel, j_check, j_max, guess.x, guess.x, guess.u, guess.p, what_ref);

plot_3DoF_trajectory(t_k, x_star, u_star, glideslope_angle_max, gimbal_max, T_min, T_max, step = 1);

figure
plot((1:numel(sol_info.what_inf_dj)) * j_check, abs(sol_info.what_inf_dj)); hold on
plot((1:numel(sol_info.zhat_inf_dj)) * j_check, abs(sol_info.zhat_inf_dj)); hold off
title("PIPG Iterations")
subtitle(sprintf("Status: %s, Primal Infeasible: %s, Dual Infeasible: %s", sol_info.solution_status, string(sol_info.primal_infeasible), string(sol_info.dual_infeasible)))
xlabel("Iteration")
ylabel("Convergence Metric")
yscale("log")
legend("w_{\infty\Delta j}", "z_{\infty\Delta j}")
grid on

%% Compare general vs custom
function [val] = obj_function_x(x, xi, u, s, x_ref, u_ref, s_ref, w_tr, w_tr_s, w_vse, w_m)
    val = 0;
    for k = 1 : size(x, 2)
        val = val + w_tr * (x(:, k) - x_ref(:, k))' * (x(:, k) - x_ref(:, k));
    end
    for k = 1 : size(u, 2) 
        val = val + w_tr * (u(:, k) - u_ref(:, k))' * (u(:, k) - u_ref(:, k));
    end
    for k = 1 : numel(s)
        val = val + w_tr_s * (s(k) - s_ref(k))' * (s(k) - s_ref(k));
    end
    for k = 1 : size(xi, 2)
        val = val + w_vse * (x(:, k) - xi(:, k))' * (x(:, k) - xi(:, k));
    end

    val = val / 2;

    val = val - w_m * x(7, end);
end


function [val] = obj_function_cvx(z, z_ref, x_i, xi_i, u_i, s_i, w_tr, w_tr_s, w_vse, w_m)
    val = 0;
    for k = 1 : numel(x_i)
        val = val + w_tr * (z(x_i(k)) - z_ref(x_i(k))) * (z(x_i(k)) - z_ref(x_i(k)));
    end
    for k = 1 : numel(u_i) 
        val = val + w_tr * (z(u_i(k)) - z_ref(u_i(k))) * (z(u_i(k)) - z_ref(u_i(k)));
    end
    for k = 1 : numel(s_i)
        val = val + w_tr_s * (z(s_i(k)) - z_ref(s_i(k))) * (z(s_i(k)) - z_ref(s_i(k)));
    end
    for k = 1 : numel(xi_i)
        val = val + w_vse * (z(x_i(k)) - z(xi_i(k))) * (z(x_i(k)) - z(xi_i(k)));
    end

    val = val / 2;

    val = val - w_m * z(x_i(end));
end

function [val] = obj_function(z, z_ref, x_i, xi_i, u_i, s_i, w_tr, w_tr_s, w_vse, w_m, Q_state, Q_u, Q_s, Q, q, c)
    val = 0;

    val = val + w_tr * norm(z(x_i) - z_ref(x_i)) .^ 2;
    val = val + w_tr * norm(z(u_i) - z_ref(u_i)) .^ 2;
    val = val + w_tr_s * norm(z(s_i) - z_ref(s_i)) .^ 2;
    val = val + w_vse * norm(z(x_i) - z(xi_i)) .^ 2;

    val2 = z([x_i, xi_i])' * Q_state * z([x_i, xi_i]) + z(u_i)' * Q_u * z(u_i) + z(s_i)' * Q_s * z(s_i) ...
        - w_tr * 2 * z(x_i)' * z_ref(x_i) - w_tr * 2 * z(u_i)' * z_ref(u_i) - w_tr_s * 2 * z(s_i)' * z_ref(s_i) ...
        + w_tr * z_ref(x_i)' * z_ref(x_i) + w_tr * z_ref(u_i)' * z_ref(u_i) + w_tr_s * z_ref(s_i)' * z_ref(s_i);

    val3 = 1/2 * z' * Q * z - q' * z + c;

    val = val / 2;

    val = val - w_m * z(x_i(end));
end
% 
% function [val] = obj_function_1(z, z_ref, x_i, xi_i, u_i, s_i, w_tr, w_tr_s, w_vse, w_m)
%     val = 0;
% 
%     val = val + w_tr * norm(z(x_i) - z_ref(x_i)) .^ 2;
%     val2 = w_tr * ;
%     val = val + w_tr * norm(z(u_i) - z_ref(u_i)) .^ 2;
%     val = val + w_tr_s * norm(z(s_i) - z_ref(s_i)) .^ 2;
%     val = val + w_vse * norm(z(x_i) - z(xi_i)) .^ 2;
% 
%     val = val / 2;
% 
%     val = val - w_m * z(x_i(end));
% end