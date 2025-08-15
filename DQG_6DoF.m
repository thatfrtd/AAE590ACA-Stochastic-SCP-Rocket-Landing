%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Values from paper
% Vehicle Parameters
g_0 = 9.81; % [m / s2]
Isp_ME = 300; % [s]
Isp_RCS = 200; % [s]
alpha_ME = 1 / (Isp_ME * g_0); % [s / m]
alpha_RCS = 1 / (Isp_RCS * g_0); % [s / m]
T_min = 600; % [kg m / s2]
T_max = 3000; % [kg m / s2] 
tau_max = 50; % [kg m2 / s2]
Tdot_max = 0.75 * (T_max - T_min); % [kg m / s3]
delta_max = deg2rad(5); % [rad]
deltadot_max = deg2rad(5); % [rad / s]
phidot_max = deg2rad(5); % [rad / s]
L = 1; % [m] Distance from CoM to nozzle
p_B = [0.5; 0; -sqrt(3) / 2]; % line of sight
m_i = 1500; % [kg]
m_f = 750; % [kg]
I = [4.2; 4.2; 0.6]; % [kg m2] ASSUMING CONSTANT MOMENT OF INERTIA

theta_max = deg2rad(90); % [rad]
w_max = deg2rad(5); % [rad / s]
v_max = 90; % [m / s]
h_min = 100; % [m]
rho_max = 1250; % [m]
rho_min = 500; % [m]

theta_STC_max = deg2rad(20); % [rad]
w_STC_max = deg2rad(1); % [rad / s]
v_STC_max = deg2rad(30); % [m / s]
mu_STC_max = deg2rad(2); % [rad]

vehicle = Vehicle(m_i, L, L * 3, delta_max, T_min, T_max, I = I);

% Problem Parameters
tf = 140; % [s] guess
N = 25; % []

g = 1.625; % [m / s2]

r_I_i = [3000; 600; 3000]; % [m]
r_I_f = [0; 0; 100]; % [m]
v_I_i = [-60; 30; -30]; % [m / s]
v_z_I_f = -2; % [m / s]
q_i = [-0.15; 0.3; -1; 1];
q_i = q_i ./ norm(q_i);
q_f = [0; 0; -1.25; 1]; 
q_f = q_f ./ norm(q_f);
w_B_i = [0; 0; 0]; % [rad / s]

m_ind = 1;
qq_ind = 2:9;
q_ind = 2:5;
r_ind = 6:9;
angvel_ind = 10:12;
v_ind = 13:15;
w_ind = 10:15;
T_ind = 1;
delta_ind = 2;
phi_ind = 3;
tau_ind = 4:6;
tf_ind = 1;

x_0 = [m_i; construct_dual_quat(q_i, r_I_i); w_B_i; quat_rot(q_i, v_I_i)];
%terminal_bc = @(x, p) [0; x(qq_ind(1:2)); (1/2 * q_mul_matrix([r_I_f; 0]) * [zeros(2, 2); eye(2)] - eye(4)) * x(qq_ind(3:8)); ... 
%                          zeros([5, 1]); v_z_I_f];
qq_f = construct_dual_quat(q_f, r_I_f);
x_f = [qq_f; zeros([5, 1]); v_z_I_f];
terminal_bc = @(x, p) [0; x(2:end) - x_f];

tspan = [0, 1];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "FOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

% PTR algorithm parameters
ptr_ops.iter_max = 30;
ptr_ops.iter_min = 2;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 2e2;
ptr_ops.w_tr = ones(1, Nu) * 3e2;
ptr_ops.w_tr_p = 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 1e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 1;

scale = false;

%% Get Dynamics
f = @(t, x, u, p) SymDynamicsDualQuat6DoF(x, u, p, g, vehicle.L, vehicle.I, alpha_ME, alpha_RCS);

%% Specify Constraints
psi = @(x) sign((rho_max - norm(2 * x(r_ind))) * (norm(2 * x(r_ind)) - rho_min));

% Convex state path constraints
z_I = [0; 0; 1];
glideslope_angle_max = deg2rad(45); % [rad] - not part of original problem
glideslope_constraint = {1:N, @(t, x, u, p) -x(qq_ind)' * [zeros(4, 4), q_mul_matrix([z_I; 0])'; q_mul_matrix([z_I; 0]), zeros(4, 4)] * x(qq_ind) + norm(2 * x(qq_ind(5:8))) * cos(glideslope_angle_max)};

min_tf_constraint = {1, @(t, x, u, p) 130 - p};
max_tf_constraint = {1, @(t, x, u, p) p - 150};

% For k = N
min_mass_constraint = {N, @(t, x, u, p) m_f - x(m_ind)};
state_convex_constraints = {min_mass_constraint, min_tf_constraint, max_tf_constraint};

% Convex control constraints
min_thrust_constraint = {1:N, @(t, x, u, p) T_min - u(T_ind)};
max_thrust_constraint = {1:N, @(t, x, u, p) u(T_ind) - T_max};
min_gimbal_constraint = {1:N, @(t, x, u, p) 0 - u(delta_ind)};
max_gimbal_constraint = {1:N, @(t, x, u, p) u(delta_ind) - delta_max};
min_phi_constraint = {1:N, @(t, x, u, p) 0 - u(phi_ind)};
max_phi_constraint = {1:N, @(t, x, u, p) u(phi_ind) - 2 * pi};
max_tau_constraint = {1:N, @(t, x, u, p) norm(u(tau_ind), Inf) - tau_max};
control_convex_constraints = {max_thrust_constraint, max_gimbal_constraint, min_gimbal_constraint, min_phi_constraint, max_phi_constraint, max_tau_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state constraints
max_angvel_constraint = {2 : (N - 1), @(t, x, u, p, x_ref, u_ref, p_ref, k) norm(x(angvel_ind), Inf) - max(-psi(x_ref(:, k)) * w_max, w_STC_max)};
max_vel_constraint = {2 : (N - 1), @(t, x, u, p, x_ref, u_ref, p_ref, k) norm(x(v_ind), 2) - max(-psi(x_ref(:, k)) * v_max, v_STC_max)};
max_STC_tilt_angle_constraint = {2 : (N - 1), @(t, x, u, p, x_ref, u_ref, p_ref, k) (psi(x_ref(:, k)) >= 0) * (x_ref(q_ind(1:2), k)' * x(q_ind(1:2)) - norm(x_ref(q_ind(1:2), k)) * sin(theta_STC_max / 2))};
max_los_constraint = {2 : (N - 1), @(t, x, u, p, x_ref, u_ref, p_ref, k) (psi(x_ref(:, k)) >= 0) * (x_ref(qq_ind, k)' * [zeros(4, 4), q_conj_mul_matrix([p_B; 0])'; q_conj_mul_matrix([p_B; 0]), zeros(4, 4)] * x(qq_ind) + norm(2 * x(qq_ind(5:8))) * cos(mu_STC_max))}; % max line of sight angle
max_tilt_angle_constraint = {2 : (N - 1), @(t, x, u, p, x_ref, u_ref, p_ref, k) (psi(x_ref(:, k)) < 0) * (x_ref(q_ind(1:2), k)' * x(q_ind(1:2)) - norm(x_ref(q_ind(1:2), k)) * sin(theta_max / 2))};
min_altitude_constraint = {2 : (N - 1), @(t, x, u, p, x_ref, u_ref, p_ref, k) (psi(x_ref(:, k)) < 0) * (h_min - x_ref(qq_ind, k)' * [zeros(4, 4), q_conj_mul_matrix([z_I; 0])'; q_conj_mul_matrix([z_I; 0]), zeros(4, 4)] * x(qq_ind))}; % Minimum altitude
state_nonconvex_constraints = {max_angvel_constraint, max_vel_constraint, max_STC_tilt_angle_constraint, max_los_constraint, max_tilt_angle_constraint, min_altitude_constraint};

% Nonconvex control constraints
min_thrust_rate_constraint = {2:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) max(T_min, -Tdot_max * p_ref(tf_ind) / (N - 1) + u_ref(T_ind, k - 1)) - u(T_ind)};
max_thrust_rate_constraint = {2:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) u(T_ind) - min(T_max, Tdot_max * p_ref(tf_ind) / (N - 1) + u_ref(T_ind, k - 1))};
min_delta_rate_constraint = {2:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) max(0, -deltadot_max * p_ref(tf_ind) / (N - 1) + u_ref(delta_ind, k - 1)) - u(delta_ind)};
max_delta_rate_constraint = {2:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) u(delta_ind) - min(delta_max, deltadot_max * p_ref(tf_ind) / (N - 1) + u_ref(delta_ind, k - 1))};
min_phi_rate_constraint = {2:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) max(0, -phidot_max * p_ref(tf_ind) / (N - 1) + u_ref(phi_ind, k - 1)) - u(phi_ind)};
max_phi_rate_constraint = {2:N, @(t, x, u, p, x_ref, u_ref, p_ref, k) u(phi_ind) - min(2 * pi, phidot_max * p_ref(tf_ind) / (N - 1) + u_ref(phi_ind, k - 1))};
control_nonconvex_constraints = {min_thrust_rate_constraint, max_thrust_rate_constraint, min_delta_rate_constraint, max_delta_rate_constraint, min_phi_rate_constraint, max_phi_rate_constraint};

% Combine nonconvex constraints
nonconvex_constraints = [];

%% Specify Objective
if u_hold == "ZOH"
    min_fuel_objective = @(x, u, p, x_ref, u_ref, p_ref) -x(m_ind, N);
elseif u_hold == "FOH"
    min_fuel_objective = @(x, u, p, x_ref, u_ref, p_ref) -x(m_ind, N);
end

%% Create Guess

%sl_guess = guess_6DoF(x_0, x_f, N, Nu, delta_t, vehicle);

%%
[y, p, r] = dcm2angle(quat2rot(q_f), "ZYX");
rad2deg(y)

%%
qq_sclerp = sclerp(x_0(qq_ind)', qq_f', N);

q = qq_sclerp(1:4, :);

scl_guess = struct();

w_B = zeros(3, N);
v_I = interp1([0; 1], [v_I_i'; [0; 0; v_z_I_f]'], t_k)';

scl_guess.x = [qq_sclerp; w_B; quat_rot_array(q, v_I)];
%scl_guess.x([1, 3], :) = -scl_guess.x([1, 3], :);

scl_guess.p = tf;

if u_hold == "ZOH"
    scl_guess.x = [interp1([0, 1], [m_i, m_f], t_k(1:(end - 1))); scl_guess.x];
elseif u_hold == "FOH"
    scl_guess.x = [interp1([0, 1], [m_i, m_f], t_k); scl_guess.x];
end
%%
T_guess = scl_guess.x(1, :) * g;
T_guess = clip(T_guess, T_min, T_max);
delta_guess = zeros([1, N]);
phi_guess = zeros([1, N]);
tau_guess = 1e-3 * ones([3, N]);
scl_guess.u = [T_guess; delta_guess; phi_guess; tau_guess]; 

%%

%CasADi_sol = CasADi_solve_6DoF(x_0, qq_f, sl_guess.x, sl_guess.u, vehicle, N, delta_t, glideslope_angle_max);
%%
guess = scl_guess;
if u_hold == "ZOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "previous","extrap")';
elseif u_hold == "FOH"
    guess.u = interp1(t_k(1:size(guess.u, 2)), guess.u', t_k(1:Nu), "linear","extrap")';
end
%%
% figure
% plot_6DoF_trajectory(t_k, sl_guess.x, sl_guess.u, glideslope_angle_max, gimbal_max, T_min, T_max)
% 
% figure
% plot_6DoF_time_histories(t_k, sl_guess.x, sl_guess.u)

figure
plot_6DoFdqg_trajectory(t_k, guess.x, guess.u, glideslope_angle_max, delta_max, T_min, T_max, step = 1)
%%
figure
plot_6DoFdqg_time_histories(t_k, guess.x, guess.u)

%% Construct Problem Object
prob_6DoF = DeterministicProblem(x_0, qq_f, N, u_hold, 1, f, guess, convex_constraints, min_fuel_objective, terminal_bc = terminal_bc, scale = scale, nonconvex_constraints = nonconvex_constraints, discretization_method = "error", integration_tolerance = 1e-8);

%%
Delta = calculate_defect(prob_6DoF, guess.x, guess.u, guess.p);
norm(Delta)

%% Test Discretization on Initial Guess
[prob_6DoF, Delta_disc] = prob_6DoF.discretize(guess.x, guess.u, guess.p);

x_disc = prob_6DoF.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob_6DoF.cont_prop(guess.u, guess.p);

%%
figure
plot_6DoFdqg_trajectory(t_cont * guess.p, x_cont, u_cont, glideslope_angle_max, delta_max, T_min, T_max, step = 300)
%%
figure
plot_6DoFdqg_time_histories(t_cont * guess.p, x_cont, u_cont)


%%
figure
comparison_plot_6DoFdqg_trajectory({guess.x, x_cont, x_disc}, ["Guess", "Continuous Propagation", "Discrete Propagation"], glideslope_angle_max, linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")
%%
figure
comparison_plot_6DoFdqg_time_histories({t_k, t_cont, t_k}, {guess.x, x_cont, x_disc}, {guess.u, u_cont, guess.u}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Initial Guess")

%%
nx = 15;
nu = 6;
np = 1;

% Linearized matrices
t_sym = sym("t");
x_sym = sym("x", [nx, 1]);
u_sym = sym("u", [nu, 1]);
p_sym = sym("p", [np, 1]);

A = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
B = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), u_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
S = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), p_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);

A_k = prob_6DoF.disc.A_k;
B_k_plus = prob_6DoF.disc.B_plus_k;
B_k_minus = prob_6DoF.disc.B_minus_k;
S_k = prob_6DoF.disc.E_k;
d_k = prob_6DoF.disc.c_k;

%% Try RK4 Discretization
N_sub = 100;
[A_k_rk4, B_k_plus_rk4, B_k_minus_rk4, S_k_rk4, d_k_rk4, Delta_rk4] = discretize_error_dynamics_FOH_RK4(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk4 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk4 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk4 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk4 - S_k), "all");
d_err = sum(pagenorm(d_k_rk4 - d_k), "all");
Delta_err = norm(Delta_rk4 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.5f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);

%% Try RKV6(5) Discretization
N_sub = 100;
[A_k_rk65, B_k_plus_rk65, B_k_minus_rk65, S_k_rk65, d_k_rk65, Delta_rk65] = discretize_error_dynamics_FOH_RKV65(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk65 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk65 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk65 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk65 - S_k), "all");
d_err = sum(pagenorm(d_k_rk65 - d_k), "all");
Delta_err = norm(Delta_rk65 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.10f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);

%% Try RKV8(7) Discretization
N_sub = 100;
[A_k_rk87, B_k_plus_rk87, B_k_minus_rk87, S_k_rk87, d_k_rk87, Delta_rk87] = discretize_error_dynamics_FOH_RKV87(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk87 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk87 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk87 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk87 - S_k), "all");
d_err = sum(pagenorm(d_k_rk87 - d_k), "all");
Delta_err = norm(Delta_rk87 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.10f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);


%% Solve Problem with PTR
ptr_sol = ptr(prob_6DoF, ptr_ops);

%%
ptr_sol.converged_i = ptr_ops.iter_max;

%%
tiledlayout(1, 3)

nexttile
plot(0:ptr_sol.converged_i, [prob_6DoF.objective(prob_6DoF.guess.x, prob_6DoF.guess.u, prob_6DoF.guess.p), [ptr_sol.info.J]]); %hold on
%yline(CasADi_sol.objective); hold off
legend("PTR Iterations", "CasADi Solution")
title("Objective vs Iteration")
grid on

nexttile
plot(ptr_sol.delta_xp)
title("Stopping Criteria vs Iteration")
grid on

nexttile
plot(0:ptr_sol.converged_i, vecnorm(ptr_sol.Delta(:, 1:(ptr_sol.converged_i + 1)), 2, 1))
title("Defect Norm vs Iteration")
grid on

%%
i = ptr_sol.converged_i;

[t_cont_sol, x_cont_sol, u_cont_sol] = prob_6DoF.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i));

figure
plot_6DoFdqg_trajectory(t_k, ptr_sol.x(:, :, i), ptr_sol.u(:, :, i), glideslope_angle_max, delta_max, T_min, T_max, step = 1)

figure
comparison_plot_6DoFdqg_trajectory({guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, ["Guess", "Continuous Propagation", "Solution Output"], glideslope_angle_max, linestyle = [":", "-", "--", "-"], title = "Continuous vs Discrete Propagation of Solution")

figure
comparison_plot_6DoFdqg_time_histories({t_k, t_cont_sol, t_k}, {guess.x, x_cont_sol, ptr_sol.x(:, :, i)}, {guess.u, u_cont_sol, ptr_sol.u(:, :, i)}, ["Guess", "Cont", "Disc"], linestyle = [":", "-", "--"], title = "Continuous vs Discrete Propagation of Solution")

%% 
function R = quat2rot(Q)
q_0 = Q(4, :);
q_1 = Q(1, :);
q_2 = Q(2, :);
q_3 = Q(3, :);
R(1,1,:) = q_0.*q_0 + q_1.*q_1 - q_2.*q_2 - q_3.*q_3;
R(1,2,:) = 2.*(q_1.*q_2 - q_0.*q_3);
R(1,3,:) = 2.*(q_0.*q_2 + q_1.*q_3);
R(2,1,:) = 2.*(q_0.*q_3 + q_1.*q_2);
R(2,2,:) = q_0.*q_0 - q_1.*q_1 + q_2.*q_2 - q_3.*q_3;
R(2,3,:) = 2.*(q_2.*q_3 - q_0.*q_1);
R(3,1,:) = 2.*(q_1.*q_3 - q_0.*q_2);
R(3,2,:) = 2.*(q_0.*q_1 + q_2.*q_3);
R(3,3,:) = q_0.*q_0 - q_1.*q_1 - q_2.*q_2 + q_3.*q_3;
end

function Q = rot2quat(R)
q_0 = 1/2*sqrt(trace(R)+1);
q_r = 1/2*[...
    sign(R(3,2) - R(2,3))*sqrt(R(1,1) - R(2,2) - R(3,3) + 1)
    sign(R(1,3) - R(3,1))*sqrt(R(2,2) - R(3,3) - R(1,1) + 1)
    sign(R(2,1) - R(1,2))*sqrt(R(3,3) - R(1,1) - R(2,2) + 1)];
Q = [q_r.', q_0];
end

function PQ = quatProduct(P, Q)
p_0 = P(1);
q_0 = Q(1);
p_r = P(2:4);
q_r = Q(2:4);
scalarPart = p_0*q_0 - dot(p_r,q_r);
vectorPart = p_0*q_r + q_0*p_r + cross(p_r,q_r);
PQ = [vectorPart, scalarPart];
end