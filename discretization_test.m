%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% HW1 Q1
% Author: Travis Hastreiter 
% Created On: 25 January, 2024
% Description: Simulate orbits with J2 and SRP perturbations
% Most Recent Change: 30 January, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Initialize
% Vehicle Parameters
Isp = 225; % [s]
alpha = 1 / (9.81e-3 * Isp); % [s / km]
T_min = 6000e-3; % [kg km / s2]
T_max = 22500e-3; % [kg km / s2]
tau_max = 100e-6; % [kg km2 / s2]
I = [19150; 13600; 13600] * (1e-3) ^ 2; % [kg km2] ASSUMING CONSTANT MOMENT OF INERTIA
L = 0.25e-3; % [km] Distance from CoM to nozzle
m_dry = 2100; % [kg]
m_wet = 1150; % [kg]
m_0 = m_dry + m_wet;
gimbal_max = deg2rad(10); % [rad]
g = 1.62e-3;
time_min_max_thrust = 3; % [s] time to throttle from min to max thrust
max_gimbal_rate = deg2rad(10); % [rad / s] max rate of change of gimbal angle
gimbal_max_final = deg2rad(5); % [rad]
gimbal_STC_trigger_height = 50 * 1e-3; % [km]

vehicle = Vehicle(m_dry, L, L * 3, gimbal_max, T_min, T_max, I = I);

% Problem Parameters
tf = 35; % [s]
N = 10; % []
r_0 = [250; -100; 433] * 1e-3; % [km]
v_0 = [0; 0; -35] * 1e-3; % [km / s]
theta_0 = [deg2rad(0); deg2rad(90); deg2rad(0)]; % [rad]
initial_roll = deg2rad(0);
R_0 = make_R(initial_roll, 3) * angle2dcm(theta_0(1), theta_0(2), theta_0(3));
q_0 = qexp(RLog(R_0));
w_0 = deg2rad([0; 0; 0]); % [rad / s]
glideslope_angle_max = deg2rad(65); % [rad]
pitch_max = deg2rad(45); % [rad]
angvel_max = deg2rad(20); % [rad]

theta_f = [0; deg2rad(90); 0]; % [rad]
R_f = angle2dcm(theta_f(1), theta_f(2), theta_f(3));
q_f = qexp(RLog(R_f));

x_0 = [r_0; v_0; q_0; w_0; m_0];
x_f = [[0; 0; 30] * 1e-3; [0; 0; -1] * 1e-3; q_f; zeros(3, 1)];

tspan = [0, tf];
t_k = linspace(tspan(1), tspan(2), N);
delta_t = t_k(2) - t_k(1);

u_hold = "FOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

initial_guess = "straight line";

parser = "CVXPY";

nx = 14; 
nu = 4;
np = 0;

scale_hint.x_max = [max(r_0) * ones([3, 1]); max(abs(v_0)) * ones([3, 1]); ones([4, 1]); max(abs([w_0; deg2rad(10)])) * ones([3, 1]); m_0];
scale_hint.x_min = [-max(r_0) * ones([3, 1]); -max(abs(v_0)) * ones([3, 1]); -ones([4, 1]); -max(abs([w_0; deg2rad(10)])) * ones([3, 1]); m_dry];
scale_hint.u_max = [T_max; gimbal_max; gimbal_max; pi / 4];
scale_hint.u_min = [T_min; -gimbal_max; -gimbal_max; -pi / 4];

% Get dynamics
f = @(t, x, u, p) SymDynamicsQuat6DoF_localrot_alphabeta(x, u, L, I, alpha, g);
%%
% Create objective
Q = diag(1 ./ scale_hint.x_max .^ 2);
R = diag(1 ./ scale_hint.u_max .^ 2);
tracking_cost = @(x, u, x_ref, u_ref) quadform(x - x_ref, Q) + quadform(u - u_ref, R);
MPC_objective = @(x, u, x_ref, u_ref, p) tracking_cost(x, u, x_ref, u_ref) - x(14, end) / m_0 + sum(abs(u(4, :)));
%%
% Linearized matrices
t_sym = sym("t");
x_sym = sym("x", [nx, 1]);
u_sym = sym("u", [nu, 1]);
p_sym = sym("p", [np, 1]);

A = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), x_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
B = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), u_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
S = matlabFunction(jacobian(f(t_sym, x_sym, u_sym, p_sym), p_sym),"Vars", [{t_sym}; {x_sym}; {u_sym}; {p_sym}]);
%%
%[text] Test discretization using multiple different methods and compare speed, discretization accuracy, propagation accuracy for timesteps and full time preiod, and solution accuracy for MPC (?)
%[text] ## Exact Discretization using ODE78
A_k = prob_6DoF.disc.A_k; %[output:3d061686]
B_k_plus = prob_6DoF.disc.B_plus_k;
B_k_minus = prob_6DoF.disc.B_minus_k;
S_k = prob_6DoF.disc.E_k;
d_k = prob_6DoF.disc.c_k;
%[text] ## Euler Discretization
N_sub = 100;
[A_k_rk4, B_k_plus_rk4, B_k_minus_rk4, S_k_rk4, d_k_rk4, Delta_rk4] = discretize_error_dynamics_FOH_RK4(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk4 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk4 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk4 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk4 - S_k), "all");
d_err = sum(pagenorm(d_k_rk4 - d_k), "all");
Delta_err = norm(Delta_rk4 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.5f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);
%[text] ## RK4 Discretization
N_sub = 100;
[A_k_rk4, B_k_plus_rk4, B_k_minus_rk4, S_k_rk4, d_k_rk4, Delta_rk4] = discretize_error_dynamics_FOH_RK4(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk4 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk4 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk4 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk4 - S_k), "all");
d_err = sum(pagenorm(d_k_rk4 - d_k), "all");
Delta_err = norm(Delta_rk4 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.5f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);
%[text] ## RK7 Discretization
N_sub = 100;
[A_k_rk4, B_k_plus_rk4, B_k_minus_rk4, S_k_rk4, d_k_rk4, Delta_rk4] = discretize_error_dynamics_FOH_RK4(f, A, B, S, N, tspan, guess.x, guess.u, guess.p, N_sub);

A_err = sum(pagenorm(A_k_rk4 - A_k), "all");
B_minus_err = sum(pagenorm(B_k_minus_rk4 - B_k_minus), "all");
B_plus_err = sum(pagenorm(B_k_plus_rk4 - B_k_plus), "all");
S_err = sum(pagenorm(S_k_rk4 - S_k), "all");
d_err = sum(pagenorm(d_k_rk4 - d_k), "all");
Delta_err = norm(Delta_rk4 - Delta_disc);
fprintf("A: %.3f, B-: %.3f, B+: %.3f, S: %.3f, d: %.3f, Delta: %.5f\n", A_err, B_minus_err, B_plus_err, S_err, d_err, Delta_err);
%%
%[text] ## Helper Functions
function [q] = qexp(tau)
    theta = norm(tau);
    u = tau / theta;
    q = [u * sin(theta / 2); cos(theta / 2)];
end

function [tau] = RLog(R)
    theta = acos((trace(R) - 1) / 2);
    u = vee(R - R') / (2 * sin(theta));

    tau = theta * u;

    function [tau] = vee(tau_hat)
        tau = [tau_hat(3, 2); tau_hat(1, 3); tau_hat(2, 1)];
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
%[output:3d061686]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unable to resolve the name 'prob_6DoF.disc.A_k'."}}
%---
