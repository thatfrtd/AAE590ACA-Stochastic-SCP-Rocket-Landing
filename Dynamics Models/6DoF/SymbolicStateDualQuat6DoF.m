syms g L alpha_ME alpha_RCS;
I = sym("I", [3; 1]);
c = [g; L; I; alpha_ME; alpha_RCS];

% states

v_I = sym("v_I", [3;1]);
q_var = sym("q", [4;1]);
w_B = sym("w_B", [3;1]);
m = sym("m");

qq = sym("qq", [8;1]); %[q; 1/2 * q_mul([r_I; 0], q)];
q = qq(1:4);
%w_tilde = sym("w_tilde", [8; 1]); %[[w_B; 0]; q_mul(q_conj(q), q_mul([v_I; 0], q))];
v_B = sym("v_B", [3;1]);%w_tilde(5:7); quat_rot(q, v_I)
w = [w_B; v_B];
w_tilde = [w_B; 0; v_B; 0];

x = [m; qq; w];

% controls
T = sym("T", [1;1]); % thrust magnitude 
delta = sym("delta", [1;1]); % gimbal deflection angle
phi = sym("delta", [1;1]); % gimbal azimuth angle
tau = sym("gamma", [3;1]); % RCS moment

u = [T; delta; phi; tau];

% parameters
tf = sym("tf");
p = [tf];

Tau_B = [T * sin(delta) * cos(phi); T * sin(delta) * sin(phi); T * cos(delta); tau];

% intermediate variables
J_mat = [zeros(3, 3), m * eye(3); m * diag(I), zeros(3, 3)];
J_inv = inv(J_mat);

Phi = [eye(3), zeros(3, 3); skew([0; 0; -L]), eye(3)];

g_I = [0; 0; -g];
g_B = quat_rot(q, g_I);

g_B_vec = [zeros([3, 1]); g_B];

% q dot
qq_dot = 1/2 * qq_mul(qq, w_tilde);

% w_dot
w_cross = [skew(w_B), zeros(3, 3); skew(v_B), skew(w_B)];
w_dot = J_inv * (-w_cross * J_mat * w + Phi * Tau_B) + g_B_vec;

% m dot
m_dot = -(alpha_ME * T + alpha_RCS * sqrt(tau(1) ^ 2 + tau(2) ^ 2 + tau(3) ^ 2) / L);

x_dot = [m_dot; qq_dot; w_dot] * tf;
%j_a = jacobian(x_dot, x);
%j_b = jacobian(x_dot, u);

vars = [{x}; {u}; {p}; {g; L; I; alpha_ME; alpha_RCS}];

% Create equations of motion function for optimizer
matlabFunction(x_dot,"File","Dynamics Models/6DoF/SymDynamicsDualQuat6DoF","Vars",vars);

% % Create equations of motion block for Simulink model
% open("EoM_Euler6DoF.slx")
% matlabFunctionBlock('EoM_Euler6DoF/SymDynamicsEuler6DoF',x_dot,'Vars',vars);
% close_system("EoM_Euler6DoF.slx", 1)
% 
% % Create Jacobian functions for Kalman filter
% matlabFunction(j_a,"File","Dynamics Models/6DoF/SymXJacobianEuler6DoF","Vars",vars);
% matlabFunction(j_b,"File","Dynamics Models/6DoF/SymUJacobianEuler6DoF","Vars",vars);
% 

%xdot_func6dof = matlabFunction(state);
%j_a_func6dof = matlabFunction(j_a);
%j_b_func6dof = matlabFunction(j_b);

%save("dynamics6dof.mat", "xdot_func6dof", "j_a_func6dof", "j_b_func6dof")
%save("SymDynamics6DoF.m", "xdot_func6dof")
%save("SymUJacobian6DoF.m", "j_b_func6dof")
%save("SymXJacobian6DoF.m", "j_a_func6dof")

