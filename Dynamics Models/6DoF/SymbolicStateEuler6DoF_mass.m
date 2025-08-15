% constants (same as 3DoF)
g = 1.62e-3;%3.7114e-3; % [km / s2]

syms L alpha;
I = sym("I", [3; 1]);
c = [g; L; I; alpha];

% states

r = sym("r", [3;1]);
v = sym("v", [3,1]);
theta = sym("theta", [3;1]);
w = sym("w", [3;1]);
m = sym("m", [1;1]);

x = [r; v; theta; w; m];

% controls

T = sym("T", [3;1]);
thrust_mag = sym("thrust_mag", 1);
gamma = sym("gamma", 1);

u = [T; thrust_mag; gamma];

% calculate theta dot
b_inverse =  (1/sin(theta(2))) .* [0 sin(theta(3)) cos(theta(3));
     0 sin(theta(2))*cos(theta(3)) -sin(theta(2))*sin(theta(3));
     sin(theta(2)) -cos(theta(2))*sin(theta(3)) -cos(theta(2))*cos(theta(3))];
thetadot = b_inverse * w;

% r dot
r_dot = v;

Rxtheta1 = [1 0 0; 0 cos(theta(1)) sin(theta(1)); 0 -sin(theta(1)) cos(theta(1))];
Rxtheta2 = [cos(theta(2)) 0 -sin(theta(2)); 0 1 0; sin(theta(2)) 0 cos(theta(2))];
Rxtheta3 = [1 0 0; 0 cos(theta(3)) sin(theta(3)); 0 -sin(theta(3)) cos(theta(3))];
C_be = Rxtheta3 * Rxtheta2 * Rxtheta1;
C_eb = C_be.';
T_e = C_eb * T;
% v dot
v_dot = T_e / m * cos(gamma) - [0; 0; g];
%
M = thrust_mag .* sin(gamma) .* [1; 0; 0] + cross([-L; 0; 0], cos(gamma) * T);
% w dot
w_dot = (M + (I([2; 3; 1]) - I([3; 1; 2])) .* w([2; 3; 1]) .* w([3; 1; 2])) ./ I([1; 2; 3;]);

% m dot
m_dot = -alpha * thrust_mag;

x_dot = [r_dot; v_dot; thetadot; w_dot; m_dot];
%j_a = jacobian(x_dot, x);
%j_b = jacobian(x_dot, u);

vars = [{x}; {u}; {L; I; alpha}];

% Create equations of motion function for optimizer
matlabFunction(x_dot,"File","Dynamics Models/6DoF/SymDynamicsEuler6DoF_mass","Vars",vars);

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