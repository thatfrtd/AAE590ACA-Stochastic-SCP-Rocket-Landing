%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 19 April, 2025
% Description: linear 2DoF rocket landing dynamics
% Most Recent Change: 19 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

syms alpha tf g;
t = sym("t");
r = sym("r", [3, 1]);
v = sym("v", [3, 1]);
m = sym("m");
x = [r;v;m];

thrust = sym("thrust", [3, 1]);
thrust_mag = sym("thrust_mag", 1);
u = [thrust; thrust_mag];

p = [tf];

A = [zeros(3), eye(3), zeros(3, 1); zeros(3), zeros(3), zeros(3, 1); zeros(1, 3), zeros(1, 3), 0];
B = [zeros(3), zeros(3, 1); eye(3), zeros(3, 1); zeros(1, 3), -alpha];
c = [zeros(3, 1); [0; 0; -g]; 0];

xdot = (A * x + B * u + c) * tf;

%j_a = jacobian(xdot, x);
%j_b = jacobian(xdot, u);

% Create equations of motion function for optimizer
matlabFunction(xdot,"File","Dynamics Models/3DoF_linear/SymDynamics3DoF_linear","Vars", [{t}; {x}; {u}; {p}; {alpha}; {g}]);

% Create equations of motion block for Simulink model
%matlabFunctionBlock('EoM_3DoF/SymDynamics3DoF',xdot,'Vars',[x; u; mass; L; I])

% Create Jacobian functions for Kalman filter
%matlabFunction(j_a,"File","3DoF/SymXJacobian3DoF","Vars",[x; u; mass; L; I]);
%matlabFunction(j_b,"File","3DoF/SymUJacobian3DoF","Vars",[x; u; mass; L; I]);