%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 19 April, 2025
% Description: linear 2DoF rocket landing dynamics
% Most Recent Change: 19 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 3.7114e-3; % [km / s2]

syms alpha;
t = sym("t");
r = sym("r", [2, 1]);
v = sym("v", [2, 1]);
m = sym("m");
x = [r;v;m];

thrust = sym("thrust", [2, 1]);
u = [thrust];

A = [zeros(2), eye(2), zeros(2, 1); zeros(2), zeros(2), zeros(2, 1); zeros(1, 2), zeros(1, 2), 0];
B = [zeros(2), zeros(2, 1); eye(2), zeros(2, 1); zeros(1, 2), -alpha];
c = [zeros(2, 1); [0; -g]; 0];

xdot = A * x + B * [u; norm(u)] + c;

%j_a = jacobian(xdot, x);
%j_b = jacobian(xdot, u);

% Create equations of motion function for optimizer
matlabFunction(xdot,"File","Dynamics Models/2DoF_linear/SymDynamics2DoF_linear_noumag","Vars", [{t}; {x}; {u}; {alpha}]);

% Create equations of motion block for Simulink model
%matlabFunctionBlock('EoM_3DoF/SymDynamics3DoF',xdot,'Vars',[x; u; mass; L; I])

% Create Jacobian functions for Kalman filter
%matlabFunction(j_a,"File","3DoF/SymXJacobian3DoF","Vars",[x; u; mass; L; I]);
%matlabFunction(j_b,"File","3DoF/SymUJacobian3DoF","Vars",[x; u; mass; L; I]);