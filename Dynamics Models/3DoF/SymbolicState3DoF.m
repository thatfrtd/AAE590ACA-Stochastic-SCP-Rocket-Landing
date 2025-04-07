%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF rocket landing dynamics
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 3.7114e-3; % [km / s2]

syms mass L I;
t = sym("t");
r = sym("r", [2, 1]);
v = sym("v", [2, 1]);
theta = sym("theta", 1);
w = sym("w", 1);
x = [r;v;theta;w];

thrust = sym("thrust", [2, 1]);
thrust_mag = sym("thrust_mag", 1);
u = [thrust; thrust_mag];

rdot = v;
thetadot = w;
M = cross([-L; 0; 0], [thrust; 0]);
wdot = M(3) / I;

rotationMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
T_e = rotationMatrix * thrust;
vdot = T_e / mass - [0; g];

xdot = [rdot; vdot; thetadot; wdot];

%j_a = jacobian(xdot, x);
%j_b = jacobian(xdot, u);

% Create equations of motion function for optimizer
matlabFunction(xdot,"File","Dynamics Models/3DoF/SymDynamics3DoF","Vars", [{t}; {x}; {u}; {mass; L; I}]);

% Create equations of motion block for Simulink model
%matlabFunctionBlock('EoM_3DoF/SymDynamics3DoF',xdot,'Vars',[x; u; mass; L; I])

% Create Jacobian functions for Kalman filter
%matlabFunction(j_a,"File","3DoF/SymXJacobian3DoF","Vars",[x; u; mass; L; I]);
%matlabFunction(j_b,"File","3DoF/SymUJacobian3DoF","Vars",[x; u; mass; L; I]);