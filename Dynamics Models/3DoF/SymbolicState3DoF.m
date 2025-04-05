%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA - AC TrajOpt
% 3 Degree of Freedom Symbolic Equations of Motion Script
% Author: Lucas Mohler, Travis Hastreiter 
% Created On: 2 November, 2024
% Description: Creates symbolic equations of motion (EoM) for 3DoF planar 
% rocket landing problem. Creates function file and Simulink block for the 
% EoMs and function files for the Jacobians.
% Most Recent Change: Travis Hastreiter 5 November, 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81;

syms mass L I;
r = sym("r", [2;1]);
v = sym("v", [2;1]);
theta = sym("theta", 1);
w = sym("w", 1);
x = [r;v;theta;w];

thrust = sym("thrust", 1);
gimbal = sym("gimbal", 1);
u = [thrust;gimbal];

rdot = v;
thetadot = w;
wdot = -L*thrust*sin(gimbal)/I;

totalTheta = theta+gimbal;
rotationMatrix = [cos(totalTheta) sin(totalTheta); sin(totalTheta) cos(totalTheta)];
T_e = rotationMatrix * [thrust; 0];
vdot = T_e/mass - [0; g];

xdot = [rdot; vdot; thetadot; wdot];

j_a = jacobian(xdot, x);
j_b = jacobian(xdot, u);

% Create equations of motion function for optimizer
matlabFunction(xdot,"File","3DoF/SymDynamics3DoF","Vars",[x; u; mass; L; I]);

% Create equations of motion block for Simulink model
matlabFunctionBlock('EoM_3DoF/SymDynamics3DoF',xdot,'Vars',[x; u; mass; L; I])

% Create Jacobian functions for Kalman filter
matlabFunction(j_a,"File","3DoF/SymXJacobian3DoF","Vars",[x; u; mass; L; I]);
matlabFunction(j_b,"File","3DoF/SymUJacobian3DoF","Vars",[x; u; mass; L; I]);