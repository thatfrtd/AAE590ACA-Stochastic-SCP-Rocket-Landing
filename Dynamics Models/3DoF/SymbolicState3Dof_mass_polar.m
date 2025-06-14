%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter/Atharva
% Created On: 6 April, 2025
% Description: 3DoF rocket landing dynamics with changing mass
% Most Recent Change: 27 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g = 9.81e-3; % [km / s2]

syms L I alpha;
t = sym("t");
r = sym("r", [2, 1]);
v = sym("v", [2, 1]);
theta = sym("theta", 1);
w = sym("w", 1);
m = sym("m", 1);
x = [r;v;theta;w;m];

thrust_mag = sym("T");
gimbal = sym("delta");
u = [thrust_mag; gimbal];

tf = sym("tf");
p = [tf];

thrust = [cos(gimbal) -sin(gimbal); sin(gimbal) cos(gimbal)] * [thrust_mag; 0];

rdot = v;
thetadot = w;
M = cross([-L; 0; 0], [thrust; 0]);
wdot = M(3) / I;

rotationMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
T_e = rotationMatrix * thrust;
vdot = T_e / m - [0; g];

mdot = -alpha * thrust_mag;

xdot = [rdot; vdot; thetadot; wdot; mdot];

j_a = jacobian(xdot, x);
j_b = jacobian(xdot, u);
j_p = jacobian(xdot, p);

% Create equations of motion function for optimizer
matlabFunction(xdot,"File","Dynamics Models/3DoF/SymDynamics3DoF_mass_polar","Vars", [{t}; {x}; {u}; {L; I; alpha}]);

% Create equations of motion block for Simulink model
%matlabFunctionBlock('EoM_3DoF/SymDynamics3DoF',xdot,'Vars',[x; u; mass; L; I])

% Create Jacobian functions
matlabFunction(j_a,"File","Dynamics Models/3DoF/A_3DoF","Vars",[{t}; {x}; {u}; {L; I; alpha}]);
matlabFunction(j_b,"File","Dynamics Models/3DoF/B_3DoF","Vars",[{t}; {x}; {u}; {L; I; alpha}]);
matlabFunction(j_p,"File","Dynamics Models/3DoF/S_3DoF","Vars",[{t}; {x}; {u}; {L; I; alpha}]);