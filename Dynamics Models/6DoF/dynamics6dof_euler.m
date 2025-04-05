% constants (same as 3DoF)
syms mass L I1 I2 I3 g; 
c = [g; mass; L; I1; I2; I3]

% states

r = sym("r", [3;1])
v = sym("v", [3,1])
theta = sym("theta", [3;1])
w = sym("w", [3;1])

x = [r; v; theta; w]

% controls

T = sym("T", [3;1])
gamma = sym("gamma", 1)

u = [T; gamma]

% calculate theta dot
b_inverse =  (1/sin(theta(2))) .* [0 sin(theta(3)) cos(theta(3));
     0 sin(theta(2))*cos(theta(3)) -sin(theta(2))*sin(theta(3));
     sin(theta(2)) -cos(theta(2))*sin(theta(3)) -cos(theta(2))*cos(theta(3))];
thetadot = b_inverse * w

% x dot
x_dot = v

Rxtheta1 = [1 0 0; 0 cos(theta(1)) -sin(theta(1)); 0 sin(theta(1)) cos(theta(1))];
Rxtheta2 = [cos(theta(2)) 0 sin(theta(2)); 0 1 0; -sin(theta(2)) 0 cos(theta(2))];
Rxtheta3 = [1 0 0; 0 cos(theta(3)) -sin(theta(3)); 0 sin(theta(3)) cos(theta(3))];
C_be = Rxtheta1 * Rxtheta2 * Rxtheta3;
C_eb = C_be.'
% v dot
firstPart = C_eb * T / mass * cos(u(4)) - [0, 0, g]
v_dot = C_eb * T / mass * cos(u(4)) - [0; 0; g]
%
M = norm(T) .* sin(u(4)) .* [1; 0; 0] + cross([-L; 0; 0], cos(u(4)) * T);
% w dot
w_dot = (M + (c([5; 6; 4]) - c([6; 4; 5])) .* w([2; 3; 1]) .* w([3; 1; 2]))

state = [x_dot; v_dot; thetadot; w_dot]
j_a = jacobian(state, x)
j_b = jacobian(state, u)


xdot_func6dof = matlabFunction(state);
j_a_func6dof = matlabFunction(j_a);
j_b_func6dof = matlabFunction(j_b);

save("dynamics6dof.mat", "xdot_func6dof", "j_a_func6dof", "j_b_func6dof")
save("SymDynamics6DoF.m", "xdot_func6dof")
save("SymUJacobian6DoF.m", "j_b_func6dof")
save("SymXJacobian6DoF.m", "j_a_func6dof")