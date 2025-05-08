function [guess] = guess_6DoF(x_0, x_f, N, Nu, t_step, vehicle)
%STATE_GUESS Summary of this function goes here
%   Detailed explanation goes here

g = 9.81; % [m/s2]

% Control guess
T_guess = 0.4; % [N] Initial thrust guess
u_guess = [T_guess; 0; 0; T_guess; 0] .* ones([5, Nu]);

% Parameter guess
v_0 = x_0(4:6);
v_f = x_f(4:6);
%tf_guess = norm(v_f-v_0 + sqrt(2 * [0; g] * x_0(2)), 2)/(g - T_guess * vehicle.max_thrust/vehicle.m);

% State guess
x_guess = zeros([13, N]);

% r_guess
x_guess(1:3, :) = straight_line_interpolate(x_0(1:3), x_f(1:3), N);

% v_guess
v_cst = (x_f(1:3) - x_0(1:3)) / (N * t_step);
x_guess(4:6, :) = repmat(v_cst, 1, N);

% θ_guess 

q_start = quaternion(x_0(7:10).');
q_end = quaternion(x_f(7:10).');
quat2 = slerp(q_start, q_end, linspace(0, 1, N));
x_guess(7:10, :) = compact(quat2).';

% ω_guess
for i = 1:N
    x_guess(11:13, i) = angvel(quat2(i), t_step, "frame").';
end

guess.x = x_guess;
guess.u = u_guess;
guess.p = [];

end