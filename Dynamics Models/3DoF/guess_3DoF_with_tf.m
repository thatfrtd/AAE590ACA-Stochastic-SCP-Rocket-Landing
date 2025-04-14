function [x_guess, u_guess, tf_guess] = guess_3DoF_with_tf(x_0, x_f, steps, vehicle)
%STATE_GUESS Summary of this function goes here
%   Detailed explanation goes here

g = 9.81; % [m/s2]

% Control guess
T_guess = 0.4; % [N] Initial thrust guess
u_guess = [T_guess, 0] .* ones([steps, 2]);

% Parameter guess
v_0 = x_0(3:4);
v_f = x_f(3:4);
tf_guess = norm(v_f-v_0 + sqrt(2 * [0; g] * x_0(2)), 2)/(g - T_guess * vehicle.max_thrust/vehicle.m);

% State guess
x_guess = zeros([steps, 6]);

% r_guess
x_guess(:, 1:2) = straight_line_interpolate(x_0(1:2), x_f(1:2), steps);

% v_guess
v_cst = (x_f(1:2) - x_0(1:2)) / tf_guess;
x_guess(:, 3:4) = repmat(v_cst, steps, 1);

% θ_guess 
x_guess(:, 5) = straight_line_interpolate(x_0(5), x_f(5), steps);

% ω_guess
w_cst = (x_0(5) - x_f(5)) / tf_guess;

x_guess(:, 6) = w_cst;

end

