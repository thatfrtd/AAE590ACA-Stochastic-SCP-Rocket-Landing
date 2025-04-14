function [guess] = guess_3DoF(x_0, x_f, N, Nu, t_step, vehicle)
%STATE_GUESS Summary of this function goes here
%   Detailed explanation goes here

g = 9.81; % [m/s2]

% Control guess
T_guess = vehicle.min_thrust; % [N] Initial thrust guess
u_guess = [T_guess; 0; T_guess] .* ones([3, Nu]);

% Parameter guess
v_0 = x_0(3:4);
v_f = x_f(3:4);

% State guess
x_guess = zeros([6, N]);

% r_guess
x_guess(1:2, :) = straight_line_interpolate(x_0(1:2), x_f(1:2), N);

% v_guess
v_cst = (x_f(1:2) - x_0(1:2)) / (N * t_step);
x_guess(3:4, :) = repmat(v_cst, 1, N);

% θ_guess 
x_guess(5, :) = straight_line_interpolate(x_0(5), x_f(5), N);

% ω_guess
w_cst = (x_0(5) - x_f(5)) / (N * t_step);

x_guess(6, :) = w_cst;

guess.x = x_guess;
guess.u = u_guess;
guess.p = [];

end

