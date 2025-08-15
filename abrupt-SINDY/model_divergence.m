function [diverged, lambda_bar, lambda_bar_measurement] = model_divergence(t_k, x_meas, x_pred, x_pred_t, A, delta_x_tol, alpha, x_err)
%MODEL_DIVERGENCE Detection of abrupt changes in model
%   Abrupt change in model is detected by comparing
%   the estimated Lyapunov time of the data with the model prediction
% t_k - time array
% x_meas - measured state trajectory
% x_pred - predicted state trajectory
% x_pred_t - predicted state trajectory starting from x_meas(:, 1)
% A - state jacobian function as a function of state (control, time, parameters plugged in)
% delta_x_tol - state fluctuation tolerance
% alpha - empirical measurement divergence multiplier accounting for finite-time statistics 
% x_err - @(x_1, x_2) function to compute error between two elements of x

% Calculate prediction horizon. When model changes, prediction horizon
% should decrease with smaller T indicating a bigger change
[T, T_ind] = prediction_horizon(t_k, x_meas, x_pred_t, delta_x_tol, x_err);
% Approximation of local Lyapunov exponent during T
lambda_bar = 0;
for k = 1 : T_ind
    lambda_bar = lambda_bar + eigs(A(x_meas(:, k)), 1);
end
lambda_bar = lambda_bar / T_ind; 

% Model and measurement have diverged if model time scale (lambda_bar) and 
% measured one differ
delta_x_T = x_err(x_pred(:, T_ind), x_meas(:, T_ind));
delta_x_t = x_err(x_pred(:, 1), x_meas(:, 1));
measurement_divergence = log(delta_x_T / delta_x_t); % paper not super clear here, probably typo??
% Measurement approximate local Lyapunov exponent
lambda_bar_measurement = measurement_divergence / T;

% Not sure on this inequality sign for non-chaotic systems with control
diverged = alpha * lambda_bar_measurement > lambda_bar;
end