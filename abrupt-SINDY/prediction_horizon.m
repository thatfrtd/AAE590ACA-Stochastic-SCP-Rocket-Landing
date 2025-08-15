function [T, T_ind] = prediction_horizon(t_k, x_meas, x_pred_t, delta_x_tol, x_err)
%PREDICTION_HORIZON Summary of this function goes here
%   The first passage time dt where prediction x_pred(t + dt) with 
%   initial condition x_meas(t) and measurement x_meas(t + dt) differ by
%   more than delta_x. Note t is assumed to correspond to the first element
%   of the arrays (t = t_k(1)) so x_meas(t) = x_meas(:, 1) = x_pred_t(:, 1).
% t_k - time array
% x_meas - measured state trajectory
% x_pred_t - predicted state trajectory starting from x_meas(:, 1)
% delta_x_tol - state fluctuation tolerance
% x_err - function to compute error between two elements of x

T_ind = find(x_err(x_pred_t, x_meas) >= delta_x_tol, 1, "first");
T = t_k(T_ind) - t_k(1);

end

