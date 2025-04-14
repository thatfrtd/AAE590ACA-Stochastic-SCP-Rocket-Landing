function [t_sub, x_array, xhat_array, Phat_array] = propagate_cont_kalman_filter(x_0, P_0, u, f, G, A_ref, B_ref, c_ref, G_ref, L_k, C_k, D_k, f_0, g_0, t_k, tspan, N_sub, w_k, v_k, tolerances)
%CONT_KALMAN_FILTER Continuous propagation of continuous-discrete Kalman
%filter
%   Note - control is constant between measurements

nx = numel(x_0);

t_sub = linspace(t_k(1), t_k(end), N_sub * numel(t_k));
if t_k(end) ~= tspan(end)
    t_sub = [t_sub(1:(end - 1)), linspace(t_k(end), tspan(end), max(ceil((t_k(end) - tspan(end)) / (t_sub(2) - t_sub(1))), 3))];
end

x_array = zeros([nx, numel(t_sub)]);
x_array(:, 1) = x_0;
xhat_array = zeros([nx, numel(t_sub)]);
xhat_array(:, 1) = x_0;
Phat_array = zeros([nx, nx, numel(t_sub)]);
Phat_array(:, :, 1) = P_0;

for k = 1:numel(t_k)
    if k ~= t_k || t_k(end) == tspan(end)
        k_0 = max((k - 1) * N_sub - 1, 1);
        y_k = [reshape(x_0(:, k_0), []); reshape(xhat_array(:, k_0), []); reshape(Phat_array(:, :, k_0), [])];
        u_k = u(t_k(k), xhat_array(:, k_0));
        % Continuously propogate over the interval between measurements (with constant control)
        [t_m, y_m] = ode45(@(t, y) cont_kalman_filter_DE(t, y, u_k, f, G, A_ref, B_ref, c_ref, G_ref, w_k, nx), t_sub((k - 1) * N_sub + (1:N_sub)), y_k, tolerances);
    
        x_m = y_m(:, 1:nx)';
        xhat_minus = y_m(:, nx + (1:nx))';
        Phat_minus = permute(reshape(y_m(:, (n+1):end), [nx, nx]), [2, 3, 1]);
    
        x_array(:, (k - 1) * N_sub + (1:N_sub)) = x_m;
        xhat_array(:, (k - 1) * N_sub + (1:N_sub)) = xhat_minus;
        Phat_array(:, :, (k - 1) * N_sub + (1:N_sub)) = Phat_minus;
    
        % Perform discrete measurement update
        y_k = f_0(t, x_m, u(t, xhat_minus(:, end))) + g_0(x_m, u(t, xhat_minus(:, end))) * v_k;
        ytilde_minus_k = innovation_process(y_k, C_k(k), xhat_minus(:, end));
        xhat_array(:, k * N_sub) = estimate_measurement_update(xhat_minus(:, end), L_k(k), ytilde_minus_k);
        Phat_array(:, :, k * N_sub) = covariance_measurement_update(L_k(k), C_k(k), Phat_minus(:, :, end), D_k(k));
    else
        [t_m, y_m] = ode45(@(t, y) cont_kalman_filter_DE(t, y, u, f, G, A_ref, B_ref, c_ref, G_ref, w_k, nx), t_sub(((k - 1) * N_sub + 1):end), y_0, tolerances);

        x_m = y_m(:, 1:nx)';
        xhat_minus = y_m(:, nx + (1:nx))';
        Phat_minus = permute(reshape(y_m(:, (n+1):end), [nx, nx]), [2, 3, 1]);
    
        x_array(:, ((k - 1) * N_sub + 1):end) = x_m;
        xhat_array(:, ((k - 1) * N_sub + 1):end) = xhat_minus;
        Phat_array(:, :, ((k - 1) * N_sub + 1):end) = Phat_minus; 
    end
end
end

function [ydot] = cont_kalman_filter_DE(t, y, u, f, G, A_ref, B_ref, c_ref, G_ref, w_k, nx)
    x = y(1:nx);
    xhat = y(nx + (1:nx));
    Phat = reshape(y((n+1):end), [nx, nx]);

    xdot = f(t, x, u(t, xhat)) + G(t, x, u(t, xhat)) * w_k;
    xhatdot = A_ref(t) * xhat + B_ref(t) * u(t, xhat) + c_ref(t);
    Phatdot = continuous_kalman_covariance_derivative(A_ref(t), Phat, G_ref(t));

    ydot = [xdot(:); xhatdot(:); Phatdot(:)];
end