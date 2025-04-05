function [A_k, B_k_plus, B_k_minus, c_k] = integrate_discrete_FOH(x0, A, B, c, f, u, tspan, tolerances)
    % Integrates STM and state with Bk+, Bk-, and ck
    %   Uses ODE45 to integrate the State Transition Matrix and the state using
    %   the given A matrix and dynamics f over the time period in tspan using
    %   the specified error tolerance.

    % Create initial condition
    nx = numel(x0);

    STM0 = eye(nx);
    B0 = zeros(size(B(0, x0, u(0))));
    c0 = zeros([nx, 1]);

    nu = size(B0, 2);

    y0 = [x0; STM0(:); B0(:); B0(:); c0];

    % Set up linear interpolation functions
    sigma_plus = @(t) (t - tspan(1)) ./ (tspan(2) -  tspan(1));
    sigma_minus = @(t) (tspan(2) - t) ./ (tspan(2) -  tspan(1));

    % Simulate    
    [~, y] = ode45(@(t, y) STM_diff_eq_FOH(t, y, A, B, c, f, u, sigma_plus, sigma_minus, nx), tspan, y0, tolerances);

    y_f = y(end, :);

    % Unpack solution
    x = y_f(:, 1:nx);
    A_k = reshape(y_f(:, (nx + 1) : (nx * (nx + 1))), nx, nx);
    B_k_plus = A_k * reshape(y_f(:, (nx * (nx + 1) + 1) : (nx * (nx + 1) + nx * nu)), nx, nu);
    B_k_minus = A_k * reshape(y_f(:, (nx * (nx + 1) + nx * nu + 1) : (nx * (nx + 1) + 2 * nx * nu)), nx, nu);
    c_k = A_k * y_f(:, (nx * (nx + 1) + 2 * nx * nu + 1) : (nx * (nx + 1) + 2 * nx * nu + nx)).';
end

function ydot = STM_diff_eq_FOH(t, y, A, B, c, f, u, sigma_plus, sigma_minus, n)
    x = y(1:n);
    STM = reshape(y((n + 1) : (n * (n + 1))), n, n);

    xdot = f(t, x, u(t));
    A_kdot = A(t, x, u(t)) * STM;
    B_k_plusdot = STM \ B(t, x, u(t)) * sigma_plus(t);
    B_k_minusdot = STM \ B(t, x, u(t)) * sigma_minus(t);
    c_kdot = STM \ c(t, x, u(t));

    ydot = [xdot; A_kdot(:); B_k_plusdot(:); B_k_minusdot(:); c_kdot];
end