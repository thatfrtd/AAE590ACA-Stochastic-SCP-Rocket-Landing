function [A_k, B_k, c_k] = integrate_discrete_ZOH(x0, A, B, c, f, u, tspan, tolerances)
    % Integrates STM and state with Bk and ck
    %   Uses ODE45 to integrate the State Transition Matrix and the state using
    %   the given A matrix and dynamics f over the time period in tspan using
    %   the specified error tolerance.

    % Create initial condition
    nx = numel(x0);

    STM0 = eye(nx);
    B0 = zeros(size(B(0, x0, u(tspan(1)))));
    c0 = zeros([nx, 1]);

    nu = size(B0, 2);

    y0 = [x0; STM0(:); B0(:); c0];

    % Simulate    
    [~, y] = ode45(@(t, y) STM_diff_eq_ZOH(t, y, A, B, c, f, u, nx), tspan, y0, tolerances);

    y_f = y(end, :); 

    % Unpack solution
    x = y_f(:, 1:nx);
    A_k = reshape(y_f(:, (nx + 1) : (nx * (nx + 1))), nx, nx);
    B_k = A_k * reshape(y_f(:, (nx * (nx + 1) + 1) : (nx * (nx + 1) + nx * nu)), nx, nu);
    c_k = A_k * y_f(:, (nx * (nx + 1) + nx * nu + 1) : (nx * (nx + 1) + nx * nu + nx)).';
end

function [ydot] = STM_diff_eq_ZOH(t, y, A, B, c, f, u, n)
    x = y(1:n);
    STM = reshape(y((n + 1) : (n * (n + 1))), n, n);

    xdot = f(t, x, u(t));
    A_kdot = A(t, x, u(t)) * STM;
    B_kdot = STM \ B(t, x, u(t));
    c_kdot = STM \ c(t, x, u(t));

    ydot = [xdot; A_kdot(:); B_kdot(:); c_kdot];
end