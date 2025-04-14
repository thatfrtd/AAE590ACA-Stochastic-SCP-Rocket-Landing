function [A_k, B_k, E_k, c_k, G_k, x_kp1] = integrate_stochastic_discrete_ZOH(x0, A, B, E, c, f, G, u, p, tspan, tolerances)
    % Integrates STM and state with Bk and ck
    %   Uses ODE45 to integrate the State Transition Matrix and the state using
    %   the given A matrix and dynamics f over the time period in tspan using
    %   the specified error tolerance.

    % Create initial condition
    nx = numel(x0);
    np = numel(p);
    nw = size(G(0, x0, 0, 0), 2);

    STM0 = eye(nx);
    B0 = zeros(size(B(0, x0, u, p)));
    E0 = zeros(nx, np);
    c0 = zeros([nx, 1]);
    G0 = zeros([nx, nw]); % Not 100% sure on dimension yet

    nu = size(B0, 2);
    
    y0 = [x0; STM0(:); B0(:); E0(:); c0, G0(:)];

    % Simulate    
    [~, y] = ode45(@(t, y) STM_diff_eq_ZOH(t, y, A, B, E, c, f, G, u, p, nx), tspan, y0, tolerances);

    y_f = y(end, :); 

    % Unpack solution
    x_kp1 = y_f(:, 1:nx)';
    A_k = reshape(y_f(:, (nx + 1) : (nx * (nx + 1))), nx, nx);
    B_k = A_k * reshape(y_f(:, (nx * (nx + 1) + 1) : (nx * (nx + 1) + nx * nu)), nx, nu);
    E_k = A_k * reshape(y_f(:, (nx * (nx + 1) + nx * nu + 1) : (nx * (nx + 1) + nx * nu + nx * np)), nx, np);
    c_k = A_k * y_f(:, (nx * (nx + 1) + nx * nu + nx * np + 1) : (nx * (nx + 1) + nx * nu + nx * np + nx)).';
    Sigma_k = A_k * reshape(y_f(:, (nx * (nx + 1) + nx * nu + nx * np + nx) + 1 : (nx * (nx + 1) + nx * nu + nx * np + nx) + nx * nw), nx, nw) * Ak'; % If not working, check this dimension first
    G_k = chol(Sigma_k, 'lower');
end

function [ydot] = STM_diff_eq_ZOH(t, y, A, B, E, c, f, G, u, p, n)
    x = y(1:n);
    STM = reshape(y((n + 1) : (n * (n + 1))), n, n);

    xdot = f(t, x, u, p);
    A_kdot = A(t, x, u, p) * STM;
    B_kdot = STM \ B(t, x, u, p);
    E_kdot = STM \ E(t, x, u(t), p);
    c_kdot = STM \ c(t, x, u, p);
    G_dot = STM\G(t, x, u, p) * (STM \ G(t, x, u, p))';

    ydot = [xdot; A_kdot(:); B_kdot(:); E_kdot(:); c_kdot, G_dot(:)];
end