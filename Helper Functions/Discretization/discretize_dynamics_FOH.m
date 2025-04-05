function [A_k, B_k_plus, B_k_minus, c_k] = discretize_dynamics_FOH(f, A, B, c, N, tspan, x_ref, u_ref, tolerances)
    % Discretization of a dynamical system assuming FOH control
    % Make c optional? c = f - Ax - Bu

    nx = numel(x_ref(0));
    nu = numel(u_ref(0));

    t_k = linspace(tspan(1), tspan(2), N + 1);
    A_k = zeros([nx, nx, N]);
    B_k_plus = zeros([nx, nu, N]);
    B_k_minus = zeros([nx, nu, N]);
    c_k = zeros([nx, 1, N]);
    
    for k = 1:N
        [A_k(:, :, k), B_k_plus(:, :, k), B_k_minus(:, :, k), c_k(:, :, k)] = integrate_discrete_FOH(x_ref(t_k(k)), A, B, c, f, u_ref, [t_k(k), t_k(k + 1)], tolerances);
    end
end