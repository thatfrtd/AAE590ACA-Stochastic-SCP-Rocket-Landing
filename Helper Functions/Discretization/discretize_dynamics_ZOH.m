function [A_k, B_k, c_k] = discretize_dynamics_ZOH(f, A, B, c, N, tspan, x_ref, u_ref, tolerances)
    % Discretization of a dynamical system assuming ZOH control'''
    % Make c optional? c = f - Ax - Bu

    nx = numel(x_ref(tspan(1)));
    nu = numel(u_ref(tspan(1)));

    t_k = linspace(tspan(1), tspan(2), N + 1);
    A_k = zeros([nx, nx, N]);
    B_k = zeros([nx, nu, N]);
    c_k = zeros([nx, 1, N]);
    
    for k = 1:N
        [A_k(:, :, k), B_k(:, :, k), c_k(:, :, k)] = integrate_discrete_ZOH(x_ref(t_k(k)), A, B, c, f, u_ref, [t_k(k), t_k(k + 1)], tolerances);
    end
end