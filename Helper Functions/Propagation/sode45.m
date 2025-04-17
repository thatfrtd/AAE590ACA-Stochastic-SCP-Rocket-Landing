function [t, x, u_applied] = sode45(f, G, u, p, w, tspan, delta_t, x0, tolerances)
%SODE45 Summary of this function goes here
%   Detailed explanation goes here

t_k = tspan(1):delta_t:tspan(end);
w_k = w(numel(t_k));

w_k_func = @(t) interp1(t_k, w_k', t, "previous", "extrap")';

[t, x] = ode45(@(t, x) f(t, x, u(t, x), p) + G(t, x, u(t, x), p) / sqrt(delta_t) * w_k_func(t), tspan, x0, tolerances);

x = x';

u_applied = zeros(size(u(0, x0), 1), numel(t));
for i = 1:numel(t)
    u_applied(:, i) = u(t(i), x(:, i));
end
end

