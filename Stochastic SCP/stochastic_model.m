function [outputArg1,outputArg2] = stochastic_model()
%STOCHASTIC_MODEL Summary of this function goes here
%   Made to be able ot interface with UQLab

m = size();

t = zeros([numel(tspan), m]);
x = zeros([prob_3DoF.n.x, numel(tspan), m]);
u = zeros([prob_3DoF.n.u, numel(tspan), m]);

parfor i = 1:m
    [t(:, i), x(:, :, i), u_applied(:, :, i)] = sode45(prob_3DoF.cont.f, G, u_func, guess.p, w, tspan, delta_t, x_0_m(:, i), tolerances);
end
end

