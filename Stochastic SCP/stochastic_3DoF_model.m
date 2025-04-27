function [Y] = stochastic_3DoF_model(X, P)
%STOCHASTIC_MODEL Summary of this function goes here
%   Made to be able ot interface with UQLab
% Inputs: [N, N_in]
% - X(:, 1:6) initial state stds
% - X(:, 7:12) final state stds
% - X(:, 13:15) disturbance stds
% - X(:, 16:21) measurement stds
% Parameters
% - stoch_prob_3DoF
% - ptr_ops
% Outputs: [N, N_out]
% - Y(:, ) Constraint margins?

%% Extract Inputs
x_0_stds = X(:, 1:6);
x_f_stds = X(:, 7:12);
g_c_stds = X(:, 13:15);
g_0_stds = X(:, 16:21);

%% Adjust StochasticProblem with Inputs
stoch_prob_3DoF = P.stoch_prob_3DoF;

stoch_prob_3DoF.P0 = diag(x_0_stds .^ 2);
stoch_prob_3DoF.Pf = diag(x_f_stds .^ 2);
stoch_prob_3DoF.cont.G = P.G(g_c_stds);
stoch_prob_3DoF.filter.g_0 = P.g_0(g_0_stds);

stoch_prob_3DoF = stoch_prob_3DoF.linearize();
[stoch_prob_3DoF, Delta] = stoch_prob_3DoF.discretize(stoch_prob_3DoF.guess.x, stoch_prob_3DoF.guess.u, stoch_prob_3DoF.guess.p);

%% Optimize
stoch_ptr_sols = batch_stoch_ptr(stoch_prob_2DoF, ptr_ops, P.thread_number);

%% MC Simulations
m = size();

t = zeros([numel(tspan), m]);
x = zeros([prob_3DoF.n.x, numel(tspan), m]);
u = zeros([prob_3DoF.n.u, numel(tspan), m]);

parfor i = 1:m
    [t(:, i), x(:, :, i), u_applied(:, :, i)] = sode45(stoch_prob_3DoF.cont.f, G, u_func, guess.p, w, tspan, delta_t, x_0_m(:, i), tolerances);
end

%% Package Outputs
% Calculate constraint margins?

end

