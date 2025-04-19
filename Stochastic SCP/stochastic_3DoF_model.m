function [Y] = stochastic_3Dof_model(X, P)
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
stoch_prob_3DoF.P0 = diag(x_0_stds .^ 2);
stoch_prob_3DoF.Pf = diag(x_f_stds .^ 2);
stoch_prob_3DoF.G = ;
stoch_prob_3DoF.g_0 = ;

%% Optimize
ptr_sol = Stochastic_ptr(prob_3DoF, ptr_ops);

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

