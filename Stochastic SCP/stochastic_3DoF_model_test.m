function [Y] = stochastic_3DoF_model_test(X, P)
%STOCHASTIC_MODEL Summary of this function goes here
%   Made to be able ot interface with UQLab (not vectorized)
% Inputs: [N_in]
% - X(1:6) initial state stds
% - X(7:12) final state stds
% - X(13:15) disturbance stds
% - X(16:21) measurement stds
% Parameters
% - stoch_prob_3DoF
% - G
% - g_0
% - K_k
% - m
% Outputs: [N_out]
% - Y(1, 1:2) Final covariance

%% Extract Inputs
x_0_stds = X(1:6);
x_f_stds = X(7:12);
g_c_stds = X(13:15);
g_0_stds = X(16:21);

%% Adjust StochasticProblem with Inputs
stoch_prob_3DoF = P.stoch_prob_3DoF;

stoch_prob_3DoF.P0 = diag(x_0_stds .^ 2);
stoch_prob_3DoF.Pf = diag(x_f_stds .^ 2);
stoch_prob_3DoF.cont.G = P.G(g_c_stds);
stoch_prob_3DoF.filter.g_0 = P.g_0(g_0_stds);

stoch_prob_3DoF = stoch_prob_3DoF.linearize();
[stoch_prob_3DoF, Delta] = stoch_prob_3DoF.discretize(stoch_prob_3DoF.guess.x, stoch_prob_3DoF.guess.u, stoch_prob_3DoF.guess.p);

K_k = P.K_k;

%% Optimize - Need to implement first :))
%ptr_sol = Stochastic_ptr(prob_3DoF, ptr_ops);

%% MC Simulations
m = P.m;

%t_fb = zeros([P.N_sub * (stoch_prob_3DoF.N - 1) + 1, m]);
x_fb = zeros([stoch_prob_3DoF.n.x, P.N_sub * (stoch_prob_3DoF.N - 1) + 1, m]);
%xhat_fb = zeros([stoch_prob_3DoF.n.x, P.N_sub * (stoch_prob_3DoF.N - 1) + 1, m]);
%Phat_fb = zeros([stoch_prob_3DoF.n.x, stoch_prob_3DoF.n.x, P.N_sub * (stoch_prob_3DoF.N - 1) + 1, m]);
%u_fb = zeros([stoch_prob_3DoF.n.u, P.N_sub * (stoch_prob_3DoF.N - 1) + 1, m]);

parfor i = 1:m
    %[t_fb(:, i), x_fb(:, :, i), xhat_fb(:, :, i), Phat_fb(:, :, :, i), u_fb(:, :, i)] = stoch_prob_3DoF.disc_prop(stoch_prob_3DoF.guess.x, stoch_prob_3DoF.guess.u, stoch_prob_3DoF.guess.p, K_k);
    [~, x_fb(:, :, i), ~, ~, ~] = stoch_prob_3DoF.disc_prop(stoch_prob_3DoF.guess.x, stoch_prob_3DoF.guess.u, stoch_prob_3DoF.guess.p, K_k);
end

%% Calculate Outputs
% Calculate constraint margins?

% Estimate final covariance and standard deviation for r_x and r_y
final_r_stds = sqrt(eigs(cov(squeeze(x_fb(1:2,end,:))' - mean(squeeze(x_fb(1:2,end,:)), 2)')));

%% Package Outputs
Y = reshape(final_r_stds, [1, 2]);

end

