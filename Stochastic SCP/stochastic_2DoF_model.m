function [Y] = stochastic_2DoF_model(X, P)
%STOCHASTIC_MODEL Summary of this function goes here
%   Made to be able ot interface with UQLab
% Inputs: [N, N_in]
% - X(:, 1:5) initial estimated state stds
% - X(:, 6:10) initial state estimation stds
% - X(:, 11:15) final state stds
% - X(:, 16:18) disturbance stds
% - X(:, 19:23) measurement stds
% Parameters
% - stoch_prob_3DoF
% - ptr_ops
% Outputs: [N, N_out]
% - Y(:, ) Objective value (99% delta V both predicted and numerical)

%% Initialize
N_evals = size(X, 1);

Y = zeros([N_evals, 1]);

%% Extract Inputs
xhat0_stds = X(:, 1:5);
xtilde0_stds = X(:, 6:10);
x_f_stds = X(:, 11:15);
g_c_stds = X(:, 16:18);
g_0_stds = X(:, 19:23);

%% Optimize and Evaluate
for e = 1:P.thread_number:N_evals
    %% Adjust StochasticProblem with Inputs
    stoch_prob_2DoFs = {};
    ptr_opss = {};

    N_e = (e - 1) * P.thread_number;
    threads = min(N_evals - N_e, P.thread_number);

    for t = 1:threads
        stoch_prob_2DoF = P.stoch_prob;
    
        stoch_prob_2DoF.Phat0 = diag(xhat0_stds(N_e + t, :) .^ 2);
        stoch_prob_2DoF.Ptilde0 = diag(xtilde0_stds(N_e + t, :) .^ 2);
        stoch_prob_2DoF.Pf = diag(x_f_stds(N_e + t, :) .^ 2);
        stoch_prob_2DoF.cont.G = P.G(g_c_stds(N_e + t, :));
        stoch_prob_2DoF.filter.g_0 = P.g_0(g_0_stds(N_e + t, :));
        
        % Update the terminal constraints with the new r_y std
        stoch_prob_2DoF.xf = stoch_prob_2DoF.xf(x_f_stds(N_e + t, :));
        stoch_prob_2DoF.terminal_bc = stoch_prob_2DoF.terminal_bc(x_f_stds(N_e + t, :));
        stoch_prob_2DoF.nonconvex_constraints = stoch_prob_2DoF.nonconvex_constraints(x_f_stds(N_e + t, :));

        stoch_prob_2DoF = stoch_prob_2DoF.linearize();
        [stoch_prob_2DoF, Delta] = stoch_prob_2DoF.discretize(stoch_prob_2DoF.guess.x, stoch_prob_2DoF.guess.u, stoch_prob_2DoF.guess.p);

        stoch_prob_2DoFs{t} = stoch_prob_2DoF;
        ptr_opss{t} = P.ptr_ops;
    end

    %% Optimize
    stoch_ptr_sols = batch_stoch_ptr(stoch_prob_2DoFs, ptr_opss, threads);
    
    format longG
    save(sprintf("sensitivity_sol_2DoF_%s", string(round(posixtime(datetime)))), "stoch_ptr_sols");

    %% MC Simulations to Calculate Output
    for s = 1:numel(stoch_ptr_sols)
        if stoch_ptr_sols{s}{1}.converged
            Y(N_e + s, :) = MC_calc_output(stoch_prob_2DoF{s}, stoch_ptr_sols{s}{1}, P.m);
        else
            Y(N_e + s, :) = nan;
        end
    end
end

end

function [Y] = MC_calc_output(stoch_prob_2DoF, stoch_ptr_sol, m)
    K_k_opt = recover_gain_matrix(stoch_ptr_sol.X(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.S(:, :, stoch_ptr_sol.converged_i));

    t_ofb = zeros([stoch_prob_2DoF.N, m]);
    x_ofb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
    xhat_ofb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
    Phat_ofb = zeros([stoch_prob_2DoF.n.x, stoch_prob_2DoF.n.x, stoch_prob_2DoF.N, m]);
    u_ofb = zeros([stoch_prob_2DoF.n.u, stoch_prob_2DoF.Nu, m]);
    
    parfor i = 1:m
        [t_ofb(:, i), x_ofb(:, :, i), xhat_ofb(:, :, i), Phat_ofb(:, :, :, i), u_ofb(:, :, i)] = stoch_prob_2DoF.disc_prop(stoch_ptr_sol.x(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.u(:, :, stoch_ptr_sol.converged_i), stoch_ptr_sol.p(:, stoch_ptr_sol.converged_i), K_k_opt);
    end
    
    % Calculate 99% Delta V
    Y = zeros([1, 2]);

    objectives = stoch_ptr_sols.info.J;
    Y(1) = objectives(end);

    Y(2) = ksdensity(squeeze(sum(vecnorm(u_ofb, 2, 1) .* diff(t_ofb, 1, 2), 2)), 0.99, "Function","icdf"); 
end