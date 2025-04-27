%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 17 April, 2025
% Description: Test of Monte Carlo and sensitivity analysis with UQLab.
% First optimize then do a number of MC runs to estimate the constraint
% satisfaction. 
% Most Recent Change: 19 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEED TO MAKE STOCHASTIC MODEL??
%% Run Deterministic 3DoF
Deterministic_2DoF_linear

nu = 2;
nr = 2;
nx = 5;

%% Stochastify Objective and Constraints
% Objective
stochastic_min_fuel_objective = @(x, u, p, X_k, S_k) einsum(@(k) norm(u(1:2, k), 2) + sigma_mag_confidence(1e-2, nu) * norm(S_k(:, (tri(k - 1, nx) + 1):tri(k, nx)), 2), 1:Nu) * (t_k(2) - t_k(1));

% Convex state path constraints
state_convex_constraints = {};

% Nonconvex state path constraints
glideslope_angle_max = deg2rad(80);
h_glideslope = @(sigma_xf) calculate_glideslope_offset(sigma_xf(1:2) * norminv(1 - 1e-3 / 2), glideslope_angle_max);
glideslope_constraint = @(sigma_xf) (@(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) (norm(x(1:2) + [0; h_glideslope(sigma_xf)] + 0*sigma_mag_confidence(1e-3, nr) * norm(X_k_ref(1:2, (tri(k - 1, nx) + 1):tri(k, nx)))) * cos(glideslope_angle_max) - (x(2) + h_glideslope(sigma_xf) - norminv(1 - 1e-3) * norm(X_k_ref(2, (tri(k - 1, nx) + 1):tri(k, nx))))));
state_nonconvex_constraints = @(sigma_xf) {glideslope_constraint(sigma_xf)};

% Convex control constraints
control_convex_constraints = {};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex control constraints
lcvx_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) norm(u(1:2)) + sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k) - thrust_magnitude_bound(S_k_ref, u_ref, k, t_k, T_max, m_0, alpha, nu);
min_thrust_constraint = @(x, u, p, X_k, S_k, x_ref, u_ref, p_ref, X_k_ref, S_k_ref, k) T_min * (exp(-(x_ref(5, k) - norminv(1 - 1e-3) * norm(X_k_ref(5, (tri(k - 1, nx) + 1):tri(k, nx)))))) - (norm(u_ref(1:2)) - sigma_mag_confidence(1e-3 / 2, nu) * norm(S_k));
control_nonconvex_constraints = {lcvx_thrust_constraint};

% Combine nonconvex constraints
nonconvex_constraints = @(sigma_xf) [state_nonconvex_constraints(sigma_xf), control_nonconvex_constraints];

%% Initialize
delta_t = 0.1;

G = @(g_c_stds) (@(t, x, u, p) sqrt(delta_t) * [zeros([2, 3]); ... % velocity
                   [g_c_stds(1); g_c_stds(2)] .* eye([2, 3]); ... % acceleration
                   g_c_stds(3) * [0, 0, 1]]); ... % mass flow 
f_0 = @(t,x,u,p) x;
g_0 = @(g_0_stds) (@(t, x, u, p) diag(g_0_stds));


f_stoch = @(t, x, u, p) SymDynamics2DoF_linear_noumag(t, x, u, vehicle.alpha);
ptr_sol_mod = ptr_sol;
ptr_sol_mod.u = ptr_sol_mod.u(1:2, :, :);
prob_2DoF.xf = @(sigma_xf) [0; sigma_xf(2) * norminv(1 - 1e-3); 0; 0];
stoch_terminal_bc = @(sigma_xf) (@(x, p) [x(1:4) - prob_2DoF.xf(sigma_xf); 0]);
stoch_prob_2DoF = StochasticProblem.stochastify_discrete_problem(prob_2DoF, G(zeros([3, 1])), f_0, g_0(zeros([6, 1])), zeros(6), zeros(6), zeros(6), f = f_stoch, sol = ptr_sol_mod, objective = stochastic_min_fuel_objective, convex_constraints = convex_constraints, nonconvex_constraints = nonconvex_constraints, terminal_bc = stoch_terminal_bc);

ptr_ops.iter_min = 3;
ptr_ops.iter_max = 10;

parameters.stoch_prob = stoch_prob_2DoF;
parameters.ptr_ops = ptr_ops;
parameters.G = G;
parameters.g_0 = g_0;
parameters.m = 200;
parameters.N_sub = 1;
parameters.thread_number = 4;

%% Create Model
uqlab

Model1Opts.mFile = 'stochastic_2DoF_model';
Model1Opts.Parameters = parameters;
Model1Opts.isVectorized = true;
Model1Opts.isStochastic = false;
%modelopts.stochasticSim.supportRep = true; % Model can optimize then
%sample many realizations using MC (need to make isStochastic)

myModel = uq_createModel(Model1Opts);

%% Create Input Distributions
% Inputs: [N_in]
% - X(:, 1:5) initial estimated state stds
% - X(:, 6:10) initial state estimation stds
% - X(:, 11:15) final state stds
% - X(:, 16:18) disturbance stds
% - X(:, 19:23) measurement stds

% Nominal uncertainties
% Initial estimated state
sigma_xhat0 = [10e-3, 100e-3; ... % r_x
            10e-3, 100e-3; ... % r_y
            1e-3, 50e-3; ... % v_x
            1e-3, 50e-3; ... % v_y
            1e-4, 1e-4]; % mass

% Initial state estimation error
sigma_xtilde0 = [1e-3, 50e-3; ... % r_x
            1e-3, 50e-3; ... % r_y
            1e-3, 10e-3; ... % v_x
            1e-3, 50e-3; ... % v_y
            2e-5, 2e-5]; % mass

% Final state
sigma_xf = [1e-3, 1e-2; ... % r_x
            1e-3, 5e-3; ... % r_y
            1e-4, 1e-3; ... % v_x
            1e-4, 1e-3; ... % v_y
            1e-4, 1e-4]; % mass

% Disturbance
sigma_accelx = 0.5e-3;
sigma_accely = 0.1e-3;
sigma_m = 1e-7;
sigma_dtsb = [0.5e-3, 1e-3; ... % a_x
              0.1e-3, 1e-3; ... % a_y
              1e-7, 1e-7]; ... % dmass

% Measurement functions
g_0_stds = [0.1e-3, 1e-3; ... % r_x
            0.5e-3, 2e-3; ... % 
            0.05e-3, 0.5e-3; ... % v_x
            0.05e-3, 0.5e-3; ... % v_y
            1e-5, 1e-5]; ... % mass

% Create input distibutions
states = ["rx", "ry", "vx", "vy", "mass"];
disturbances = ["ax", "ay", "dmass"];

for i_xhat0 = 1:nx
    InputOpts.Marginals(i_xhat0).Name = states(i_xhat0) + "hat0";  % Initial state std
    InputOpts.Marginals(i_xhat0).Type = 'Uniform';
    InputOpts.Marginals(i_xhat0).Parameters = sigma_xhat0(i_xhat0, :);  
end
for i_xtilde0 = 1:nx
    InputOpts.Marginals(i_xtilde0 + i_xhat0).Name = states(i_xtilde0) + "tilde0";  % Initial state std
    InputOpts.Marginals(i_xtilde0 + i_xhat0).Type = 'Uniform';
    InputOpts.Marginals(i_xtilde0 + i_xhat0).Parameters = sigma_xtilde0(i_xtilde0, :);  
end
for i_xf = 1:nx
    InputOpts.Marginals(i_xf + i_xhat0 + i_xtilde0).Name = states(i_xf) + "f";  % Final state std
    InputOpts.Marginals(i_xf + i_xhat0 + i_xtilde0).Type = 'Uniform';
    InputOpts.Marginals(i_xf + i_xhat0 + i_xtilde0).Parameters = sigma_xf(i_xf, :);                                                                                                                                              
end
for i_dtsb = 1:3
    InputOpts.Marginals(i_dtsb + i_xhat0 + i_xtilde0 + i_xf).Name = disturbances(i_dtsb);  % Disturbance std
    InputOpts.Marginals(i_dtsb + i_xhat0 + i_xtilde0 + i_xf).Type = 'Uniform';
    InputOpts.Marginals(i_dtsb + i_xhat0 + i_xtilde0 + i_xf).Parameters = sigma_dtsb(i_dtsb, :);
end
for i_msr = 1:nx
    InputOpts.Marginals(i_msr + i_xhat0 + i_xtilde0 + i_xf + i_dtsb).Name = states(i_msr) + "m";  % Measurment std
    InputOpts.Marginals(i_msr + i_xhat0 + i_xtilde0 + i_xf + i_dtsb).Type = 'Uniform';
    InputOpts.Marginals(i_msr + i_xhat0 + i_xtilde0 + i_xf + i_dtsb).Parameters = g_0_stds(i_msr, :); 
end

myInput = uq_createInput(InputOpts);

%%
X = uq_getSample(200);

% Manually do iterations so we can see what iteration it is on and handle
% when the optimization doesn't converge
Y = uq_evalModel(myModel, X);

format longG
save(sprintf("sensitivity_99p_dv_2DoF_%s", string(round(posixtime(datetime)))), "Y");

%%
scatter(X(:, 14), Y(:, 2))

%% Fit Polynomial Chaos Expansion Surrogate Model
PCEOpts.Type = 'Metamodel';
PCEOpts.MetaType = 'PCE';
%PCEOpts.FullModel = myModel;

PCEOpts.Degree = 1:15;

MetaOpts.Method = 'LARS';

%PCEOpts.TruncOptions.qNorm = 0.7;
%PCEOpts.TruncOptions.MaxInteraction = 2;

% If wanting UQLab to run iterations for us (no loading bar...)
%PCEOpts.ExpDesign.NSamples = 200;
%PCEOpts.ExpDesign.Sampling = 'LHS';

% Can also specify results from existing data
PCEOpts.ExpDesign.Sampling = "user";
PCEOpts.ExpDesign.X = X;
PCEOpts.ExpDesign.Y = Y;
% Or from a data file
%PCEOpts.DataFile = "data.mat"; % Has to have only .X and .Y as variables

myPCE = uq_createModel(PCEOpts);

uq_display(myPCE, 1:2)

%% Look at Surrogate Model Error Visually

Y_test = uq_evalModel(myPCE, myPCE.ExpDesign.X);

% Percent error
scatter(myPCE.ExpDesign.X(:, 3), 100 * (myPCE.ExpDesign.Y(:, 1) - Y_test(:, 1)) ./ myPCE.ExpDesign.Y(:, 1)); hold on
scatter(myPCE.ExpDesign.X(:, 3), 100 * (myPCE.ExpDesign.Y(:, 2) - Y_test(:, 2)) ./ myPCE.ExpDesign.Y(:, 2)); hold off

%% Sobol Indices Sensitivity Analysis
% Set up analysis
SobolOpts.Type = 'Sensitivity';
SobolOpts.Method = 'Sobol';
SobolOpts.Sobol.Order = 1;

% Run analysis
mySobolAnalysisPCE = uq_createAnalysis(SobolOpts);

%% Display Sensitivity Analysis Results
uq_display(mySobolAnalysisPCE)
%%
uq_figure()
uq_bar(1:21, mySobolAnalysisPCE.Results.Total, 0.25)
% Set axes limits
ylim([0 1])
% Set labels
xlabel('Variable name')
ylabel('Total Sobol'' indices')
set(gca, 'XTick', 1:length(InputOpts.Marginals),...
    'XTickLabel', mySobolAnalysisPCE.Results.VariableNames)
% Set Title
title(sprintf('PCE-based Sobol Indice Sensitivity Analysis (%d simulations)', myPCE.ExpDesign.NSamples))
% Set legend
uq_legend({...
    "Final $\sigma_{r_x}$", "Final $\sigma_{r_y}$"},...
    'Location', 'northeast', "Interpreter", "latex")
