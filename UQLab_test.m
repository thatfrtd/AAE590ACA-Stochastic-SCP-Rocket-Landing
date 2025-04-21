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
Deterministic_3DoF

%% Initialize

G = @(g_c_stds) (@(t, x, u, p) [zeros([2, 3]); ... % velocity
                   [g_c_stds(1); g_c_stds(2)] .* eye([2, 3]); ... % acceleration
                   zeros([1, 3]); ... % angular velocity
                   g_c_stds(3) * [0, 0, 1]]); ... % angular acceleration
f_0 = @(t,x,u,p) x;
g_0 = @(g_0_stds) (@(t, x, u, p) diag(g_0_stds));

stoch_prob_3DoF = StochasticProblem.stochastify_discrete_problem(prob_3DoF, G(zeros([3, 1])), f_0, g_0(zeros([6, 1])), zeros(6), zeros(6), sol = ptr_sol);

K_k = -3*repmat([0, 1, 0, 0, 0, 0; 0, 0, 0, 0, -1, 0; 0, 1, 0, 0, 1, 0], 1, 1, stoch_prob_3DoF.Nu);

m = 20;

parameters.stoch_prob_3DoF = stoch_prob_3DoF;
parameters.G = G;
parameters.g_0 = g_0;
parameters.K_k = K_k;
parameters.m = m;
parameters.N_sub = 1;

%% Create Model
uqlab

Model1Opts.mFile = 'stochastic_3DoF_model_test';
Model1Opts.Parameters = parameters;
Model1Opts.isVectorized = false;
Model1Opts.isStochastic = false;
%modelopts.stochasticSim.supportRep = true; % Model can optimize then
%sample many realizations using MC (need to make isStochastic)

myModel = uq_createModel(Model1Opts);

%% Create Inpute Distributions
% Inputs: [N_in]
% - X(1:6) initial state stds
% - X(7:12) final state stds
% - X(13:15) disturbance stds
% - X(16:21) measurement stds

states = ["rx", "ry", "vx", "vy", "theta", "w"];
disturbances = ["ax", "ay", "alpha"];

for i_x0 = 1:6
    InputOpts.Marginals(i_x0).Name = states(i_x0) + "0";  % Initial state std
    InputOpts.Marginals(i_x0).Type = 'Uniform';
    InputOpts.Marginals(i_x0).Parameters = [1e-3 50e-3];  
end
for i_xf = 1:6
    InputOpts.Marginals(i_xf + i_x0).Name = states(i_xf) + "f";  % Final state std
    InputOpts.Marginals(i_xf + i_x0).Type = 'Uniform';
    InputOpts.Marginals(i_xf + i_x0).Parameters = [1e-4 5e-4];                                                                                                                                              
end
for i_dtsb = 1:3
    InputOpts.Marginals(i_dtsb + i_x0 + i_xf).Name = disturbances(i_dtsb);  % Disturbance std
    InputOpts.Marginals(i_dtsb + i_x0 + i_xf).Type = 'Uniform';
    InputOpts.Marginals(i_dtsb + i_x0 + i_xf).Parameters = [1e-4 1e-3];
end
for i_msr = 1:6
    InputOpts.Marginals(i_msr + i_x0 + i_xf + i_dtsb).Name = states(i_msr) + "m";  % Measurment std
    InputOpts.Marginals(i_msr + i_x0 + i_xf + i_dtsb).Type = 'Uniform';
    InputOpts.Marginals(i_msr + i_x0 + i_xf + i_dtsb).Parameters = [1e-3 5e-1] ; 
end

myInput = uq_createInput(InputOpts);

%%
X = uq_getSample(100);

Y = zeros(size(X, 1), 2);
% Manually do iterations so we can see what iteration it is on
for i = 1:size(X, 1)
    Y(i, :) = uq_evalModel(myModel, X(i, :));
    fprintf("Iteration %g out of %g Completed\n", i, size(X, 1))
end

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
