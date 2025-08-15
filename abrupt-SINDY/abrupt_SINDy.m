function [Xi, residual] = abrupt_SINDy(x_meas, numer_deriv, Xi0, Theta, lambda, max_iter_sparse, max_iter_abrupt, residual_tol)
%ABRUPT_SINDY Sparse Identification of Nonlinear Dynamics for Rapid Model Recovery
%   Based on arXiv:1803.00894v2
%   Excerpt from abstract: In this work, we propose
% a conceptual framework to recover parsimonious models of a system in response to
% abrupt changes in the low-data limit. First, the abrupt change is detected by comparing
% the estimated Lyapunov time of the data with the model prediction. Next, we apply the
% sparse identification of nonlinear dynamics (SINDy) regression to update a previously identified
% model with the fewest changes, either by addition, deletion, or modification of existing
% model terms
arguments
    x_meas % measured state trajectory
    numer_deriv % @(x) numerical derivative function
    Xi0 % current coefficients for model
    Theta % candidate library which are evaluated at x_meas and u
    lambda % minimum coefficient magnitude for sparse regression
    max_iter_sparse = 10 % maximum iterations for sequentially thresholded least squares regression
    max_iter_abrupt = 10 % maximum iterations for rapid model discovery procedure
    residual_tol = 1e-3 % updated model residual tolerance
end

%% Detect Model Divergence - TRIGGER abrupt_SINDy UPDATE AFTER model_divergence DETECTS DIVERGENCE - put in class?
% x_pred_t = predict(x_meas(:, 1), t_k);
% [diverged, lambda_bar, lambda_bar_measurement] = model_divergence(t_k, x_meas, x_pred, x_pred_t, A, delta_x_tol, alpha, x_err);

%%
nx = size(x_meas, 1);

% Numerically estimate state derivative
xdot = numer_deriv(x_meas);

% Initial sparse terms
sparse_terms = Xi0 ~= 0;

% Calculate initial residual (for testing purposes)
residual = -1 * ones([size(x_meas, 2), size(x_meas, 1), 1 + max_iter_abrupt * 3]);
residual(:, :, 1) = Theta * Xi0 - xdot;

%% Adaptive Model Fitting
if diverged
    for i = 1 : max_iter_abrupt
        %% Identify varying parameters 
        % by regressing new data onto existing sparse structure
        % note that regular regression is used
        for ind = 1 : nx
            Xi(sparse_terms, ind) = Theta(:, sparse_terms) \ xdot(:, ind);
        end

        % Calculate residual
        residual(:, :, 2 + (i - 1) * 3) = Theta * Xi - xdot;
        if sum(vecnorm(residual(:, :, 2 + (i - 1) * 3), 2, 2)) < residual_tol
            break;
        end

        %% Identify deletions 
        % by sparsly regressing on the sparse columns of
        % Theta that correspond to nonzero rows in Xi0
        for ind = 1 : nx
            [Xi(sparse_terms, ind), residual(:, :, 3 + (i - 1) * 3)] = SINDy(x_meas, xdot, Xi(sparse_terms, ind), Theta(:, sparse_terms), lambda, max_iter_sparse);
        end

        if sum(vecnorm(residual(:, :, 3 + (i - 1) * 3), 2, 2)) < residual_tol
            break;
        end

        %% Identify additions 
        % by sparsly regressing on
        % the inactive columns of Theta that correspond to zero rows in Xi0
        for ind = 1 : nx
            [Xi(~sparse_terms, ind), residual(:, :, 4 + (i - 1) * 3)] = SINDy(x_meas, residual(:, :, 3 + (i - 1) * 3), Xi(~sparse_terms, ind), Theta(:, ~sparse_terms), lambda, max_iter_sparse);
        end

        if sum(vecnorm(residual(:, :, 4 + (i - 1) * 3), 2, 2)) < residual_tol
            break;
        end
    end
end

end

