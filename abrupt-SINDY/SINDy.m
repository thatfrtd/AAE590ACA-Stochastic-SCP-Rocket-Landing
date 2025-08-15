function [Xi, residual] = SINDy(x, xdot, Xi0, Theta, lambda, max_iter)
%SINDY Sparse Identification of Nonlinear Dynamics
%   Sparse Identification of Nonlinear Dynamics possibly with control.
%   Solves the sparse regression problem of 
%       xdot = Theta(x, u) * Xi
%   Should also probably try DMDc (should be worse than SINDYc because it's just linear
arguments
    x % state trajectory
    xdot % numerical derivative of state trajectory
    Xi0 % initial model coefficient matrix
    Theta % candidate library which are evaluated at x and u
    lambda = 0.025 % sparsity promoting knob - min coefficient magnitude
    max_iter = 10 % maximum iterationsend
end

% TECHNICALLY SHOULD TRANSFORM THETA TO BE UNIT VARIANCE AND ZERO MEAN
% - would have to do a least-squares fit with the untransformed resulting
%   coeffiencients to recover updated fitted model in untransformed
%   variables

nx = size(x, 1);

if ~isempty(Xi0)
    Xi = Xi0;
else
    Xi = Theta \ xdot;
end

for k = 1 : max_iter
    small_inds = abs(Xi) < lambda; % Find small coefficients
    Xi(small_inds) = 0; % Threshold

    % Regress dynamics on remaining terms to find sparse Xi
    for i = 1 : nx
        big_inds = ~small_inds(:, i); % Get remaining terms
        Xi(big_inds, i) = Theta(:, big_inds) \ xdot(:, i); % Regress dynamics
    end
end

% Calculate residual
residual = Theta * Xi - xdot;

end

