function [K] = recover_gain_matrix(X, S)
%RECOVER_GAIN_MATRIX Summary of this function goes here
%   Detailed explanation goes here
nx = size(X, 1);
nu = size(S, 1);
nk = (-1 + sqrt(1 + 8 * size(X, 2) / nx)) / 2;

tri = @(k) k * (k + 1) / 2 * nx;

K = zeros([nu, nx, nk]);

for k = 1:nk
    X_k = X(:, (tri(k - 1) + 1):tri(k));
    S_k = S(:, (tri(k - 1) + 1):tri(k));
    K_k = S_k * (pinv(X_k' * X_k) * X_k');
    %K_k2 = S_k * ((X_k' * X_k) \ X_k')
    K(:, :, k) = K_k;
end
end