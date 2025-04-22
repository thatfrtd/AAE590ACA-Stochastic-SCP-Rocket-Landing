function [P_k, Pu_k] = recover_est_covariances(X_k, S_k)
%RECOVER_EST_COVARIANCES Summary of this function goes here
%   Detailed explanation goes here
nx = size(X_k, 1);
nu = size(S_k, 1);
nk = (-1 + sqrt(1 + 8 * size(X_k, 2) / nx)) / 2;

tri = @(k) k * (k + 1) / 2 * nx;

P_k = zeros([nx, nx, nk]);
Pu_k = zeros([nu, nu, nk]);

for k = 1:nk
    P_k(:, :, k) = X_k(:, (tri(k - 1) + 1):tri(k)) * X_k(:, (tri(k - 1) + 1):tri(k))';
    Pu_k(:, :, k) = S_k(:, (tri(k - 1) + 1):tri(k)) * S_k(:, (tri(k - 1) + 1):tri(k))';
end
end

