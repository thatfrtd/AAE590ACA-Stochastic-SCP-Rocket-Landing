function [qhat, Hhat, L, L_inv] = L2_hypersphere_preconditioning(Q, q, H, lambda)
%L2_HYPERSPHERE_PRECONDITIONING Conditions bivariate quadratic function to
%have condition number 1
%   Makes bivariate quadratic function, which is in the objective function,
%   have level sets which are L2-norm hyperspheres and have condition
%   number of 1. Should make cut the number of solver iterations by 5 to 10
%   times.
% Q must be symmetric positive definite

L = chol(Q);
L_inv = inv(L);
Hhat = H * L_inv;
qhat = lambda * L_inv' * q;

end

