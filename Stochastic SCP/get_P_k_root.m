function P_k_root = get_P_k_root(P, C_k, D_k)
% get_P_k_root finds the P_{yk^-}^{1/2}
%   This is needed in the convexified covariance propogation 
%   in the stochastic subproblem solution

    P_k_root = zeros(size(P));

    P_yk = find_P_yk(P, C_k, D_k);
    for i = 1:numel(P, 4)
        P_k_root(:, :, i) = chol(P_yk(:, :, i), 'lower');
    end
end