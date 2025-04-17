function P_yk = find_P_yk(P, C_k, D_k)
% This function can find and return P_yk for 
%   the convexified covariance propogation conditions in 
%   a cvx solver
    P_yk = zeros(size(P));

    for i = 1:size(P, 3) % P is a square matrix, with a third dimension for each discretized element
        Ck = C_k(:, :, i);
        Dk = D_k(:, :, i);
        Pk = P(:, :, i);

        P_yk(:, :, i) = Ck * Pk * Ck' + Dk * Dk';
    end
end