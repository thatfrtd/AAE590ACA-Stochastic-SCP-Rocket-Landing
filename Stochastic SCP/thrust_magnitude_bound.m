function [Gamma_k] = thrust_magnitude_bound(S_k_ref, u_ref, k, t_k, T_max, m_0, alpha, nu)
%GAMMA_K Summary of this function goes here
%   Detailed explanation goes here
tri = @(k) k * (k + 1) / 2 * 5;

delta_t = diff(t_k);

TWR_0 = T_max / m_0;

TWR_multiplier_k = exp(einsum(@(i) alpha * delta_t(i) ...
    * max(0, norm(u_ref(1:2, i)) ...
           - sigma_mag_confidence(1e-3 / (2 * k), nu) * norm(S_k_ref(:, (tri(i - 1) + 1):tri(i)))), 1:(k - 1)));

Gamma_k = TWR_0 * TWR_multiplier_k;
end

