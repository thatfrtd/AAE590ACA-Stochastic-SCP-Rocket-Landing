function [sigma] = power_iteration(Hhat, z, options)
%POWER_ITERATION Computes max eigenvalue squared of matrix
%   For general matrix
arguments
    Hhat
    z % should be randomly generated
    options.tol_abs = 1e-3
    options.tol_rel = 1e-3
    options.eps_buff = 0.1
    options.j_max = 5000
end

sigma = norm(z);

for j = 1 : options.j_max
    w = 1 / sigma * Hhat * z;
    z = Hhat.' * w;
    sigma_star = norm(z);
    
    if abs(sigma_star - sigma) <= options.tol_abs + options.tol_rel * max(sigma_star, sigma)
        break;
    elseif j < options.j_max
        sigma = sigma_star;
    end
end

sigma = (1 + options.eps_buff) * sigma_star;

end