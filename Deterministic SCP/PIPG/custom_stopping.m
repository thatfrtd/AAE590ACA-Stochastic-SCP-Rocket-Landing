function [terminate, zhat_inf_dj, what_inf_dj] = custom_stopping(xhat_jp1, xihat_jp1, uhat_jp1, shat_jp1, what_jp1, xhat_j, xihat_j, uhat_j, shat_j, what_j, tol_abs, tol_rel)
%CUSTOM_STOPPING Summary of this function goes here
%   Detailed explanation goes here

zhat_inf_jp1 = shat_jp1;
zhat_inf_j = shat_j;
zhat_inf_dj = abs(shat_jp1 - shat_j);

zhat_inf_jp1 = max(zhat_inf_jp1, max([vecnorm(xhat_jp1, Inf, 1), vecnorm(xihat_jp1, Inf, 1), vecnorm(uhat_jp1, Inf, 1)], [], "all"));
zhat_inf_j = max(zhat_inf_j, max([vecnorm(xhat_j, Inf, 1), vecnorm(xihat_j, Inf, 1), vecnorm(uhat_j, Inf, 1)], [], "all"));
zhat_inf_dj = max(zhat_inf_dj, max([vecnorm(xhat_jp1 - xhat_j, Inf, 1), vecnorm(xihat_jp1 - xihat_j, Inf, 1), vecnorm(uhat_jp1 - uhat_j, Inf, 1)], [], "all"));

what_inf_jp1 = 0;
what_inf_j = 0;
what_inf_dj = 0;

what_inf_jp1 = max(what_inf_jp1, max(vecnorm(what_jp1, Inf, 1), [], "all"));
what_inf_j = max(what_inf_j, max(vecnorm(what_j, Inf, 1), [], "all"));
what_inf_dj = max(what_inf_dj, max(vecnorm(what_jp1 - what_j, Inf, 1), [], "all"));

if zhat_inf_dj <= tol_abs + tol_rel * max(zhat_inf_jp1, zhat_inf_j) && what_inf_dj <= tol_abs + tol_rel * max(what_inf_jp1, what_inf_j)
    terminate = true;
else
    terminate = false;
end

end

