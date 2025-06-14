function [terminate, zhat_inf_dj, what_inf_dj] = stopping(zhat_jp1, what_jp1, zhat_j, what_j, tol_abs, tol_rel)
%STOPPING Summary of this function goes here
%   Detailed explanation goes here

zhat_inf_jp1 = norm(zhat_jp1, Inf);
zhat_inf_j = norm(zhat_j, Inf);
zhat_inf_dj = norm(zhat_jp1 - zhat_j, Inf);

what_inf_jp1 = norm(what_jp1, Inf);
what_inf_j = norm(what_j, Inf);
what_inf_dj = norm(what_jp1 - what_j, Inf);

if zhat_inf_dj <= tol_abs + tol_rel * max(zhat_inf_jp1, zhat_inf_j) && what_inf_dj <= tol_abs + tol_rel * max(what_inf_jp1, what_inf_j)
    terminate = true;
else
    terminate = false;
end

end

