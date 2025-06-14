function [qhat, Ahat_1, Ahat_2, Bhat_minus, Bhat_plus, Shat, l, l_inv] = custom_L2_hypersphere_preconditioning(w_vse, w_tr, w_tr_s, q, A, B_minus, B_plus, S, lambda)
%L2_HYPERSPHERE_PRECONDITIONING Conditions bivariate quadratic function to
%have condition number 1
%   Makes bivariate quadratic function, which is in the objective function,
%   have level sets which are L2-norm hyperspheres and have condition
%   number of 1. Should make cut the number of solver iterations by 5 to 10
%   times. Specifically customized for optimal control
% w_vse, w_tr, w_tr_s > 0

l.x1 = sqrt(w_tr + w_vse);
l.x2 = -w_vse / l.x1;
l.xi = sqrt(w_tr * w_vse) / l.x1;
l.u = sqrt(w_tr);
l.s = sqrt(w_tr_s);

l_inv.x1 = 1 / l.x1;
l_inv.x2 = -l.x2 / l.x1 / l.xi;
l_inv.xi = 1 / l.xi;
l_inv.u = 1 / l.u;
l_inv.s = 1 / l.s;

Ahat_1 = l_inv.x1 * A;
Ahat_2 = l_inv.x2 * A;
Bhat_minus = l_inv.u * B_minus;
Bhat_plus = l_inv.u * B_plus;
Shat = l_inv.s * S;

qhat.x = lambda * l_inv.x1 * q.x;
qhat.xi = lambda * (l_inv.x2 * q.x + l_inv.xi * q.xi);
qhat.u = lambda * l_inv.u * q.u;
qhat.s = lambda * l_inv.s * q.s;

end

