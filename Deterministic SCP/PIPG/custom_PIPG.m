function [x_star, xi_star, u_star, s_star, what_star, sol_info] = custom_PIPG(qhat, Ahat_1, Ahat_2, Bhat_minus, Bhat_plus, Shat, d, D, l, l_inv, lambda, sigma, rho, tol_abs, tol_rel, j_check, j_max, x_ref, xi_ref, u_ref, s_ref, what_ref)
%CUSTOM_PIPG Proportional Integral Projected Gradient convex conic solver specifically for trajectory optimization
%   Extrapolated Proportional Integral Projected Gradient (xPIPG)
%   customized for optimal control. This uses virtual state instead of 
%   virtual control where state constraints are imposed on a virtual state 
%   and the distance between the virtual state and state is penalized making 
%   dynamical and control constraints satisfied at every iteration and state 
%   constraints being satisfied asymptotically through PTR.
% rho in [1.5, 1.9] is usually a good choice - should converge about twice
% as fast as not doing extrapolation (rho = 1)

sol_info = [];
check_infeasibility = true;
sol_info.primal_infeasible = false;
sol_info.dual_infeasible = false;
sol_info.solution_status = "Unsolved";

N = size(qhat.x, 2);
qhat.x = reshape(qhat.x, [], 1, N);
qhat.xi = reshape(qhat.xi, [], 1, N);
qhat.u = reshape(qhat.u, [], 1, N);

t1 = tic;

%% Transform previous primal solution
xhat_ref = reshape(l.x1 * x_ref + l.x2 * xi_ref, [], 1, N);
xihat_ref = reshape(l.xi * xi_ref, [], 1, N);
uhat_ref = reshape(l.u * u_ref, [], 1, N);
shat_ref = reshape(l.s * s_ref, [], 1, 1);

%% Initialize transformed primal variables
x_zeta_j = xhat_ref;
xi_zeta_j = xihat_ref;
u_zeta_j = uhat_ref;
s_zeta_j = shat_ref;

%% Intialize transformed dual variable
eta_j = reshape(what_ref, [], 1, N - 1);
what_j = reshape(what_ref, [], 1, N - 1);

%% Calculate step-size (sigma from power iteration)
alpha = 2 / (lambda + sqrt(lambda^2 + 4 * sigma));

%% Initialize
xhat_j = xhat_ref;
xihat_j = xihat_ref;
uhat_j = uhat_ref;
shat_j = shat_ref;

xhat_jp1 = zeros(size(xhat_j));
xihat_jp1 = zeros(size(xihat_j));
uhat_jp1 = zeros(size(uhat_j));

for j = 1 : j_max
    xhat_jp1(:, :, 1) = x_zeta_j(:, :, 1) - alpha * (lambda * x_zeta_j(:, :, 1) + qhat.x(:, :, 1) + Ahat_1(:, :, 1)' * eta_j(:, :, 1));
    xihat_jp1(:, :, 1) = l.xi * D.ic_xi.project_onto_cones(1, l_inv.xi * (xi_zeta_j(:, :, 1) - alpha * (lambda * xi_zeta_j(:, :, 1) + qhat.xi(:, :, 1) + Ahat_2(:, :, 1)' * eta_j(:, :, 1))));
    uhat_jp1(:, :, 1) = l.u * D.ic_u.project_onto_cones(1, l_inv.u * (u_zeta_j(:, :, 1) - alpha * (lambda * u_zeta_j(:, :, 1) + qhat.u(:, :, 1) + Bhat_minus(:, :, 1)' * eta_j(:, :, 1))));

    S = sum(pagemtimes(pagetranspose(Shat), eta_j), 3);

    %% Projected gradient step
    xhat_jp1(:, :, 2 : (N - 1)) = x_zeta_j(:, :, 2 : (N - 1)) - alpha * (lambda * xhat_j(:, :, 2 : (N - 1)) + qhat.x(:, :, 2 : (N - 1)) + pagemtimes(pagetranspose(Ahat_1(:, :, 2 : (N - 1))), eta_j(:, :, 2 : (N - 1))) - l_inv.x1 * eta_j(:, :, 1 : (N - 2)));
    xihat_jp1(:, :, 2 : (N - 1)) = l.xi * D.xi.project_onto_cones(2 : (N - 1), l_inv.xi * (xi_zeta_j(:, :, 2 : (N - 1)) - alpha * (lambda * xihat_j(:, :, 2 : (N - 1)) + qhat.xi(:, :, 2 : (N - 1)) + pagemtimes(pagetranspose(Ahat_2(:, :, 2 : (N - 1))), eta_j(:, :, 2 : (N - 1))) - l_inv.x2 * eta_j(:, :, 1 : (N - 2)))));
    uhat_jp1(:, :, 2 : (N - 1)) = l.u * D.u.project_onto_cones(2 : (N - 1), l_inv.u * (u_zeta_j(:, :, 2 : (N - 1)) - alpha * (lambda * uhat_j(:, :, 2 : (N - 1)) + qhat.u(:, :, 2 : (N - 1)) + pagemtimes(pagetranspose(Bhat_minus(:, :, 2 : (N - 1))), eta_j(:, :, 2 : (N - 1))) + pagemtimes(pagetranspose(Bhat_plus(:, :, 1 : (N - 2))), eta_j(:, :, 1 : (N - 2))))));

    xhat_jp1(:, :, N) = x_zeta_j(:, :, N) - alpha * (lambda * x_zeta_j(:, :, N) + qhat.x(:, :, N) - l_inv.x1 * eta_j(:, :, N - 1));
    xihat_jp1(:, :, N) = l.xi * D.tc_xi.project_onto_cones(N, l_inv.xi * (xi_zeta_j(:, :, N) - alpha * (lambda * xi_zeta_j(:, :, N) + qhat.xi(:, :, N) - l_inv.x2 * eta_j(:, :, N - 1))));
    uhat_jp1(:, :, N) = l.u * D.u.project_onto_cones(N, l_inv.u * (u_zeta_j(:, :, N) - alpha * (lambda * u_zeta_j(:, :, N) + qhat.u(:, :, N) - pagemtimes(pagetranspose(Bhat_plus(:, :, N - 1)), eta_j(:, :, N - 1)))));
    shat_jp1 = l.s * D.s.project_onto_cones(1, l_inv.s * (s_zeta_j - alpha * (lambda * s_zeta_j + qhat.s + S)));

    %% Proportional-integral feedback of affine equality constraint violation
    what_jp1 = eta_j + alpha * (pagemtimes(Ahat_1, 2 * xhat_jp1(:, :, 1 : (N - 1)) - x_zeta_j(:, :, 1 : (N - 1))) ...
                              + pagemtimes(Ahat_2, 2 * xihat_jp1(:, :, 1 : (N - 1)) - xi_zeta_j(:, :, 1 : (N - 1))) ...
                              + pagemtimes(Bhat_minus, 2 * uhat_jp1(:, :, 1 : (N - 1)) - u_zeta_j(:, :, 1 : (N - 1))) ...
                              + pagemtimes(Bhat_plus, 2 * uhat_jp1(:, :, 2 : N) - u_zeta_j(:, :, 2 : N)) ...
                              + Shat * (2 * shat_jp1 - s_zeta_j) ...
                              - l_inv.x1 * (2 * xhat_jp1(:, :, 2 : N) - x_zeta_j(:, :, 2 : N)) - l_inv.x2 * (2 * xihat_jp1(:, :, 2 : N) - xi_zeta_j(:, :, 2 : N)) + d);

    %% Extrapolate transformed primal variables
    x_zeta_jp1 = (1 - rho) * x_zeta_j + rho * xhat_jp1;
    xi_zeta_jp1 = (1 - rho) * xi_zeta_j + rho * xihat_jp1;
    u_zeta_jp1 = (1 - rho) * u_zeta_j + rho * uhat_jp1;
    s_zeta_jp1 = (1 - rho) * s_zeta_j + rho * shat_jp1;

    %% Extrapolate transformed dual variable
    eta_jp1 = (1 - rho) * eta_j + rho * what_jp1;

    %% Check stopping criterion every j_check iterations
    if mod(j, j_check) == 0
        [terminate, zhat_inf_dj, what_inf_dj] = custom_stopping(xhat_jp1, xihat_jp1, uhat_jp1, shat_jp1, what_jp1, xhat_j, xihat_j, uhat_j, shat_j, what_j, tol_abs, tol_rel);
        sol_info.zhat_inf_dj(j / j_check) = zhat_inf_dj;
        sol_info.what_inf_dj(j / j_check) = what_inf_dj;

       if terminate
           sol_info.iterations = j;
           sol_info.solution_status = "Optimal";

           break;
       elseif check_infeasibility % Try to check infeasibility
           % if vecnorm(xhat_jp1 - xhat_j, 2, 1) > tol_infeas
           %     sol_info.primal_infeasibility = true;
           % end
           % if vecnorm([xihat_jp1 - xihat_j; uhat_jp1 - uhat_j; shat_jp1 - shat_j], 2, 1) > tol_infeas 
           %     sol_info.dual_infeasibility = true;
           % end
       end
    end

    %% Update values for next iteration
    xhat_j = xhat_jp1;
    xihat_j = xihat_jp1;
    uhat_j = uhat_jp1;
    shat_j = shat_jp1;

    x_zeta_j = x_zeta_jp1;
    xi_zeta_j = xi_zeta_jp1;
    u_zeta_j = u_zeta_jp1;
    s_zeta_j = s_zeta_jp1;

    eta_j = eta_jp1;
end

%% Recover original primal variables
x_star = squeeze(l_inv.x1 * xhat_jp1() + l_inv.x2 * xihat_jp1);
xi_star = squeeze(l_inv.xi * xihat_jp1);
u_star = squeeze(l_inv.u * uhat_jp1);
s_star = squeeze(l_inv.s * shat_jp1);

%% Retain transformed dual variable
what_star = squeeze(what_jp1);

t2 = toc(t1);

fprintf("Custom PIPG Time: %.3f ms\n", t2 * 1000 )

end

