function [Delta] = calculate_defect(prob, x_ref, u_ref, p_ref)
%CALCULATE_DEFECT Summary of this function goes here
%   If reference state, control, and parameter vectors satisfy the
%   nonlinear dynamics, then the defects will all be zero. This means the
%   defects can probide a necessary and sufficient measure of dynamic
%   feasibility
    
    x_prop = zeros([prob.n.x, prob.N - 1]);
    t_k = linspace(0, prob.tf, prob.N);

    % Propagate trajectory with actual control and dynamics
    for k = 1:(prob.N - 1)
        if prob.u_hold == "ZOH"
            u_ref_k = @(t) u_ref(:, k);
        elseif prob.u_hold == "FOH"
            u_ref_k = @(t) interp1([t_k(k), t_k(k + 1)], [u_ref(:, k), u_ref(:, k + 1)], t);
        end
        [~, x_prop_k] = ode45(@(t, x) prob.cont.f(t, x, u_ref_k(t), p_ref), [t_k(k), t_k(k + 1)], x_ref(:, k), prob.tolerances);
        x_prop(:, k) = x_prop_k(end, :)';
    end
   
    Delta = vecnorm(x_prop - x_ref(:, 2:end), 2, 1);
end

