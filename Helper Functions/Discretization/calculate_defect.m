function [Delta] = calculate_defect(prob, i)
%CALCULATE_DEFECT Summary of this function goes here
%   If reference state, control, and parameter vectors satisfy the
%   nonlinear dynamics, then the defects will all be zero. This means the
%   defects can probide a necessary and sufficient measure of dynamic
%   feasibility

    x_ref = prob.x_ref(:, :, i);
    x_prop = zeros([prob.nx, prob.N - 1]);

    % Propagate trajectory with actual control and dynamics
    for k = 1:(prob.N - 1)
        [~, x_prop(:, k)] = ode45(@(t, x) prob.f(t, x), [prob.t(k), prob.t(k + 1)], x_ref(:, k, i), prob.tolerances);
    end
   
    Delta = vecnorm(x_prop - x_ref(:, 2:end, i), 2, 2);
end

