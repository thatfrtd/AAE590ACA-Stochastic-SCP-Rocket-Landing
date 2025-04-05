function [prob] = convexify_along_reference(prob, x_ref, u_ref, p_ref)
%CONVEXIFY_ALONG_REFERENCE Summary of this function goes here
%   Detailed explanation goes here
% Normalize - should already be right?

% Linearize - not needed, done symbolically

% Discretize
prob = prob.discretize(x_ref, u_ref, p_ref);

% Constraint Approximation - chosen per constraint

% BC Approximation

end