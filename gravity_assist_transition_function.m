function [x_plus, delta_v_GA, delta_v_GA_mag, r_pfb] = gravity_assist_transition_function(x_minus, x_k, x_kp1, x_p, mu, r_p_min, nr, v_scale_ref)
%GRAVITY_ASSIST_TRANSITION_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
r_i = 1:nr;
v_i = nr + (1:nr);

v_minus = x_minus(v_i);
v_p = x_p(v_i);

%% Use optimization results to calculate GA parameters
v_k = x_k(v_i);
v_kp1 = x_kp1(v_i);

% Predicted excess velocity
v_inf_k = norm(v_k - v_p) * v_scale_ref;

% Calculate predicted turn angle
delta_pred = acos(dot((v_kp1 - v_p), (v_k - v_p)) / norm(v_k - v_p) ^ 2);
if nr == 2
    delta_sign = [0 0 1] * sign(cross([v_k - v_p; 0], [v_kp1 - v_p; 0]));
elseif nr == 3
    delta_sign = sign(cross(v_k - v_p, v_kp1 - v_p));
end

% Calculate the periapsis radius of flyby (which will be used for propagation)
r_pfb = mu / v_inf_k ^ 2 * (1 / sin(delta_pred / 2) - 1);

% Check if valid flyby
r_pfb_valid = r_pfb > r_p_min;
if ~r_pfb_valid
    warning("Periapsis radius of flyby, %.5f, is less than min of %.5f", r_pfb, r_p_min)
end

%% Apply GA to propagated state
v_inf_minus = v_minus - v_p;
v_inf = norm(v_inf_minus) * v_scale_ref;

% Actual turn angle for propagation
delta_prop = 2 * asin(1 / (1 + r_pfb * v_inf ^ 2 / mu)) * delta_sign;

delta_v_GA_mag = 2 * v_inf * sin(delta_prop / 2);

% Rotate v_inf_minus by turn angle to get v_inf_plus
if nr == 2
    R = make_R2(delta_prop);
elseif nr == 3
    R = make_R(delta_prop, 3); % Needs more frame conversions 
end
v_inf_plus = R * v_inf_minus;

v_plus = v_inf_plus + v_p;

delta_v_GA = v_plus - v_minus;

x_plus = [x_minus(r_i); v_plus; x_minus(end)];

end

