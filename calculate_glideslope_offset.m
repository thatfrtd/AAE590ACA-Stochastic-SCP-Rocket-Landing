function [h] = calculate_glideslope_offset(sigma_r, glideslope_angle)
%CALCULATE_GLIDESLOPE_OFFSET Summary of this function goes here
%   Detailed explanation goes here

h = norm(sigma_r(1:(end - 1))) / tan(glideslope_angle) - sigma_r(end);

end

