%% Constants
R_E = 6378.1363; % [km]
mu_E = 398600.4415; % [km3 / s2]

%% Q2
alt = 560; % [km]

a = alt + R_E;
r = a;
e = 0;
p = a * (1 - e ^ 2); % or just a

energy = -mu_E / (2 * a);
v = sqrt(2 * (energy + mu_E / r)); % [km / s]

% Check with formula for circular velocity
v_c = sqrt(mu_E / r);

%% Q3
a = 12 * R_E; % [km]
b = 10 * R_E; % [km]

e = sqrt(1 - (b / a) ^ 2);
p = a * (1 - e ^ 2);
r_p = p / (1 + e);
r_a = p / (1 - e);
theta_bend = 180 - atan2d(b, (r_a - a));
r_bend = sqrt(b ^ 2 + (r_a - a) ^ 2);
theta = acosd((p / r_bend - 1) / e);
r_ = p / (1 + e * cosd(theta)) / R_E;

%% Q4
e = 0.5;
r_1 = 5 * R_E; % [km]
v_1 = 5; % [km / s]
v_2 = 4; % [km / s]

energy = 1/2 * v_1 ^ 2 - mu_E / r_1;
r_2 = mu_E / (1/2 * v_2 ^ 2 - energy) / R_E;

%% Q5
e = 1.2;
theta_star = deg2rad(-6);
energy = +6.64; % [km2 / s2]

a = -mu_E / (2 * energy);
p = a * (1 - e ^ 2);
r = p / (1 + e * cos(theta_star)) / R_E;
hit_earth = r < R_E; % Hits earth
r_p = a * (1 - e) / R_E;
r_a = a * (1 + e) / R_E;