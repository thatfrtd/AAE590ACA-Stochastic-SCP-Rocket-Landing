%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project 

% Author: Travis Hastreiter 
% Created On: 6 April, 2025
% Description: 3DoF landing of rocket using PTR SCP algorithm
% Most Recent Change: 6 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize
% Vehicle Parameters
u_max = 0.1; % []
alpha = 1; % [?]
 
AU = 149597898; % [km]

a_venus = 108207284 / AU; % [AU]
nu0_venus = deg2rad(-50); % [rad]
R_venus = 6351.0; % [km]
mu_venus = 324858.59882646; % [km3 / s2] 

a_earth = 1; % [AU]
nu0_earth = 0; % [rad]

a_mars = 1.524; % [AU]
nu0_mars = pi; % [rad]

m_0 = 1;

% All times are nondimensionalized by Earth's mean motion
tf = 8; % [] arrival time (nondimensionalized)

mu = 1; % [] Sun's gravitaional parameter (nondimensionalized)
mu_sun = 132712440017.99; % [km3 / s2]

% Problem Parameters
N = 60; % []

nr = 2;

%% Create initial conditions
P_earth_over_P_mars = (a_earth / a_mars) ^ (3/2);
P_earth_over_P_venus = (a_earth / a_venus) ^ (3/2);

earth_pos = @(t) a_earth .* [cos(t + nu0_earth); ...
                             sin(t + nu0_earth)];
mars_pos = @(t) a_mars .* [cos(t * P_earth_over_P_mars + nu0_mars); ...
                           sin(t * P_earth_over_P_mars + nu0_mars)];
venus_pos = @(t) a_venus .* [cos(t * P_earth_over_P_venus + nu0_venus); ...
                           sin(t * P_earth_over_P_venus + nu0_venus)];

nuf_mars = @(tf) tf * P_earth_over_P_mars + nu0_mars;
nu_venus = @(t) t * P_earth_over_P_venus + nu0_venus;

x_0 = [earth_pos(0); v_circ(earth_pos(0), nu0_earth, mu); m_0];
x_f = @(tf) [mars_pos(tf); v_circ(mars_pos(tf), nuf_mars(tf), mu)];

x_venus = @(t) [venus_pos(t); v_circ(venus_pos(t), nu_venus(t), mu)];

v_scale_ref = norm(v_circ(earth_pos(0) * AU, nu0_earth, mu_sun));

%% Finish setting up problem
venus_GA_transition_k = round(N / 3);
delta_t = tf / (N - 2);
venus_GA_t = delta_t * (venus_GA_transition_k - 1);

tspan = [0, 2];
t_k = [linspace(0, 1, venus_GA_transition_k), linspace(1, 2, (N - venus_GA_transition_k))];
t_k(venus_GA_transition_k + 1) = 1 + 1e-12;

u_hold = "ZOH";
Nu = (u_hold == "ZOH") * (N - 1) + (u_hold == "FOH") * N;

nx = 5;
nu = 2;
np = 1;

initial_guess   = "straight line"; % "CasADi" or "straight line"

% PTR algorithm parameters
ptr_ops.iter_max = 10;
ptr_ops.iter_min = 5;
ptr_ops.Delta_min = 5e-5;
ptr_ops.w_vc = 1e3;
ptr_ops.w_tr = ones(1, Nu) * 1e-2;
ptr_ops.w_tr_p = ones(1, np) * 1e-1;
ptr_ops.update_w_tr = false;
ptr_ops.delta_tol = 2e-3;
ptr_ops.q = 2;
ptr_ops.alpha_x = 1;
ptr_ops.alpha_u = 1;
ptr_ops.alpha_p = 0;

scale = false;

%% Specify Constraints
% Convex state path constraints
final_time_constraint = {1, @(t, x, u, p, k) sum(p) - 8};
state_convex_constraints = {};

% Convex control constraints
max_control_constraint = {1:N, @(t, x, u, p, k)  norm(u(1:2)) - u_max};
control_convex_constraints = {max_control_constraint};

% Combine convex constraints
convex_constraints = [state_convex_constraints, control_convex_constraints];

% Nonconvex state path constraints
state_nonconvex_constraints = {};

% Nonconvex control constraints
control_nonconvex_constraints = {};

nonconvex_constraints = [state_nonconvex_constraints, control_nonconvex_constraints];

% Terminal boundary conditions
terminal_bc = @(t, x, u, p) [x(1:4) - x_f(p(1)); 0];
terminal_bc_linearized = linearize_constraint(terminal_bc, nx, nu, np, "p", 1:2);
terminal_bc_linearized = @(x, p, x_ref, p_ref) terminal_bc_linearized(0, x, 0, p, x_ref, 0, p_ref);

%% Phases
% Phase transitions
phase_transition_k = [venus_GA_transition_k];

[venus_GA_phase_transition.convex_constraints, venus_GA_phase_transition.nonconvex_constraints, original_constraints] = construct_GA_constraints_fixedtf(venus_GA_transition_k, t_k, mu_venus, R_venus, x_venus, nr, v_scale_ref);

venus_GA_phase_transition.transition_func = @(t, x, x_k, x_kp1) gravity_assist_transition_function(x, x_k, x_kp1, x_venus(t), mu_venus, R_venus, nr, v_scale_ref);

phase_transition = {venus_GA_phase_transition};

%% Get Dynamics
f = @(t, x, u, p, k) state_equation(x, u, mu, alpha) * (p(1) * (k < phase_transition_k) + p(2) * (k > phase_transition_k));
f = multidynamic_f(f, 1:(N - 1));

%% Specify Objective
%min_fuel_objective = @(x, u, p) sum(norms(u(1:2, [1:(venus_GA_transition_k - 1), (venus_GA_transition_k + 1) : Nu]), 2, 1)) * tf / (N - 2);
min_fuel_objective = @(x, u, p, x_ref, u_ref, p_ref) einsum(@(k) norm(u(1:2, k)) * p_ref(1) / (N - 1) + norm(u_ref(1:2, k)) * (p(1) - p_ref(1)) / (venus_GA_transition_k - 1), 1:(venus_GA_transition_k - 1)) ...
                                                   + einsum(@(k) norm(u(1:2, k)) * p_ref(2) / (N - 1) + norm(u_ref(1:2, k)) * (p(2) - p_ref(2)) / (Nu - (venus_GA_transition_k + 1)), (venus_GA_transition_k + 1) : Nu);

%% Create Guess
AU_guess = [interp1([0, t_k(venus_GA_transition_k)], [a_earth, a_venus]', t_k(1 : venus_GA_transition_k))'- 1e-6; interp1([t_k(venus_GA_transition_k), tf], [a_venus, a_mars]', t_k((venus_GA_transition_k + 1) : end))' + 1e-6];
nu_guess = [interp1([0, t_k(venus_GA_transition_k)], [nu0_earth, nu_venus(t_k(venus_GA_transition_k))]', t_k(1 : venus_GA_transition_k))'; interp1([t_k(venus_GA_transition_k), tf], [nu_venus(t_k(venus_GA_transition_k)), nuf_mars(tf)]', t_k((venus_GA_transition_k + 1) : end))' + 1e-4];
r_guess = [AU_guess .* cos(nu_guess), AU_guess .* sin(nu_guess)]';
v_guess = v_circ_array(r_guess, nu_guess', mu);
m_guess = m_0 * ones([1, N]);

guess.x = [r_guess; v_guess; m_guess];
guess.u = interp1(tspan, zeros(nu, 2)', t_k(1:Nu))' + 1e-3;
phase_1_s = venus_GA_t;
phase_2_s = tf - venus_GA_t;
guess.p = [phase_1_s; phase_2_s];

%% Construct Problem Object
prob = DeterministicProblem(x_0, x_f, N, u_hold, 1, f, guess, convex_constraints, min_fuel_objective, scale = scale, terminal_bc = terminal_bc_linearized, nonconvex_constraints = nonconvex_constraints, phase_transition=phase_transition, phase_transition_k=phase_transition_k, t_k = t_k);

%%
Delta = calculate_defect(prob, guess.x, guess.u, guess.p);
norm(Delta)

%% Test Discretization on Initial Guess

[prob, Delta_disc] = prob.discretize(guess.x, guess.u, guess.p);

x_disc = prob.disc_prop(guess.u, guess.p);

[t_cont, x_cont, u_cont] = prob.cont_prop(guess.u, guess.p, x_k = guess.x);

%% Solve Problem with PTR
%ptr_sol_nophase = ptr(prob, ptr_ops);

%% Solve Phased Problem with PTR
ptr_sol = ptr_phase(prob, ptr_ops);

%%
if ~ptr_sol.converged
    ptr_sol.converged_i = ptr_ops.iter_max;
end

x = ptr_sol.x(:, :, ptr_sol.converged_i + 1);
u = ptr_sol.u(:, :, ptr_sol.converged_i + 1);
p = ptr_sol.p(:, ptr_sol.converged_i + 1);

%%
[prob, Delta_disc] = prob.discretize(x, u, p);
    
x_disc = prob.disc_prop(u, p);

[t_cont, x_cont, u_cont] = prob.cont_prop(u, p, x_k = x);



%%

i = ptr_sol.converged_i + 1;
[t_cont_sol, x_cont_sol, u_cont_sol] = prob.cont_prop(ptr_sol.u(:, :, i), ptr_sol.p(:, i), x_k = x);
%%
earth_trajectory = earth_pos(t_cont_sol');
mars_trajectory = mars_pos(t_cont_sol');
venus_trajectory = venus_pos(t_cont_sol');
x_venus_GA = x_venus(t_k(venus_GA_transition_k));

x_at_GA = x_cont_sol(:, find(t_cont_sol <= venus_GA_t,1,'last'));

[~, delta_v_GA, delta_v_GA_mag, r_pfb] = gravity_assist_transition_function(x_at_GA, x(:, venus_GA_transition_k), x(:, venus_GA_transition_k + 1), x_venus_GA, mu_venus, R_venus, nr, v_scale_ref);
r_pfb / R_venus

figure
plot(x(1, :), x(2, :), DisplayName="Minimum Fuel Transfer"); hold on
quiver(x(1, 1:Nu)',x(2, 1:Nu)',u(1, :)',u(2, :)', 1, DisplayName="Thrust"); hold on
plot(earth_trajectory(1, :), earth_trajectory(2, :), DisplayName="Earth Orbit"); hold on;
plot(mars_trajectory(1, :), mars_trajectory(2, :), DisplayName="Mars Orbit"); hold on
plot(venus_trajectory(1, :), venus_trajectory(2, :), DisplayName="Venus Orbit"); hold on
scatter(x_0(1), x_0(2), DisplayName="Earth at t = 0"); hold on
scatter(x_venus_GA(1), x_venus_GA(2), DisplayName=sprintf("Venus at t = %.4g", t_k(venus_GA_transition_k))); hold on
scatter(x_at_GA(1), x_at_GA(2), DisplayName=sprintf("Spacecraft at t = %.4g", t_k(venus_GA_transition_k)), Marker="x"); hold on
scatter(x_f(1), x_f(2), DisplayName=sprintf("Mars at t = %.4g", tf)); hold on;
quiver(x_venus_GA(1), x_venus_GA(2), delta_v_GA(1), delta_v_GA(2), 10, "filled", DisplayName=sprintf("GA Delta V %.3f km/s", delta_v_GA_mag), Color="m"); hold off
title("Minimum Fuel Transfer Between Earth and Mars")
subtitle(sprintf("For final time of %.4g", tf))
xlabel("X [AU]")
ylabel("Y [AU]")
legend(Location="eastoutside")
grid on
axis equal

%% Plot optimal Control
fuel = sum(squeeze(vecnorm(u(1:2, 1:end))) .* tf / (N - 1), 2);

figure
tiledlayout(1, 1, "TileSpacing","compact")

nexttile
if u_hold == "ZOH"
    stairs(t_k(1:Nu), u'); hold on
elseif u_hold == "FOH"
    plot(t_k(1:Nu), u); hold on
end
xline(venus_GA_t)
title("Optimal Control History", Interpreter="latex")
subtitle(sprintf("Total Fuel Used of %.3f units", fuel(end)))
xlabel("Time []")
ylabel("Control Value")
xlim([0, tf])
legend("u_1", "u_2", "Venus GA", location = "south", Orientation="horizontal")
grid on

%%
N_anim = 1 / 0.01;
t_cont = linspace(0, tf, N_anim);
x_cont = interp1(t_cont_sol, x_cont_sol', t_cont)';

guess_trajectory = [interp1(t_k(1 : venus_GA_transition_k), x(:, (1 : venus_GA_transition_k))', t_cont(t_cont < venus_GA_t))', interp1(t_k((venus_GA_transition_k + 1) : end), x(:, (venus_GA_transition_k + 1) : end)', t_cont(venus_GA_t <= t_cont))'];
earth_trajectory = earth_pos(t_cont);
mars_trajectory = mars_pos(t_cont);
venus_trajectory = venus_pos(t_cont);

figure
t = t_cont;

r1 = x_cont(1,:)';
r2 = x_cont(2,:)';
r3 = x_cont(1,:)' * 0;

r = sqrt(r1.^2 + r2.^2 + r3.^2);

spacecraft = animatedline('LineWidth',2);
spacecraft_guess = animatedline('LineWidth', 1.5);
earth = animatedline('lineWidth', 1);
mars = animatedline('lineWidth', 1);
venus = animatedline('lineWidth', 1);

hold on;

spacecraft.Color = 'blue';
spacecraft.Visible = 'on';
spacecraft_guess.Color = 'cyan';
spacecraft_guess.Visible = 'on';
earth.Color = 'green';
earth.Visible = 'on';
mars.Color = 'red';
mars.Visible = 'on';
venus.Color = 'm';
venus.Visible = 'on';

title(sprintf('Spacecraft 1\nTime: %0.2f sec | Radius: %0.2f', t(1),r(1)), 'Interpreter', 'latex','Color','k');

j_0 = int16((t(1) / t(length(t)))*(length(t))) + 1; % stating time index (integer)

for j = j_0:length(t) % this foor loop generates the animated plot with position vectors with indexes
    addpoints(spacecraft_guess,guess_trajectory(1, j),guess_trajectory(2, j),0);
    addpoints(spacecraft,r1(j),r2(j),r3(j));
    addpoints(earth,earth_trajectory(1, j), earth_trajectory(2, j), 0);
    addpoints(mars,mars_trajectory(1, j), mars_trajectory(2, j), 0);
    addpoints(venus,venus_trajectory(1, j), venus_trajectory(2, j), 0);
    guess_head = scatter3(guess_trajectory(1, j), guess_trajectory(2, j), 0, "cyan", "filled");
    head = scatter3(r1(j),r2(j),r3(j),'filled','MarkerFaceColor','b','LineWidth', 0.8);
    earth_head = scatter3(earth_trajectory(1, j), earth_trajectory(2, j), 0, "green", "filled");
    mars_head = scatter3(mars_trajectory(1, j), mars_trajectory(2, j), 0, "red", "filled");
    venus_head = scatter3(venus_trajectory(1, j), venus_trajectory(2, j), 0, "m", "filled");

    drawnow
    pause(0.01);
    delete(head);
    delete(guess_head);
    delete(earth_head);
    delete(mars_head);
    delete(venus_head);

    title(sprintf('Earth-Mars Min Fuel Fixed tf Transfer with Venus GA at %.2f \nTime: %.2f | Mass: %0.2f', t_k(venus_GA_transition_k), t(j),x_cont(5, j)),'Interpreter','Latex','Color','k');
    axis equal
    grid on
end

%% Helper Functions
function [xdot] = state_equation(x, u, mu, alpha)
    rvec = x(1:2);
    vvec = x(3:4);
    m = x(5);
    r = sqrt(rvec(1) ^ 2 + rvec(2) ^ 2);
    
    f0 = [vvec; -mu / r ^ 3 * rvec];

    B = [zeros(2); eye(2)];

    u_mag = sqrt(u(1) ^ 2 + u(2) ^ 2);

    m_dot = -alpha * u_mag;

    xdot = [f0 + B * u(1:2) / m; m_dot];
end

function [vvec] = v_circ(rvec, nu, mu)
    r = norm(rvec);
    v = sqrt(mu / r);
    vvec = v .* [-sin(nu); cos(nu)];
end

function [vvec] = v_circ_array(rvec, nu, mu)
    r = vecnorm(rvec, 2, 1);
    v = sqrt(mu ./ r);
    vvec = v .* [-sin(nu); cos(nu)];
end