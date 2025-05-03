function [] = plot_6DoF_MC_time_histories(t_mean, x_mean, u_mean, t_fb, x_MC_fb, u_MC_fb, t_k, X_k, S_k, t_no_fb, x_MC_no_fb, T_max, T_min, max_gimbal, use_legend)
%PLOT_3DOF_MC_TIME_HISTORIES Summary of this function goes here
%   Detailed explanation goes here

tri = @(k) k * (k + 1) / 2 * 13;

m = size(x_MC_fb, 3);

[P_k, Pu_k] = recover_est_covariances(X_k, S_k);

%% States
titles = ["X Position", "Y Position", "Z Position", "X Velocity", "Y Velocity", "Z Velocity", "Theta 1", "Theta 2", "Theta 3", "X Angular Velocity", "Y Angular Velocity", "Z Angular Velocity", "Mass"];
ylabels = ["r_x [km]", "r_y [km]", "r_z [km]", "v_x [m / s]", "v_y [m / s]", "v_z [m / s]", "\theta_1 [deg]", "\theta_2 [deg]", "\theta_3 [deg]", "\omega_x [deg / s]", "\omega_y [deg / s]", "\omega_z [deg / s]", "m [kg]"];

ops = {@(x) x, @(x) x, @(x) x, @(x) 1000 * x, @(x) 1000 * x, @(x) 1000 * x, @(x) rad2deg(x), @(x) rad2deg(x), @(x) rad2deg(x), @(x) rad2deg(x), @(x) rad2deg(x), @(x) rad2deg(x), @(x) exp(x)};

for x = 1:13
    %x_3sigbound = 3 * sqrt(squeeze(project_ellipsoid(P_k, x)))';
    x_3sigbound = 3 * sqrt(squeeze(P_k(x, x, :)))';
    x_3sigbound_cont = interp1(t_k, x_3sigbound, t_mean);
    
    figure
    tiledlayout(1, 2)
    
    nexttile
    for i = 1:m
        plot(t_fb, ops{x}(x_MC_fb(x, :, i)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
    end
    plot(t_mean, ops{x}(x_mean(x, :)), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
    plot(t_mean, ops{x}(x_mean(x, :) + x_3sigbound_cont), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
    plot(t_mean, ops{x}(x_mean(x, :) - x_3sigbound_cont), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, HandleVisibility='off'); hold off
    title(titles(x) + " vs Time with Optimized Feedback Policies")
    xlabel("Time [s]")
    ylabel(ylabels(x))
    if use_legend
        legend("Location","best")
    end
    grid on
    
    nexttile
    for i = 1:m
        plot(t_no_fb, ops{x}(x_MC_no_fb(x, :, i)), Color = [192, 192, 192] / 256); hold on
    end
    plot(t_mean, ops{x}(x_mean(x, :)), Color = "k",LineWidth=1); hold off
    title(titles(x) + " vs Time without Trajectory Corrections")
    xlabel("Time [s]")
    ylabel(ylabels(x))
    grid on
end

%% Thrust Magnitude
figure
%Pu_k = pagemtimes(S_k, pagetranspose(S_k));
%thrust_3sigbound = 3 * sqrt(squeeze(project_ellipsoid(Pu_k, 3)))';
nu = 2;
epsilon = 1e-3; % 99.9% 
S_k_norm_T = zeros([size(u_MC_fb, 2), 1]);
for j = 1:size(u_MC_fb, 2)
    S_k_norm_T(j) = norm(S_k(1:3, (tri(j - 1) + 1):tri(j)));
end
thrust_3sigbound = squeeze(sigma_mag_confidence(epsilon / 2, nu) * S_k_norm_T);

for i = 1:m
    stairs(t_fb(1:size(u_MC_fb, 2)), vecnorm(u_MC_fb(1:3, :, i), 2, 1) .* exp(x_MC_fb(13, 1:size(u_MC_fb, 2))), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
end

u_mean_full = interp1(t_k(1:size(u_mean, 2)), u_mean', t_k, "previous", "extrap")';
thrust_3sigbound_full = interp1(t_k(1:size(u_mean, 2)), thrust_3sigbound, t_k, "previous", "extrap");

stairs(t_k, vecnorm(u_mean_full(1:3, :),2,1) .* exp(x_mean(13, 1:size(u_mean_full, 2))), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
stairs(t_k, (vecnorm(u_mean_full(1:3, :),2,1) + thrust_3sigbound_full) .* exp(x_mean(13, 1:size(u_mean_full, 2))), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
stairs(t_k, (vecnorm(u_mean_full(1:3, :),2,1) - thrust_3sigbound_full) .* exp(x_mean(13, 1:size(u_mean_full, 2))), Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1,HandleVisibility='off'); hold on
yline(T_max, LineWidth = 1, LineStyle="--", Color="k", DisplayName = "Constraint"); hold on
yline(T_min, LineWidth = 1, LineStyle="--", Color="k", HandleVisibility='off'); hold off
title("Thrust Magnitude vs Time with Optimized Feedback Policies")
legend(Location="southeast")
xlabel("Time [s]")
ylabel("Thrust [kN]")
grid on

%% Gimbal Angle
figure

epsilon = 1e-3; % 99.9% 
% Best guess at how to evaluate 99.9% bound... how to turn into actual
% angle bound??
u1_3sigbound = norminv(1 - epsilon / 2) * S_k_norm_T;
u2_3sigbound = sigma_mag_confidence(epsilon / 2, 1) * S_k_norm_T;

for i = 1:m
    stairs(t_fb(1:size(u_MC_fb, 2)), acosd(u_MC_fb(1, :, i) ./ vecnorm(u_MC_fb(1:3, :, i), 2, 1)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
end

u_mean_full = interp1(t_k(1:size(u_mean, 2)), u_mean', t_k, "previous", "extrap")';
u1_3sigbound_full = interp1(t_k(1:size(u_mean, 2)), u1_3sigbound, t_k, "previous", "extrap");
u2_3sigbound_full = interp1(t_k(1:size(u_mean, 2)), u2_3sigbound, t_k, "previous", "extrap");
gimbal_99p9bound_up = atan2d(u_mean_full(2,:) + u2_3sigbound_full .* sign(u_mean_full(2,:)), u_mean_full(1, :) - u1_3sigbound_full);
gimbal_99p9bound_dn = atan2d(u_mean_full(2,:) - u2_3sigbound_full .* sign(u_mean_full(2,:)), u_mean_full(1, :) + u1_3sigbound_full);

stairs(t_k, acosd(u_mean_full(1, :) ./ vecnorm(u_mean_full(1, :), 2, 1)), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
stairs(t_k, gimbal_99p9bound_up, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
stairs(t_k, gimbal_99p9bound_dn, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, HandleVisibility='off'); hold on
yline(rad2deg(max_gimbal), LineWidth = 1, LineStyle="--", Color="k", DisplayName = "Constraint"); hold on
yline(-rad2deg(max_gimbal), LineWidth = 1, LineStyle="--", Color="k", HandleVisibility='off'); hold off
title("Gimbal Angle vs Time with Optimized Feedback Policies")
legend(Location="southeast")
xlabel("Time [s]")
ylabel("Gimbal Angle [deg]")
grid on

%% Roll Angle
figure

epsilon = 1e-3; % 99.9% 
S_k_norm_gamma = zeros([size(u_MC_fb, 2), 1]);
for j = 1:size(u_MC_fb, 2)
    S_k_norm_gamma(j) = norm(S_k(4, (tri(j - 1) + 1):tri(j)));
end
gamma_3sigbound = squeeze(norminv(1 - epsilon / 2) * S_k_norm_gamma);

% Best guess at how to evaluate 99.9% bound... how to turn into actual
% angle bound??

for i = 1:m
    stairs(t_fb(1:size(u_MC_fb, 2)), rad2deg(u_MC_fb(4, :, i)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
end

stairs(t_k, u_mean_full(4, :), Color = "k",LineWidth=1, DisplayName="Nominal"); hold on
stairs(t_k, gamma_3sigbound, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, DisplayName="99.9% Bound"); hold on
stairs(t_k, -gamma_3sigbound, Color = [100, 100, 100] / 256, LineStyle=":", LineWidth=1, HandleVisibility='off'); hold off
title("Roll Angle vs Time with Optimized Feedback Policies")
legend(Location="southeast")
xlabel("Time [s]")
ylabel("Roll Angle [deg]")
grid on

 
end

