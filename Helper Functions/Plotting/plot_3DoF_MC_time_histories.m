function [] = plot_3DoF_MC_time_histories(t, x_mean, u_mean, x_MC_fb, u_MC_fb, X_k, S_k, x_MC_no_fb)
%PLOT_3DOF_MC_TIME_HISTORIES Summary of this function goes here
%   Detailed explanation goes here

m = size(x_MC_fb, 3);

figure
tiledlayout(1, 2)

nexttile
for i = 1:m
    plot(t, x_MC_fb(:, 1, i), Color = [192, 192, 192] / 256); hold on
end
plot(t, x_mean(1, :), Color = "k",LineWidth=1); hold off
title("X Position vs Time with Optimized Feedback Policies")
ylabel("Time [s]")
xlabel("X [km]")

nexttile
for i = 1:m
    plot(t, x_MC_no_fb(:, 1, i), Color = [192, 192, 192] / 256); hold on
end
plot(t, x_mean(1, :), Color = "k",LineWidth=1); hold off
title("X Position vs Time without Trajectory Corrections")
ylabel("Time [s]")
xlabel("X [km]")

end

