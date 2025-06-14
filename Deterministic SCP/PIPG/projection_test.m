%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSP ASA
% AC Optimal Control
% Author: Travis Hastreiter 
% Created On: 14 May, 2025
% Description: Test of analytical projections on convex sets. Also test if
% non-branching versions are equivalent
% Most Recent Change: 14 May, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Test Type 1 SOC Projection
n = 3;
m = 1000;

alpha = ones([1, m]) * tand(20);
soc = SOC(1:n, alpha);
soc.branching = true;
soc_nbrnch = soc;
soc_nbrnch.branching = false;

z = randn(n, m) * 10;

[~, sort_i] = sort(z(1, :));
z = z(:, sort_i);

cutoff = round(m / 1);
proj_z = [soc.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];
proj_z_nbrnch = [soc_nbrnch.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];

assert(sum(proj_z ~= proj_z_nbrnch, "all") == 0, "SOC branch ~= no branch")

scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "+"); hold on
scatter3(proj_z_nbrnch(1, :), proj_z_nbrnch(2, :), proj_z_nbrnch(3, :), "x"); hold off
axis equal

soc_cval = soc.eval_constraint_functions(proj_z);
assert(sum(soc_cval{1} > 1e-12, "all") == 0, "SOC constraint not meet")

%% Test Type 2 SOC Projection
n = 3;
m = 1000;

theta = ones([1, m]) * deg2rad(60);
u = ones([n, m]) .* [sqrt(2) / 2; sqrt(2) / 2; 0];
soc2 = SOC2(1:n, theta, u);
soc2.branching = true;
soc2_nbrnch = soc2;
soc2_nbrnch.branching = false;

z = randn(n, m) * 10;

[~, sort_i] = sort(z(1, :));
z = z(:, sort_i);

cutoff = round(m / 1);
proj_z = [soc2.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];
proj_z_nbrnch = [soc2_nbrnch.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];

assert(sum(vecnorm(proj_z - proj_z_nbrnch, 2, 1), "all") < 1e-10, "SOC2 branch ~= no branch")

scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "+"); hold on
scatter3(proj_z_nbrnch(1, :), proj_z_nbrnch(2, :), proj_z_nbrnch(3, :), "x"); hold off
axis equal

soc2_cval = soc2.eval_constraint_functions(proj_z);
assert(sum(soc2_cval{1} > 1e-12, "all") == 0, "SOC2 constraint not meet")

%% Test Halfspace Projection
n = 3;
m = 1000;

zeta = ones([1, m]) * 10;
u = ones([n, m]) .* [-sqrt(2) / 2; sqrt(2) / 2; 0];
halfspace = Halfspace(1:n, u, zeta);
halfspace.branching = true;
halfspace_nbrnch = halfspace;
halfspace_nbrnch.branching = false;

z = randn(n, m) * 10;

[~, sort_i] = sort(z(1, :));
z = z(:, sort_i);

cutoff = round(m / 1);
proj_z = [halfspace.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];
proj_z_nbrnch = [halfspace_nbrnch.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];

assert(sum(proj_z ~= proj_z_nbrnch, "all") == 0, "Halfspace branch ~= no branch")

figure
scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "+"); hold on
scatter3(proj_z_nbrnch(1, :), proj_z_nbrnch(2, :), proj_z_nbrnch(3, :), "x"); hold off
axis equal

halfspace_cval = halfspace.eval_constraint_functions(proj_z);
assert(sum(halfspace_cval{1} > 1e-12, "all") == 0, "Halfspace constraint not meet")

figure
scatter(z(1, :), z(2, :)); hold on
scatter(proj_z(1, :), proj_z(2, :), "+"); hold on
scatter(proj_z_nbrnch(1, :), proj_z_nbrnch(2, :), "x"); hold on
plot(-20:0.1:20, zeta(1) + 1 / (sqrt(2) / 2) * (sqrt(2) / 2 * (-20:0.1:20))); hold off
axis equal


%% Test Halfspace2 Projection
n = 3;
m = 1000;

zeta_1 = ones([1, m]) * 1;
u_1 = ones([n, m]) .* [sqrt(2) / 2; sqrt(2) / 2; 0];
zeta_2 = ones([1, m]) * 1;
u_2 = ones([n, m]) .* [-sqrt(2) / 2; sqrt(2) / 2; 0];
halfspace2 = Halfspace2(1:n, u_1, u_2, zeta_1, zeta_2);
halfspace2.branching = true;
halfspace2_nbrnch = halfspace2;
halfspace2_nbrnch.branching = false;

z = randn(n, m) * 10;

[~, sort_i] = sort(z(1, :));
z = z(:, sort_i);

cutoff = round(m / 1);
proj_z = [halfspace2.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];
proj_z_nbrnch = [halfspace2_nbrnch.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];

assert(sum(vecnorm(proj_z - proj_z_nbrnch, 2, 1), "all") < 1e-10, "Halfspace2 branch ~= no branch")

figure
scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "+"); hold on
scatter3(proj_z_nbrnch(1, :), proj_z_nbrnch(2, :), proj_z_nbrnch(3, :), "x"); hold off
axis equal

%%
figure
scatter(z(1, :), z(2, :)); hold on
scatter(proj_z(1, :), proj_z(2, :), "+"); hold on
scatter(proj_z_nbrnch(1, :), proj_z_nbrnch(2, :), "x"); hold on
plot(-20:0.1:20, zeta_1(1) - 1 / (sqrt(2) / 2) * (sqrt(2) / 2 * (-20:0.1:20))); hold on
plot(-20:0.1:20, zeta_2(1) + 1 / (sqrt(2) / 2) * (sqrt(2) / 2 * (-20:0.1:20))); hold off
axis equal

%%

halfspace2_cval = halfspace2.eval_constraint_functions(proj_z);
%assert(sum(halfspace2_cval{1} > 1e-12, "all") == 0, "Halfspace2 constraint not meet")

%% Test Ball Projection
n = 3;
m = 1000;

r = ones([1, m]) * 9;
ball = Ball(1:n, r);

z = randn(n, m) * 10;

[~, sort_i] = sort(z(1, :));
z = z(:, sort_i);

cutoff = round(m / 1);
proj_z = [ball.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];

figure
scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "+"); hold off
axis equal

ball_cval = ball.eval_constraint_functions(proj_z);
assert(sum(ball_cval{1} > 1e-12, "all") == 0, "Ball constraint not meet")

%% Test Box Projection
n = 3;
m = 1000;

l = ones([n, m]) .* [-5; -2; -1];
u = ones([n, m]) .* [9; 9; 5];
box = Box(1:n, l, u);

z = randn(n, m) * 10;

[~, sort_i] = sort(z(1, :));
z = z(:, sort_i);

cutoff = round(m / 1);
proj_z = [box.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];

figure
scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "+"); hold off
axis equal

box_cval = box.eval_constraint_functions(proj_z);
assert(sum(box_cval{1} > 1e-12, "all") == 0, "Box constraint not meet")

%% Test Type 1 SOC-Ball Projection
n = 3;
m = 1000;

alpha = ones([1, m]) * tand(60);
r = ones([1, m]) * 12;
soc_ball = SOC_Ball(1:n, r, alpha);
soc_ball.soc.branching = true;
soc_ball_nbrnch = soc_ball;
soc_ball_nbrnch.soc.branching = false;

z = randn(n, m) * 10;

[~, sort_i] = sort(z(1, :));
z = z(:, sort_i);

cutoff = round(m / 1);
proj_z = [soc_ball.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];
proj_z_nbrnch = [soc_ball_nbrnch.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];

assert(sum(proj_z ~= proj_z_nbrnch, "all") == 0, "SOC-Ball branch ~= no branch")

scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "+"); hold on
scatter3(proj_z_nbrnch(1, :), proj_z_nbrnch(2, :), proj_z_nbrnch(3, :), "x"); hold off
axis equal

soc_ball_cval = soc_ball.eval_constraint_functions(proj_z);
assert(sum(soc_ball_cval{1} > 1e-12, "all") == 0, "SOC-Ball constraint not meet")

%% Test Type 2 SOC-Ball Projection
n = 3;
m = 1000;

theta = ones([1, m]) * deg2rad(60);
u = ones([n, m]) .* [sqrt(2) / 2; sqrt(2) / 2; 0];
r = ones([1, m]) * 12;
soc2_ball = SOC2_Ball(1:n, r, theta, u);
soc2_ball.soc2.branching = true;
soc2_ball_nbrnch = soc2_ball;
soc2_ball_nbrnch.soc2.branching = false;

z = randn(n, m) * 10;

[~, sort_i] = sort(z(1, :));
z = z(:, sort_i);

cutoff = round(m / 1);
proj_z = [soc2_ball.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];
proj_z_nbrnch = [soc2_ball_nbrnch.project(1 : cutoff, z(:, 1 : cutoff)), z(:, (cutoff + 1) : m)];

assert(sum(vecnorm(proj_z - proj_z_nbrnch, 2, 1), "all") < 1e-10, "SOC2-Ball branch ~= no branch")

scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "+"); hold on
scatter3(proj_z_nbrnch(1, :), proj_z_nbrnch(2, :), proj_z_nbrnch(3, :), "x"); hold off
axis equal

soc2_ball_cval = soc2_ball.eval_constraint_functions(proj_z);
assert(sum(soc2_ball_cval{1} > 1e-12, "all") == 0, "SOC2-Ball constraint not meet")


%%

D_lcvx = SOC_Ball(1:3, ones([1, m]) * 3 * sqrt(2), ones([1, m]) * 1); % LCVX
D_minmaxT = Box(3, ones([1, m]) * 1, ones([1, m]) * 3); % Min and max thrust
D_gimbal = SOC([2, 1], ones([1, m]) * tan(gimbal_max)); % Gimbal

proj_z_1 = D_lcvx.project(1:m, z);
proj_z_2 = D_minmaxT.project(1:m, proj_z_1);
proj_z = D_gimbal.project(1:m, proj_z_2);

D = [D_lcvx, D_minmaxT, D_gimbal];

proj_z = D.project_onto_cones(1:m, z);

cvals = D.eval_constraint_functions(proj_z);
cvals = cvals{:};

assert(sum(cvals > 1e-12, "all") == 0, "LCVX Min Max Thrust projection not tight")

scatter3(z(1, :), z(2, :), z(3, :)); hold on
scatter3(proj_z_1(1, :), proj_z_1(2, :), proj_z_1(3, :), "+"); hold on
scatter3(proj_z_2(1, :), proj_z_2(2, :), proj_z_2(3, :), "x"); hold on
scatter3(proj_z(1, :), proj_z(2, :), proj_z(3, :), "square"); hold off
axis equal