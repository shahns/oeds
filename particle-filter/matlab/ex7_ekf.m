%%Script Excercise 7-EKF: Nimish Shah s2088894
close all; clear; clc;

%% load data file
load('z_yacht2.mat');
stop = length(fi0);
load('pf_est_5000.mat');
%% initializations
% system model
u = [400 * ones(1, length(fi0)); fi0];

sigma_w_xi = 10^-2;
sigma_w_v = 10^-2;
sigma_w_a = 10^-2;
sigma_w_phi = 0.5;
sigma_w_t = 8;

Cw = diag([sigma_w_xi^2, sigma_w_xi^2, ...
    sigma_w_v^2, sigma_w_v^2, ...
    sigma_w_a^2, sigma_w_a^2, ...
    sigma_w_t^2, sigma_w_phi^2]);

% measurement model
sigma_v_n1 = 1; % compass
sigma_v_n2 = 0.3; % log-based
sigma_v_n3 = 1; % beacon

Cn = diag([sigma_v_n1^2, sigma_v_n2^2, sigma_v_n3^2]);

x0 = [5000; 10000];

% initial estimate of cov matrix and means
sigma_x_xi = 100;
sigma_x_v = 2;
sigma_x_a = 0.04;
sigma_x_t = 300;
sigma_x_phi = 10;

mu_x_xi = 0;
mu_x_v = 0;
mu_x_a = 0;
mu_x_t = 400;
mu_x_phi = 0;

Cx0 = diag([sigma_x_xi^2, sigma_x_xi^2, ...
    sigma_x_v^2, sigma_x_v^2, ...
    sigma_x_a^2, sigma_x_a^2, ...
    sigma_x_t^2, sigma_x_phi^2]);

%% EKF initializations
C = NaN(8, 8, stop);
C(:, :, 1) = Cx0;

S = NaN(3, 3, stop);
K = NaN(8, 3, stop-1);

x_pred = NaN(8, stop);
x_pred(:, 1) = [mu_x_xi, mu_x_xi, mu_x_v, ...
    mu_x_v, mu_x_a, mu_x_a, mu_x_t, mu_x_phi]';


cov_pred = NaN(8, 8, stop);
cov_pred(:, :, 1) = Cx0;

x_update = NaN(8, stop-1);
cov_update = NaN(8, 8, stop-1);

nis = NaN(1, length(imeas));

%% EKF implementation
for i = 2:stop
    % update
    isMeasAvailable = find(imeas == i-1);
    if isMeasAvailable % we do have a meas
        S(:, :, i-1) = Hjacobian(x_pred(:, i-1), x0) * cov_pred(:, :, i-1) * Hjacobian(x_pred(:, i-1), x0)' + Cn;
        K(:, :, i-1) = cov_pred(:, :, i-1) * Hjacobian(x_pred(:, i-1), x0)' * pinv(S(:, :, i-1));
        z_samples = hmeas(x_pred(:, i-1), x0);
        residual(1:2:3, :) = wrapTo180(z(1:2:3, isMeasAvailable)-z_samples(1:2:3, :));
        residual(2, :) = z(2, isMeasAvailable) - z_samples(2, :);
        x_update(:, i-1) = x_pred(:, i-1) + K(:, :, i-1) * residual;
        cov_update(:, :, i-1) = cov_pred(:, :, i-1) - K(:, :, i-1) * S(:, :, i-1) * K(:, :, i-1)';
        % nis calculation
        nis(isMeasAvailable) = residual' * pinv(S(:, :, i-1)) * residual;
    else % we do NOT have a meas
        S(:, :, i-1) = S(:, :, i-2);
        K(:, :, i-1) = K(:, :, i-2);
        x_update(:, i-1) = x_pred(:, i-1);
        cov_update(:, :, i-1) = cov_pred(:, :, i-1);
    end
    % predict
    x_pred(:, i) = fsys(x_update(:, i-1), u(:, i));
    cov_pred(:, :, i) = Fjacobian(x_update(:, i-1)) * cov_update(:, :, i-1) * Fjacobian(x_update(:, i-1))' + Cw;
end
xRange = 1:length(x_update(1, :));

%% plot for pos vs. time
figure; %xi
plot(xRange, x_update(1, :), xRange, x_update(2, :), 'LineWidth', 2);
hold on
plot(1:83401, x_est_mean(1, :), 1:83401, x_est_mean(2, :), 'LineWidth', 2);
legend(["\xi_x--EKF", "\xi_y--EKF","\xi_x--PF", "\xi_y--PF"],'Location','southwest');
xlabel('time');
ylabel('\xi in m');
axis([0 max(xRange) -Inf inf])
%% plot for vel vs. time
figure; %v
plot(xRange, x_update(3, :), xRange, x_update(4, :), 'LineWidth', 2);
hold on
plot(1:83401, x_est_mean(3, :), 1:83401, x_est_mean(3, :), 'LineWidth', 2);
legend(["v_x--EKF", "v_y--EKF","v_x--PF", "v_y--PF"],'Location','southwest');
xlabel('time');
ylabel('v in ms^{-1}');
axis([0 max(xRange) -Inf inf])
%% plot for acc vs. time
figure; %a
plot(xRange, x_update(5, :), xRange, x_update(6, :), 'LineWidth', 2);
hold on
plot(1:83401, x_est_mean(5, :), 1:83401, x_est_mean(6, :), 'LineWidth', 2);
legend(["a_x--EKF", "a_y--EKF","a_x--PF", "a_y--PF"],'Location','northeast');
xlabel('time');
ylabel('a in ms^{-2}');
axis([0 max(xRange) -Inf inf])
%% plot for t vs. time
figure; %t
plot(xRange, x_update(7, :), 'LineWidth', 2);
hold on
% plot(1:83401, x_est_mean(7, :),'LineWidth', 2);
legend(["EKF", "PF"],'Location','southeast');
xlabel('time');
ylabel('Thrust in N');
axis([0 max(xRange) -Inf inf])
%% plot for phi vs. time
figure; %phi
plot(xRange, x_update(8, :), 'LineWidth', 2);
hold on
plot(1:83401, x_est_mean(8, :), 'LineWidth', 2);
legend(["EKF", "PF"],'Location','southeast');
xlabel('time');
ylabel('\phi in degrees');
axis([0 max(xRange) -Inf inf])
%% plot for locations with uncertainity
figure;
set(gcf,'renderer','Painters');
plot(x_update(1, :), x_update(2, :), 'LineWidth', 2);
hold on
for i = 1:500:max(xRange)
    plot_cov_ellipse(x_update(1:2, i), cov_update(1:2, 1:2, i));
end
legend(["positions", "uncertainty ellipse", "\mu"],'Location','southwest')
xlabel('\xi_x')
ylabel('\xi_y')

%% NIS consistency chk
nis95 = prctile(nis, 95);
figure;
plot(1:length(nis), nis, 0:0.5:length(nis)+2, nis95*ones(1, 283), '--', 'LineWidth', 2)
axis([0, 139, 0, 25])
legend('NIS', 'NIS_{95}')
xlabel('Measurement locations $$i_m$$', 'Interpreter', 'Latex', 'FontSize', 14);
ylabel('$$\hat{\mathbf z}(i_m)$$', 'Interpreter', 'Latex', 'FontSize', 14);