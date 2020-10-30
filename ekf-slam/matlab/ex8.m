%%Script Excercise 8 Part 1: Nimish Shah s2088894
close all; clear; clc;

%% load data file
load('SLAM.mat');

%% initializations
v = u(1, :);
phi = u(2, :);

no_iters = size(u, 2);
p = NaN(size(u, 1), no_iters+1);
p(:, 1) = [0; 0];

Cw_tilde = diag([std_velocity^2, std_heading^2]);

Cpredonly = NaN(size(u, 1), size(u, 1), no_iters);
Cpredonly(:, :, 1) = zeros(2);

%% Part I
F = eye(2);
for iter = 2:no_iters
    G = delta * (pi / 180) * [cosd(phi(iter-1)), -v(iter-1) * sind(phi(iter-1)); sind(phi(iter-1)), v(iter-1) * cosd(phi(iter-1))];
    p(:, iter) = F * p(:, iter-1) + delta * v(iter-1) * [cosd(phi(iter-1)); sind(phi(iter-1))];
    Cw = G * Cw_tilde * G';
    Cpredonly(:, :, iter) = F * Cpredonly(:, :, iter-1) * F' + Cw;
end

plot(p(1, :), p(2, :), '.', 'Color', '#0072BD');
hold on

% create animation showing evaluation of cov mat at each 100
for time_instant = 1:100:no_iters
    plot_cov_ellipse(p(:, time_instant), Cpredonly(:, :, time_instant), 'showMean', true);
end
xlabel("x(m)");
ylabel("y(m)");
legend("Predicted Location", 'Location', 'northwest')

%% Part II

%% initialize state vector
xest = cell(no_iters, 1); % pre-allocation of the estimated state vector
Cest = cell(no_iters, 1); % error covariance matrix
xpred = cell(no_iters+1, 1); % predicted state vector
Cpred = cell(no_iters+1, 1); % prediction covariance matrix

xpred{1} = zeros(2, 1);
Cpred{1} = zeros(2);

clear G;
G{1} = zeros(2);
unique_landmarks = [];
for time_instant = 1:no_iters
    unique_landmarks = vertcat(unique_landmarks, Z{time_instant}.id(:));
end

max_unique_landmarks = length(unique(unique_landmarks));

%% initialize bookkeeping structure

LMBOOK.state_vector_ind = zeros(max_unique_landmarks+2, 1); % index of the landmark in the state vector
LMBOOK.state_vector_dim = 2; % dim of state vector
LMBOOK.total_visible = zeros(no_iters, 1); % no of total landmatks visible at specific instant
LMBOOK.visible = zeros(max_unique_landmarks, no_iters); % visibility of landmarks over time

%% start SLAM loop

for time_instant = 1:no_iters
    current_meas = Z{time_instant};
    LMBOOK.total_visible(time_instant) = length(current_meas.id);
    current_visible_lms = current_meas.id;
    is_new_lm_available = false;
    is_known_lm_available = false;
    is_no_landmark_visible = false;
    
    %% landmark identification
    % list known landmarks from visible
    [current_known_lms, ~, current_known_lm_ind] = intersect(current_visible_lms, LMBOOK.state_vector_ind, 'stable');
    if ~isempty(current_known_lms)
        is_known_lm_available = true;
    end
    no_current_known_lms = length(current_known_lms);
    
    % list new landmarks from visible
    if ~is_known_lm_available
        current_new_lms = current_visible_lms;
    else
        current_new_lms = current_visible_lms(~ismember(current_visible_lms, current_known_lms));
    end
    no_current_new_lm = length(current_new_lms);
    
    if ~isempty(current_new_lms)
        is_new_lm_available = true;
    end
    
    % no landmarks visible
    if ~(is_new_lm_available || is_known_lm_available)
        is_no_landmark_visible = true;
    end
    
    %% agumentation
    if is_new_lm_available
        LMBOOK.state_vector_ind(LMBOOK.state_vector_dim/2+1:LMBOOK.state_vector_dim/2+no_current_new_lm) = current_new_lms;
        LMBOOK.state_vector_dim = LMBOOK.state_vector_dim + 2 * no_current_new_lm;
        % augment state vector
        xpred{time_instant} = [xpred{time_instant}; zeros(2*no_current_new_lm, 1)]; %% add sensor readings instead of zeros
        % augment pred covariance matrix
        temp_cpred = zeros(LMBOOK.state_vector_dim);
        old_size = size(Cpred{time_instant}, 1);
        temp_cpred(1:old_size, 1:old_size) = Cpred{time_instant};
        temp_cpred(old_size+1:end, old_size+1:end) = 1000^2 * eye(2*no_current_new_lm); % very high uncertainity about current prediction
        temp_cpred(1:2, old_size+1:end) = repmat(temp_cpred(1:2, 1:2), 1, no_current_new_lm);
        temp_cpred(old_size+1:end, 1:2) = repmat(temp_cpred(1:2, 1:2), no_current_new_lm, 1);
        temp_cpred(3:old_size, old_size+1:end) = repmat(temp_cpred(3:old_size, 1:2), 1, no_current_new_lm);
        temp_cpred(old_size+1:end, 3:old_size) = repmat(temp_cpred(1:2, 3:old_size), no_current_new_lm, 1);
        Cpred{time_instant} = temp_cpred;
    end
    [~, ~, current_visible_lm_state_ind] = intersect(current_visible_lms, LMBOOK.state_vector_ind, 'stable');
    
    %% update
    if ~is_no_landmark_visible
        
        %% measurement function and measurement vector
        z = current_meas.zpos;
        H = zeros(2*LMBOOK.total_visible(time_instant), LMBOOK.state_vector_dim);
        H(1:end, 1:2) = repmat(-eye(2), LMBOOK.total_visible(time_instant), 1);
        for lm = 1:LMBOOK.total_visible(time_instant)
            ind = current_visible_lm_state_ind(lm) - 1;
            H(2*lm-1:2*lm, 2*ind+1:2*ind+2) = eye(2);
        end
        Cv = eye(2*LMBOOK.total_visible(time_instant)) * stdn^2;
        zpred = H * xpred{time_instant};
        S = H * Cpred{time_instant} * H' + Cv;
        K = Cpred{time_instant} * H' / S;
        
        xest{time_instant} = xpred{time_instant} + K * (z(:) - zpred);
        Cest{time_instant} = Cpred{time_instant} - K * S * K';
    else
        xest{time_instant} = xpred{time_instant};
        Cest{time_instant} = Cpred{time_instant};
    end
    
    %% predict
    xpred{time_instant+1} = xest{time_instant};
    xpred{time_instant+1}(1:2) = F * xest{time_instant}(1:2) + delta * v(time_instant) * [cosd(phi(time_instant)), sind(phi(time_instant))]';
    G{time_instant+1} = delta * (pi / 180) * [cosd(phi(time_instant)), -v(time_instant) * sind(phi(time_instant)); sind(phi(time_instant)), v(time_instant) * cosd(phi(time_instant))];
    
    Cpred{time_instant+1} = Cest{time_instant};
    Cpred{time_instant+1}(1:2, 1:2) = F * Cest{time_instant}(1:2, 1:2) * F' + G{time_instant} * Cw_tilde * G{time_instant}';
    
end

%% plot
figure;
hold on
for time_instant = 1:no_iters
    xplotx(time_instant) = xest{time_instant}(1);
end

for time_instant = 1:no_iters
    xploty(time_instant) = xest{time_instant}(2);
end
plot(xplotx, xploty, '.', 'Color', '#0072BD');

for time_instant = 1:100:no_iters
    plot_cov_ellipse([xplotx(time_instant); xploty(time_instant)], Cest{time_instant}(1:2, 1:2))
end
xlabel("x(m)");
ylabel("y(m)");
legend("Estimated Location w/o smoothing", 'Location', 'northwest')

%% Part III
xrts{no_iters, 1} = xest{no_iters};
Crts{no_iters, 1} = Cest{no_iters};

for time_instant = no_iters - 1:-1:1
    Q = size(xest{time_instant+1}, 1);
    P = size(xest{time_instant}, 1);
    del = Q - P;
    
    F = eye(Q, P);
    F(end-del+1:end, :) = 0;
    
    C = Cest{time_instant} * F' / Cpred{time_instant+1};
    xrts{time_instant} = xest{time_instant} + C * (xrts{time_instant+1} - xpred{time_instant+1});
    Crts{time_instant} = Cest{time_instant} + C * (Crts{time_instant+1} - Cpred{time_instant+1}) * C';
end

%% plot
figure;
hold on
for time_instant = 2:no_iters
    xplotx_rts(time_instant) = xrts{time_instant}(1);
end

for time_instant = 2:no_iters
    xploty_rts(time_instant) = xrts{time_instant}(2);
end
plot(xplotx_rts, xploty_rts, '.', 'Color', '#0072BD');

for time_instant = 2:100:no_iters
    plot_cov_ellipse([xplotx_rts(time_instant); xploty_rts(time_instant)], Crts{time_instant}(1:2, 1:2))
end
xlabel("x(m)");
ylabel("y(m)");
legend("Estimated Location w/ RTS smoothing", 'Location', 'northwest')

%% Comparisions
% plot predicted + SLAM w/o smoothing
figure;
hold on
plot(p(1, :), p(2, :), '.', 'Color', '#0072BD');
plot(xplotx, xploty, '.', 'Color', '#D95319');
xlabel("x(m)")
ylabel("y(m)")
legend("Prediction Only", "SLAM w/o smoothing", 'Location', 'northwest');

% plot predicted + SLAM w/o smoothing + SLAM w/ smoothing
figure;
hold on
plot(p(1, :), p(2, :), '.', 'Color', '#0072BD');
plot(xplotx, xploty, '.', 'Color', '#D95319');
plot(xplotx_rts, xploty_rts, '.', 'Color', '#EDB120');
xlabel("x(m)")
ylabel("y(m)")
legend("Prediction Only", "SLAM w/o smoothing", "SLAM w/ RTS smoothing", 'Location', 'northwest');


%% performance plots

eigen_predonly = nan(2, no_iters);
for i = 1:no_iters
    eigen_predonly(:, i) = eig(Cpredonly(:, :, i));
    sqrt_max_eigen(i) = sqrt(max(eigen_predonly(:)));
end

eigen_slam = nan(2, no_iters);
for i = 1:no_iters
    eigen_slam(:, i) = eig(Cest{i}(1:2, 1:2));
    sqrt_max_eigen_slam(i) = sqrt(max(eigen_slam(:)));
end

eigen_rts = nan(2, no_iters);
for i = 2:no_iters
    eigen_rts(:, i) = eig(Crts{i}(1:2, 1:2));
    sqrt_max_eigen_rts(i) = sqrt(max(eigen_rts(:)));
end

figure;
plot(1:20:no_iters, sqrt_max_eigen(1:20:end), '-x')
hold on
plot(1:20:no_iters, sqrt_max_eigen_slam(1:20:end), '-x')
plot(1:20:no_iters, sqrt_max_eigen_rts(1:20:end), '-x')
xlabel("time instant")
ylabel("sqrt. of max eigen value (m)")
legend("Prediction only","SLAM w/o smoothing", "SLAM w/ RTS smoothing","Bearing only predition", 'Location', 'northwest');
%% Videos
% SLAM w/o smoothing
v = VideoWriter('without-smoothing.avi', 'Motion JPEG AVI');
open(v);
figure;

for time_instant=1:no_iters
    clf
    title('Estimated Vehicle Trajectory using KF-SLAM w/o Smoothing')
    ylim([-150,350])
    xlim([-50,450])
    hold on
    plot(xplotx(1:time_instant), xploty(1:time_instant), '.', 'Color', '#0072BD');
    plot_cov_ellipse([xplotx(time_instant); xploty(time_instant)], Cest{time_instant}(1:2, 1:2), 'showMean', true, 'color', 'b', 'labels', ["x (m)", "y (m)"])
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
% Loop closing only
v = VideoWriter('with-smoothing.avi', 'Motion JPEG AVI');
open(v);
figure;

for time_instant=2:no_iters
    clf
    title('Estimated Vehicle Trajectory using KF-SLAM w/ RTS Smoothing')
    ylim([-150,350])
    xlim([-50,450])
    hold on
    plot(xplotx_rts(1:time_instant), xploty_rts(1:time_instant), '.', 'Color', '#0072BD');
    plot_cov_ellipse([xplotx_rts(time_instant); xploty_rts(time_instant)], Crts{time_instant}(1:2, 1:2), 'showMean', true, 'color', 'b', 'labels', ["x (m)", "y (m)"])
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

% SLAM+loop closing video 
v = VideoWriter('loop-closing.avi', 'Motion JPEG AVI');
open(v);
figure;


for time_instant=1:no_iters
    clf
    title('SLAM Loopclosing')
    ylim([-150,350])
    xlim([-50,450])
    hold on
    plot(xplotx(1:time_instant), xploty(1:time_instant),'.', 'Color', '#0072BD')
    plot_cov_ellipse([xplotx(time_instant); xploty(time_instant)], Cest{time_instant}(1:2, 1:2), 'showMean', true, 'color', 'b', 'labels', ["x (m)", "y (m)"])

    legend("Estimated Location w/o Smoothing")
    frame = getframe(gcf);
    writeVideo(v,frame);
end

for time_instant=no_iters:-1:2
    clf
    hold on
    plot(xplotx, xploty,'.', 'Color', '#0072BD')
    title('SLAM Loopclosing')
    ylim([-150,350])
    xlim([-50,450])
    
    plot(xplotx_rts(no_iters:-1:time_instant), xploty_rts(no_iters:-1:time_instant), '.', 'Color', '#D95319');
    plot_cov_ellipse([xplotx_rts(time_instant); xploty_rts(time_instant)], Crts{time_instant}(1:2, 1:2), 'showMean', true, 'color', 'b', 'labels', ["x (m)", "y (m)"])

    legend("Estimated Location w/o Smoothing","Estimated Location w/ Smoothing")
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);