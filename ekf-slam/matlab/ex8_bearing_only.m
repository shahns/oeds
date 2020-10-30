%%Script Excercise 8 Part 1: Nimish Shah s2088894
close all; clear; clc;

%% load data file
load('SLAM.mat');

%% initializations
v = u(1, :);
phi = u(2, :);
no_iters = size(u, 2);
F = eye(2);

xest{1} = zeros(1, 1);
Cest{1} = zeros(1);
G{1} = zeros(2);

std_meas_noise_mtr = 5;
std_new_landmark_mtr = 22;

Cw_tilde = diag([std_velocity^2, std_heading^2]);


xpred = zeros(2, no_iters);
no_visible_landmarks = 0;

%% Part IV

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
        r = 22;
        [~, ~, current_new_lms_id] = intersect(current_new_lms, current_meas.id, 'stable');
        for lm = 1: no_current_new_lm
            xpred{time_instant} = [xpred{time_instant}; [r*cosd(current_meas.zbearing(current_new_lms_id(lm))); r*sind(current_meas.zbearing(current_new_lms_id(lm)))]+xpred{time_instant}(1:2)]; 
        end
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
        z = current_meas.zbearing;
        Cv = eye(LMBOOK.total_visible(time_instant)) * 2^2;
        [zpred, Hjacobian] = hmeas_bearing_only(xpred{time_instant}, current_visible_lm_state_ind);
        S = Hjacobian * Cpred{time_instant} * Hjacobian' + Cv;
        K = Cpred{time_instant} * Hjacobian' / S;
        
        xest{time_instant} = xpred{time_instant} + K * wrapTo180(z - zpred);
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
    xplotxx(time_instant) = xest{time_instant}(1);
end

for time_instant = 1:no_iters
    xplotyy(time_instant) = xest{time_instant}(2);
end
plot(xplotxx, xplotyy, '.', 'Color', '#0072BD');

for time_instant = 1:100:no_iters
    plot_cov_ellipse([xplotxx(time_instant); xplotyy(time_instant)], Cest{time_instant}(1:2, 1:2))
end
xlabel("x(m)");
ylabel("y(m)");
legend("Estimated Location", 'Location', 'northwest')

%% performance plots

eigen_slam = nan(2, no_iters);
for i = 1:no_iters
    eigen_slam(:, i) = eig(Cest{i}(1:2, 1:2));
    sqrt_max_eigen_slam(i) = sqrt(max(eigen_slam(:)));
end

figure;
hold on
plot(1:20:no_iters, sqrt_max_eigen_slam(1:20:end), '-x')
xlabel("time instant")
ylabel("sqrt. of max eigen value (m)")
legend("SLAM w/o smoothing", 'Location', 'northwest');
%% video
v = VideoWriter('bearing-only-loop-closing.avi', 'Motion JPEG AVI');
open(v);
figure;

for time_instant=1:no_iters
    clf
    title('Bearing Only SLAM')
    ylim([-150,350])
    xlim([-50,450])
    hold on
    plot(xplotx(1:time_instant), xploty(1:time_instant),'.', 'Color', '#0072BD')
    plot_cov_ellipse([xplotx(time_instant); xploty(time_instant)], Cest{time_instant}(1:2, 1:2), 'showMean', true, 'color', 'b', 'labels', ["x (m)", "y (m)"])

    legend("Estimated Location", 'Location', 'northwest')
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);