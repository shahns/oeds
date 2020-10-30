%%Script Excercise 7-PF: Nimish Shah s2088894
close all; clear; clc;

%% load data file
load('z_yacht2.mat');
stop = length(fi0);
num_meas = length(z);
%% initializations
number_of_particles_list = [32 200, 1200, 5000, 10000, 20000];%1000;
for i = 1:length(number_of_particles_list)
    var_book{i} = nan(8,num_meas); %for dyn samples
end
for iter = 1:length(number_of_particles_list) %for dyn samples
    number_of_particles = number_of_particles_list(iter); % for dyn samples
    clear residual; % for dyn samples
    % system model
    u = [400 * ones(1, length(fi0)); fi0];
    number_of_states = 8;
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
    invCn = inv(Cn);
    
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
    
    
    x_pred = [mu_x_xi, mu_x_xi, mu_x_v, ...
        mu_x_v, mu_x_a, mu_x_a, mu_x_t, mu_x_phi]';
    
    
    K_eff = NaN(1, length(z));
    test_var = NaN(size(z));
    x_est_mean = NaN(number_of_states, length(z));
    cov_est = NaN(number_of_states, number_of_states, length(z));
    
    % particle_log = NaN(number_of_states, number_of_particles, stop+1);
    
    %% particle filter: generation of samples
    x_samples = x_pred + sqrt(diag(Cx0)) .* randn(number_of_states, number_of_particles); % as randn has unit var
    % particle_log(:, :, 1) = x_samples;
    
    %% particle filter: condensation algorithm
    for i = 1:stop
        isMeasAvailable = find(imeas == i);
        if isMeasAvailable % we do have a meas
            fprintf('Measurement available. Meas#%i\n', isMeasAvailable);
            z_samples = hmeas(x_samples, x0, Cn); %%is Cn reqd??
            residual(1:2:3, :) = wrapTo180(z(1:2:3, isMeasAvailable)-z_samples(1:2:3, :));
            residual(2, :) = z(2, isMeasAvailable) - z_samples(2, :);
            W = exp(-0.5*sum(residual.*(invCn * residual)))'; % weights
            if (sum(W) == 0)
                error('process did not converge');
            end
            W = W / sum(W);
            CumW = cumsum(W);
            
            test_var(:, isMeasAvailable) = sum((z_samples < z(:, isMeasAvailable)), 2) / number_of_particles;
            x_est_mean(:, i) = x_samples(:, :) * W; % Sample mean
            
            dummy_cov = zeros(8);
            for p = 1:number_of_particles
                dummy_cov = dummy_cov + W(p) * (x_samples(:, p) - x_est_mean(:, i)) * (x_samples(:, p) - x_est_mean(:, i))';
            end
            cov_est(:, :, i) = dummy_cov;
            var_book{iter}(:,isMeasAvailable) = diag(dummy_cov);% for dyn samples
            % Find an index permutation using golden rule root finding
            for j = 1:number_of_particles
                R = rand;
                ja = 1;
                jb = number_of_particles;
                while (ja < jb - 1)
                    jx = floor(jb-0.382*(jb - ja));
                    fa = R - CumW(ja);
                    fb = R - CumW(jb);
                    fxx = R - CumW(jx);
                    if (fb * fxx < 0), ja = jx;
                    else, jb = jx;
                    end
                end
                ind(j) = jb;
            end
            
            % Resample
            for j = 1:number_of_particles
                Ys(:, j) = x_samples(:, ind(j));
            end
            %         particle_log(:, :, i+1) = Ys;
            K_eff(isMeasAvailable) = 1 ./ (sum(W.^2)); % effective particles
        else % we do NOT have a meas
            Ys = x_samples;
            x_est_mean(:, i) = x_samples(:, :) * W; % Sample mean
            dummy_cov = zeros(8);
            for p = 1:number_of_particles
                dummy_cov = dummy_cov + W(p) * (x_samples(:, p) - x_est_mean(:, i)) * (x_samples(:, p) - x_est_mean(:, i))';
            end
            cov_est(:, :, i) = dummy_cov;
        end
        x_samples = fsys(Ys, u(:, i), Cw);
    end
    fprintf('##########complete#########%i\n', iter);
end %for dyn samples
xRange = 1:length(x_est_mean(1, :));
%% Plot for path and uncertainity ellipse
figure;
set(gcf,'renderer','Painters')
plot(x_est_mean(1, :), x_est_mean(2, :), 'LineWidth', 2);
hold on
plot(x_est_mean(1, 1), x_est_mean(2, 1), 'r>');
plot(x_est_mean(1, end), x_est_mean(2, end), 'r*');

for p = 1:length(z)
    plot_cov_ellipse(x_est_mean(1:2, p), cov_est(1:2, 1:2, p));
end

legend(["Est. Path", "Start Point", "End Point", "Uncertainty Ellipse"], 'Location', 'southwest');
xlabel('x in meters');
ylabel('y in meters');

%% plot for pos vs. time
figure; %xi
plot(xRange, x_est_mean(1, :), xRange, x_est_mean(2, :), 'LineWidth', 2);
legend('\xi_x', '\xi_y');
xlabel('time');
ylabel('\xi in m');
axis([0 max(xRange) -Inf inf])
%% plot for vel vs. time
figure; %v
plot(xRange, x_est_mean(3, :), xRange, x_est_mean(4, :), 'LineWidth', 2);
legend('v_x', 'v_y');
xlabel('time');
ylabel('v in ms^{-1}');
axis([0 max(xRange) -Inf inf])
%% plot for acc vs. time
figure; %a
plot(xRange, x_est_mean(5, :), xRange, x_est_mean(6, :), 'LineWidth', 2);
legend('a_x', 'a_y');
xlabel('time');
ylabel('a in ms^{-2}');
axis([0 max(xRange) -Inf inf])
%% plot for t vs. time
figure; %t
plot(xRange, x_est_mean(7, :), 'LineWidth', 2);
xlabel('time');
ylabel('Thrust in N');
axis([0 max(xRange) -Inf inf])
%% plot for phi vs. time
figure; %phi
plot(xRange, x_est_mean(8, :), 'LineWidth', 2);
xlabel('time');
ylabel('\phi in degrees');
axis([0 max(xRange) -Inf inf])
%% Plot for test variables
figure;
set(gcf,'renderer','Painters')
plot(1:length(test_var(1,:)), test_var(1,:),'r.')
hold on
plot(1:length(test_var(2,:)), test_var(2,:),'b.')
plot(1:length(test_var(3,:)), test_var(3,:),'g.')
% plot(1:length(z), test_var(1, :), 'r.', 1:length(z), test_var(2, :), 'b.', 1:length(z), test_var(3, :), 'g.');
xlabel('Measurement instance (m)');
ylabel('Test variable value');
legend('U_1(i_m)', 'U_2(i_m)', 'U_3(i_m)')

%% Plot for K_eff
figure;
set(gcf,'renderer','Painters')
stem(K_eff,':+');
xlabel('Measurement instance (m)');
ylabel('Effective Particles');

%%
% close all
% AniLaTeX = [];
% frameContent = '';
% figureDirectory = 'animation-200'; % dir name
% figcounter = 0;
% filename = 'frame'; % file name
% frmps = 6;
%
% mkdir(figureDirectory);
% frames = 100; %number of frames
% for i = 1:150:stop
%
%     figure(1),
%     plot(x_est(1, :), x_est(2, :), 'LineWidth', 2);
%     hold on
%     plot(x_est(1, 1), x_est(2, 1), 'r>');
%     plot(x_est(1, end), x_est(2, end), 'r*');
%
%
%     if sum(i == imeas) == 1
%         plot(particle_log(1, :, i), particle_log(2, :, i), 'g.')
%     else
%         plot(particle_log(1, :, i), particle_log(2, :, i), 'r.')
%     end
%     hold off
%
%     set(gcf, 'PaperPositionMode', 'auto')
%
%     figcounter = figcounter + 1;
%     print(gcf, '-depsc2', '-loose', [figureDirectory, '/', filename, 'No', num2str(figcounter), '.eps']);
%     AniLaTeX = animategraphicsLaTeX(AniLaTeX, [filename, 'No', num2str(figcounter), '.eps'], figureDirectory, frameContent, frmps);
%     %
% end


% figure(1), close;
% animategraphicsLaTeX(AniLaTeX, [figureDirectory, '/', filename], 'controls=all,loop,width=3.5in');
%
% clear AniLaTeX filename figureDirectory includegraphicsOptions frmps;

%%
