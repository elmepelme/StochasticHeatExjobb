%% Inference SHE drift and diffussion using paths from covariance
clear all
clc
diffusion = pi;
drift = exp(1);

abs_moment = @(v) 2^(v/2) * gamma(v/2 + 1/2) / sqrt(pi); %% E( \vert Z \vert^v ), Z \in N(0,\sigma^2)

r_x = @(t, x1, x2) 1/sqrt(drift) * diffusion^2*sqrt(t/(2*pi*drift)) * ... 
    exp((-(x1-x2)^2) / (8*drift*t)) + diffusion^2*((x1-x2)/ ... 
    (sqrt(drift)*4*sqrt(drift)))*erf((x1-x2)/(2*sqrt(2*t*drift))) - ...
    diffusion^2 * abs(x1-x2)/(4 * sqrt(drift)* sqrt(drift));

r_t = @(t1,t2) diffusion^2/sqrt(4*pi*drift) * (sqrt(t1 + t2) - sqrt(abs(t1 - t2)));

N = 1000;
x_L = 0;
x_R = 1;
t = 1;
s = 0;
t_points = linspace(s, t, N);
x_points = linspace(x_L,x_R, N);

% Covariance matrices
C_x = zeros(N,N);
C_t = zeros(N,N);

fprintf('Computing covariance function...')
for i = 1:N
    i
    for j = 1:N
        C_x(i,j) = r_x(t, x_points(i), x_points(j));
        C_t(i,j) = r_t(t_points(i), t_points(j));
    end
end
C_t = C_t(2:end, 2:end); % Remove zeros %
fprintf('Done. Computing Choletsky...')
R_x = chol(C_x);
R_t = zeros(N,N);
R_t(2:N, 2:N) = chol(C_t);
%chol(C_t);

fprintf('Done. Please save matrices')

%% Simulating independent paths
K = 1000;
U_x = zeros(N, K);
U_t = zeros(N, K);
for k = 1:K
    k
    z = normrnd(zeros(N,1),1);
    U_x(:, k) = R_x * z;
    U_t(:, k) = R_t * z;
end
U_x = U_x(2:end, :); % första värdet är konstigt
fprintf('Done')

% kan bara skala om variationen så slipper köra om...
%%
clc
% TESTA ATT ÄNDRA PÅ DRIFT OCH DIFFUSION RANDOMLY
if drift == 1 && diffusion == 1
    fprintf('lol')
    drift = 2;
    diffusion = 1;
    factor = diffusion*drift^(-1/4);
    two_variations = factor^2*drift^(-1/2)*sum(diff(U_x).^2);
    four_variations = factor^4*sum(diff(U_t).^4);
else
    two_variations = sum(diff(U_x).^2);
    four_variations = sum(diff(U_t).^4);
end
%%

drift_est_time =  (t) * diffusion^4 * abs_moment(4) ./ (pi * four_variations);
drift_est_space = diffusion^2 * (x_R - x_L) ./ (2 * two_variations);

diffusion_est_time = ((pi * drift * four_variations) ./ (abs_moment(4) * (t))).^(1/4);
diffusion_est_space = (2*drift.*two_variations ./ (x_R - x_L)).^(1/2);

fprintf('Mean of drift_est_time: %.4f\n', mean(drift_est_time));
fprintf('Mean of drift_est_space: %.4f\n', mean(drift_est_space));
fprintf('Mean of diffusion_est_time: %.4f\n', mean(diffusion_est_time));
fprintf('Mean of diffusion_est_space: %.4f\n', mean(diffusion_est_space));

rho2_term = @(H, q, z) factorial(q)/(2^q) * (abs(z + 1)^(2*H) + abs(z-1)^(2*H) - 2*abs(z)^(2*H))^q;
D = 1000;
sum1 = 0;
sum2 = 0;
for z = -D:D
    sum1 = sum1 + rho2_term(1/4,4,z);
    sum2 = sum2 + rho2_term(1/2,2,z);
end
rho2_4 = sum1;
rho2_2 = sum2;%rho2_term(1/2, 2,0);
%

s2_4 = rho2_4 * drift^(2) * abs_moment(4)^(-2) / (N); % Real variance drift
s2_4_est = std(drift_est_time)^2 ;% Estimated variance drift
s2_2 = rho2_2 * drift^2 * abs_moment(2)^(-2) / N ;% Real variance drift
s2_2_est = std(drift_est_space)^2 ;% Estimated variance drift

s2_diff_4 = rho2_4 / (16) * diffusion^2 * abs_moment(4)^(-2) / (N); % Real variance diffusion
s2_diff_4_est = std(diffusion_est_time)^2; % estimated variance diffusion
s2_diff_2 = rho2_2 / (4) * diffusion^2 * abs_moment(2)^(-2) / (N); % Real variance diffusion  
s2_diff_2_est = std(diffusion_est_space)^2; % Estimated variance diffusion

fprintf('Drift - time Real Variance (s2_4): %.4f\n', s2_4);
fprintf('Drift - time Estimated variance (s2_4_est): %.4f\n', s2_4_est);
fprintf('Drift - space Real Variance (s2_2): %.4f\n', s2_2);
fprintf('Drift - space Estimated Variance (s2_2_est): %.4f\n', s2_2_est);
fprintf('Diffusion - time Real variance (s2_diff_4): %.4f\n', s2_diff_4);
fprintf('Diffusion - time Estimated variance (s2_diff_4_est): %.4f\n', s2_diff_4_est);
fprintf('Diffusion - space Real variance (s2_diff_2): %.4f\n', s2_diff_2);
fprintf('Diffusion - space Estimated variance (s2_diff_2_est): %.4f\n', s2_diff_2_est);
s2_4_est/s2_4
s2_diff_4_est/s2_diff_4

sig_2 = 0;
sig_4 = 0;
sig_2_space = 0;
sig_4_space = 0;
D = 1000;
r_l = @(k, H) 1/2 * (abs(k + 1)^(2*H) + abs(k - 1)^(2*H) - 2*abs(k)^(2*H));
for i = 1:D
    for j = 1:D
        sig_2 = sig_2 + r_l(abs(i-j), 1/4)^2;
        sig_4 = sig_4 + r_l(abs(i-j), 1/4)^4;
        sig_2_space = sig_2_space + r_l(abs(i-j), 1/2)^2;
        sig_4_space = sig_4_space + r_l(abs(i-j), 1/2)^4;
    end
end
sig = (72*sig_2 + 24*sig_4)/D;
sig_space = (2*sig_2_space)/D;

s2_diff_4 = sig / (16) * diffusion^2 * abs_moment(4)^(-2) / (N);
s2_4 = sig * drift^(2) * abs_moment(4)^(-2) / (N); % Real variance drift
fprintf('Drift - time Real Variance (s2_4): %.4f\n', s2_4);
fprintf('Drift - time Estimated variance (s2_4_est): %.4f\n', s2_4_est);
fprintf('Drift - space Real Variance (s2_2): %.4f\n', s2_2);
fprintf('Drift - space Estimated Variance (s2_2_est): %.4f\n', s2_2_est);
fprintf('Diffusion - time Real variance (s2_diff_4): %.4f\n', s2_diff_4);
fprintf('Diffusion - time Estimated variance (s2_diff_4_est): %.4f\n', s2_diff_4_est);
fprintf('Diffusion - space Real variance (s2_diff_2): %.4f\n', s2_diff_2);
fprintf('Diffusion - space Estimated variance (s2_diff_2_est): %.4f\n', s2_diff_2_est);
s2_4_est/s2_4
s2_diff_4_est/s2_diff_4

%%
close all
figure;
% MAN KAN KÖRA HISTFIT OCKSÅ! DÅ KÖR DEN PÅ STD OCH MEAN AV DATAN

num_bins = round(sqrt(K));

subplot(2, 2, 1);
histogram(drift_est_time, num_bins, 'Normalization', 'pdf');
hold on
x = linspace(min(drift_est_time), max(drift_est_time), 10000);
pdf_normal = normpdf(x, drift, sqrt(s2_4));
plot(x, pdf_normal, 'r', 'LineWidth', 2);
title('Drift Estimate - Time', 'FontSize', 14);
xlabel('Value', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
legend('Histogram', 'Asymptotic Normal Distribution')

subplot(2, 2, 2);
histogram(drift_est_space, num_bins, 'Normalization', 'pdf');
hold on
x = linspace(min(drift_est_space), max(drift_est_space), 10000);
pdf_normal = normpdf(x, drift, sqrt(s2_2));
plot(x, pdf_normal, 'r', 'LineWidth', 2);
title('Drift Estimate - Space', 'FontSize', 14);
xlabel('Value', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
legend('Histogram', 'Asymptotic Normal Distribution')


subplot(2, 2, 3);
histogram(diffusion_est_time, num_bins, 'Normalization', 'pdf');
hold on
x = linspace(min(diffusion_est_time), max(diffusion_est_time), 10000);
pdf_normal = normpdf(x, diffusion, sqrt(s2_diff_4));
plot(x, pdf_normal, 'r', 'LineWidth', 2);
title('Diffusion Estimate - Time', 'FontSize', 14);
xlabel('Value', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
legend('Histogram', 'Asymptotic Normal Distribution')

subplot(2, 2, 4);
histogram(diffusion_est_space, num_bins, 'Normalization', 'pdf');
hold on
x = linspace(min(diffusion_est_space), max(diffusion_est_space), 10000);
pdf_normal = normpdf(x, diffusion, sqrt(s2_diff_2));
plot(x, pdf_normal, 'r', 'LineWidth', 2);
title('Diffusion Estimate - Space', 'FontSize', 14);
xlabel('Value', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
legend('Histogram', 'Asymptotic Normal Distribution')

sgtitle(sprintf(['Asymptotic Normality Illustration for' ...
    ' Estimates\nDrift: $\\alpha$ = %.2f, Diffusion: ' ...
    '$\\sigma$ = %.2f'], drift, diffusion), 'FontSize', 30, ...
    'Interpreter', 'latex');
%%

close all
u_x = U_x(:, randi(K));
u_t = U_t(:, randi(K));

subplot(2, 1, 1);
plot(t_points, u_t, 'LineWidth', 1);
xlabel('Time')
title('$t \mapsto \sigma u_{\alpha}(t,x)$','FontSize', 25, 'interpreter','latex')
grid on;
xlim([min(t_points), max(t_points)]);

xpoints = x_points(2:end);
subplot(2, 1, 2);
plot(xpoints, u_x, 'LineWidth', 1);
xlabel('Space')
title('$x \mapsto \sigma u_{\alpha}(t,x)$','FontSize', 25, 'interpreter','latex')
grid on;
xlim([min(xpoints), max(xpoints)]);

title_str = sprintf('$\\sigma u_{\\alpha}(t,x)$ for $\\alpha = %g$, $\\sigma = %g$', drift, diffusion);
sgtitle(title_str, 'FontSize', 30, 'interpreter','latex');


