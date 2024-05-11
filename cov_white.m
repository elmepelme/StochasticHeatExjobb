%% Inference SHE drift and diffussion
clear all
clc
diffusion = exp(1);
drift = pi;

r_x = @(t, x1, x2) 1/sqrt(drift) * diffusion^2*sqrt(t/(2*pi*drift)) * ... 
    exp((-(x1-x2)^2) / (8*drift*t)) + diffusion^2*((x1-x2)/ ... 
    (sqrt(drift)*4*sqrt(drift)))*erf((x1-x2)/(2*sqrt(2*t*drift))) - ...
    diffusion^2 * abs(x1-x2)/(4 * sqrt(drift)* sqrt(drift));
t = 1;

r_t = @(t1,t2) diffusion^2/sqrt(4*pi*drift) * (sqrt(t1 + t2) - sqrt(abs(t1 - t2)));

N = 3000;
x_L = 0;
x_R = 1;
t_points = linspace(0, t, N);
x_points = linspace(x_L,x_R, N);

% Covariance matrices
C_x = zeros(N,N);
C_t = zeros(N,N);

for i = 1:N
    for j = 1:N
        C_x(i,j) = r_x(t, x_points(i), x_points(j));
        C_t(i,j) = r_t(t_points(i), t_points(j));
    end
end
C_t = C_t(2:end, 2:end); % Remove zeros

R_x = chol(C_x);
R_t = chol(C_t);

%
K = 1000;
U_x = zeros(N, K);
U_t = zeros(N - 1, K);
for k = 1:K
    z = normrnd(zeros(N,1),1);
    U_x(:, k) = R_x * z;
    U_t(:, k) = R_t * z(2:end);
end
U_x = U_x(2:end, :); % första värdet är konstigt

%%
close all
clc
two_variations = sum(diff(U_x).^2);
four_variations = sum(diff(U_t).^4);

abs_moment = @(v, s) s^v * 2^(v/2) * gamma(v/2 + 1/2) / sqrt(pi); %% E( \vert Z \vert^v ), Z \in N(0,\sigma^2)

drift_est_time = 3*diffusion^4 * t ./ (pi * four_variations);
drift_est_space = diffusion^2 * (x_R - x_L) ./ (2 * two_variations);

diffusion_est_space = (2*drift.*two_variations ./ (x_R - x_L)).^(1/2);
diffusion_est_time = ((pi * drift * four_variations) ./ (3 * t)).^(1/4);

fprintf('mean of drift_est_time, mean of drift_est_space ,mean of diffusion_est_time, and mean of diffusion_est_time.')
mean(drift_est_time)
mean(drift_est_space)
mean(diffusion_est_time)
mean(diffusion_est_space)

fprintf('std of drift_est_time, std of drift_est_space ,std of diffusion_est_time, and std of diffusion_est_time.')

std(drift_est_time)
std(drift_est_space)
std(diffusion_est_time)
std(diffusion_est_space)


%%
close all
figure;
num_bins = round(sqrt(K));

subplot(2, 2, 1);
histfit(drift_est_time, num_bins, 'normal');
title('Drift Estimate - Time', 'FontSize', 14);
xlabel('Value', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);

subplot(2, 2, 2);
histfit(drift_est_space, num_bins, 'normal');
title('Drift Estimate - Space', 'FontSize', 14);
xlabel('Value', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);

subplot(2, 2, 3);
histfit(diffusion_est_space, num_bins, 'normal');
title('Diffusion Estimate - Space', 'FontSize', 14);
xlabel('Value', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);

subplot(2, 2, 4);
histfit(diffusion_est_time, num_bins, 'normal');
title('Diffusion Estimate - Time', 'FontSize', 14);
xlabel('Value', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);

sgtitle(sprintf(['Normality Illustration for' ...
    ' Estimates\nDrift: $\\alpha$ = %.2f, Diffusion: ' ...
    '$\\sigma$ = %.2f'], drift, diffusion), 'FontSize', 30, ...
    'Interpreter', 'latex');
%%

close all
u_x = U_x(:, randi(K));
u_t = U_t(:, randi(K));
u_t = [0; u_t];

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
