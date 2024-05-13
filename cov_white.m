%% Inference SHE drift and diffussion
clear all
clc
diffusion = exp(1);
drift = pi;
abs_moment = @(v, s) s^v * 2^(v/2) * gamma(v/2 + 1/2) / sqrt(pi); %% E( \vert Z \vert^v ), Z \in N(0,\sigma^2)
r_x = @(t, x1, x2) 1/sqrt(drift) * diffusion^2*sqrt(t/(2*pi*drift)) * ... 
    exp((-(x1-x2)^2) / (8*drift*t)) + diffusion^2*((x1-x2)/ ... 
    (sqrt(drift)*4*sqrt(drift)))*erf((x1-x2)/(2*sqrt(2*t*drift))) - ...
    diffusion^2 * abs(x1-x2)/(4 * sqrt(drift)* sqrt(drift));
t = 1;

r_t = @(t1,t2) diffusion^2/sqrt(4*pi*drift) * (sqrt(t1 + t2) - sqrt(abs(t1 - t2)));

N = 5000;
x_L = 0;
x_R = 1;
t_points = linspace(0, t, N);
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
C_t = C_t(2:end, 2:end); % Remove zeros
fprintf('Done. Computing Choletsky...')
R_x = chol(C_x);
R_t = zeros(N,N);
R_t(2:N, 2:N) = chol(C_t);
%chol(C_t);

fprintf('Done. Please save matrices')
%%
%
K = 10000;
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



%%
close all
clc
two_variations = sum(diff(U_x).^2);
four_variations = sum(diff(U_t).^4);


%drift_est_time = 3*diffusion^4 * (t-t_points(2)) ./ (pi * four_variations);
drift_est_time =  (t) * diffusion^4 * abs_moment(4,1) ./ (pi * four_variations);
drift_est_space = diffusion^2 * (x_R - x_L) ./ (2 * two_variations);

diffusion_est_space = (2*drift.*two_variations ./ (x_R - x_L)).^(1/2);
%diffusion_est_time = ((pi * drift * four_variations) ./ (3 * (t-t_points(2)))).^(1/4);
diffusion_est_time = ((pi * drift * four_variations) ./ (abs_moment(4,1) * (t))).^(1/4);

fprintf('mean of drift_est_time, mean of drift_est_space ,mean of diffusion_est_time, and mean of diffusion_est_time.')
mean(drift_est_time)
mean(drift_est_space)
mean(diffusion_est_time)
mean(diffusion_est_space)

%fprintf('std of drift_est_time, std of drift_est_space ,std of diffusion_est_time, and std of diffusion_est_time.')

%std(drift_est_time)
%std(drift_est_space)
%std(diffusion_est_time)
%std(diffusion_est_space)

%%
clc
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
s2_4 = rho2_4 * drift^(2) * abs_moment(4,1)^(-2) / (N-1)
s2_4_est = std(drift_est_time)^2
s2_2 = rho2_2 * drift^2 * abs_moment(2,1)^(-2) / N
s2_2_est = std(drift_est_space)^2

s2_diff_4 = rho2_4 / (16) * diffusion^2 * abs_moment(4,1)^(-2) / (N-1)
s2_diff_4_est = std(diffusion_est_time)^2

s2_diff_2 = rho2_2 / (4) * diffusion^2 * abs_moment(2,1)^(-2) / (N)
s2_diff_2_est = std(diffusion_est_space)^2
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
x = linspace(min(diffusion_est_space), max(diffusion_est_space), 10000);
pdf_normal = normpdf(x, diffusion, sqrt(s2_diff_2));
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
