%% Simulating SHE on [x_L,x_R] (White Noise) Working!
% u_t - drift * u_xx = diffusion * \dot{W} 
% u(t, x_L) = u(t,x_R) = 0
% u(0,x) = u_0(x)
clear all
clc
%%%%%%% Defining step sizes and matrices %%%%%%
drift = 1;
diffusion = 1;

abs_moment = @(v)  2^(v/2) * gamma(v/2 + 1/2) / sqrt(pi); %% E( \vert Z \vert^v ), Z \in N(0,1)

M = 1000; % Time points
N = 999;  % Inner Space points, N + 2 points including boundary

x_L = 0;
x_R = 1;

c = 1/(pi - 2); % CFL Number c = drift * dt / (dx^2)
dx = (x_R - x_L) / (N + 1);
dt = c/drift * (dx^2); 
T = M * dt; % Stopping time is not decided by us, still sol. is self similar
%T = 0.05;
%dt = T/M;
%drift * dt / (dx^2)
%
x_points = linspace(x_L, x_R, N + 2);
t_points = linspace(0, T, M);

%%%%% Initial Conditions %%%%%%%
%u0 = @(x) sin(2*pi*x / (x_R - x_L));
%U_White_Noise(1, :) = u0(x_points);

% Solving the systems using one-step \Theta finite differences
Theta = [1]; % Finite-Diff theta
%%%%%%%% Init. Noise Field %%%%%%%%%%%
K = 1000; % How many fields
quartic_variations = zeros(length(Theta),K);
quadratic_variations = zeros(length(Theta),K);

for k = 1:K
    disp("Itter: "+ num2str(k))
    Z = normrnd(zeros(M, N), 1);

    %%%%% White Noise
    W = diffusion * sqrt(dt*dx) * Z;
    U_White_Noise = zeros(M, N + 2, length(Theta));
    %U_White_Noise(1, :, :) = u0(x_points);
    for l = 1:length(Theta) % For each choice of theta

        r_1 = drift * dt * Theta(l) /(dx^2);
        r_2 = drift * dt * (1 - Theta(l)) / (dx^2);

        A = diag(2*ones(1,N)) + diag(-1*ones(1,N-1),1) + ... 
            diag(-1*ones(1,N-1),-1);
        Ar2 = (eye(N)-r_2*A);
        Ar1 = (eye(N) + r_1*A);

        for m = 1:M-1
    
            b = Ar2*U_White_Noise(m, 2:end - 1, l)' + W(m,:)'/(dx);
            U_White_Noise(m + 1, 2:end-1, l) = (Ar1\b)';
        end

    u_t = U_White_Noise(:, round((N + 2)/2), l);
    u_x = U_White_Noise(M, :, l);
    
    % Variations for White Noise Sol.
    sum_4 = 0;

    u_time = u_t;
    for m = 1:M-1
        sum_4 = sum_4 + (u_time(m + 1) - u_time(m))^4;
    end
    quartic_variations(l,k) = sum_4;
    sum_2 = 0;
    for j = 1:N
        sum_2 = sum_2 + (u_x(j + 1) - u_x(j)).^2;
    end
    quadratic_variations(l,k) = sum_2;
    end
end
%
close all
figure('units','normalized','outerposition',[0 0 1 1])
h = surf(x_points, t_points, U_White_Noise(:,:,1));
%grid off;
%set(h, 'LineStyle', 'none');
grid on
set(h, 'LineStyle', 'none')
colormap(spring);
title('Temperature Distribution $u(t,x)$', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Space', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Time', 'FontSize', 20, 'FontWeight', 'bold');
zlabel('Temperature', 'FontSize', 20, 'FontWeight', 'bold');
view(35, 20);
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on', 'TickDir', 'out');
ax = gca;
ax.XColor = 'black';
ax.YColor = 'black';
ax.ZColor = 'black';
lighting phong;
%shading interp;
%light('Position', [1, 1, 1], 'Style', 'infinite');
set(gcf, 'Color', 'w');
%%
% Estimations
abs_moment = @(v)  2^(v/2) * gamma(v/2 + 1/2) / sqrt(pi); %% E( \vert Z \vert^v ), Z \in N(0,1)
fprintf('drift time and space')
drift_est_time = 3*diffusion^4 * T ./ (pi * quartic_variations);
drift_est_space = diffusion^2 * (x_R - x_L) ./ (2 * quadratic_variations);

fprintf('diffusion time and space')
diffusion_est_time = ((pi * drift * quartic_variations) ./ (3 * T)).^(1/4);
diffusion_est_space = (2*drift.*quadratic_variations ./ (x_R - x_L)).^(1/2);

% joint est

alpha_hat = quartic_variations * pi * abs_moment(2)^2*(x_R - x_L)^2 ./ ...
    (4 * quadratic_variations.^2*abs_moment(4) * T);

sigma_sq_hat = quartic_variations * pi * abs_moment(2)*(x_R - x_L) ./ ...
    (2 * quadratic_variations*abs_moment(4) * T);
sigma_hat = sqrt(sigma_sq_hat);
%
%%
clc
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

s2_4 = rho2_4 * drift^(2) * abs_moment(4)^(-2) / (M); % Real variance drift
s2_4_est = std(drift_est_time)^2 ;% Estimated variance drift
s2_2 = rho2_2 * drift^2 * abs_moment(2)^(-2) / N ;% Real variance drift
s2_2_est = std(drift_est_space)^2 ;% Estimated variance drift

s2_diff_4 = rho2_4 / (16) * diffusion^2 * abs_moment(4)^(-2) / (M); % Real variance diffusion
s2_diff_4_est = std(diffusion_est_time)^2; % estimated variance diffusion
s2_diff_2 = rho2_2 / (4) * diffusion^2 * abs_moment(2)^(-2) / (N); % Real variance diffusion  
s2_diff_2_est = std(diffusion_est_space)^2; % Estimated variance diffusion


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

u_t = U_White_Noise(:, round((N + 2)/2), 1);
u_x = U_White_Noise(M, :, 1);

figure;
subplot(2, 1, 1);
plot(t_points, u_t);
xlim([0 T])
xlabel('Time')
title('$t \mapsto u(t,x)$','FontSize', 16, 'interpreter','latex')
grid on;
grid minor;
box on;
set(gca, 'FontSize', 20);
set(gca, 'LineWidth', 1.5);

subplot(2, 1, 2);
plot(x_points, u_x);
xlabel('Space')
title('$x \mapsto u(t,x)$','FontSize', 16, 'interpreter','latex')
title_str = sprintf('$\\sigma u_{\\alpha}(t,x)$ for $\\alpha = %g$, $\\sigma = %g$ and $\\Theta = 0.5$', drift, diffusion);
sgtitle(title_str, 'FontSize', 30, 'interpreter','latex');

grid on;
grid minor;
box on;
set(gca, 'FontSize', 20);
set(gca, 'LineWidth', 1.5);

%%
close all
h = surf(T - t_points, x_points, flipud(fliplr(U_White_Noise(:,:,1)))');
grid off;
set(h, 'LineStyle', 'none');
title('$u(t,x)$ Simulation', 'Interpreter', 'latex');
xlabel('Time');
ylabel('Space');
%view(45, 45);
box on;
set(gca, 'FontSize', 20);
set(gca, 'LineWidth', 2);

%% Plots (if 3 diff thetas)
close all
figure
plot(t_points, u_times(:, 1));
hold on;
plot(t_points, u_times(:, 2));
plot(t_points, u_times(:, 3));
hold off;
xlabel('Time');
ylabel('u(t,x)');
title('$t \mapsto u(t,x)$ for $\Theta = 0.25, 0.5$ and 1','FontSize', 16, 'interpreter','latex')
%title_str = sprintf('$\\sigma u_{\\alpha}(t,x)$ for $\\alpha = %g$, $\\sigma = %g$', drift, sigma);
%sgtitle(title_str, 'FontSize', 30, 'interpreter','latex');

grid on;
grid minor;
box on;
set(gca, 'FontSize', 20);
set(gca, 'LineWidth', 1.5);
legend('$\Theta = 0.25$', '$\Theta = 0.5$', '$\Theta = 1$','FontSize', 16, 'interpreter','latex');

figure
plot(x_points, U_White_Noise(M, :, 1));
hold on;
plot(x_points, U_White_Noise(M, :, 2));
plot(x_points, U_White_Noise(M, :,  3));
hold off;
xlabel('Space');
ylabel('u(t,x)');
title('$x \mapsto u(t,x)$ for $\Theta = 0.25, 0.5$ and 1','FontSize', 16, 'interpreter','latex')
%title_str = sprintf('$\\sigma u_{\\alpha}(t,x)$ for $\\alpha = %g$, $\\sigma = %g$', drift, sigma);
%sgtitle(title_str, 'FontSize', 30, 'interpreter','latex');

grid on;
grid minor;
box on;
set(gca, 'FontSize', 20);
set(gca, 'LineWidth', 1.5);
legend('$\Theta = 0.25$', '$\Theta = 0.5$', '$\Theta = 1$','FontSize', 16, 'interpreter','latex');



%%
q2 = @(c, th) 0.5 ./ (sqrt(1+2*c*(2*th - 1)));
q4 = @(c, th, time) 3.*c.*time.*( (1-2.*th) ./ (sqrt(1+2*c*(2.*th - 1))) + 2.*th ./ sqrt(1+4.*c.*th)).^2;
real_quadratic_sim_01 = q2(c, Theta);
real_quadratic_sim = q2(c, Theta)*(x_R-x_L)/(drift) * diffusion^2 % Notera ändring med sigma
real_quartic_sim = q4(c, Theta, T);
real_quartic_sim = q4(c, Theta, T)/(drift) * diffusion^4 % Notera ändring med sigma
