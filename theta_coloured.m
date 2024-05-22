%% Simulating SHE on [x_L,x_R]
% u_t - drift * u_xx = \dot{F}
% u(x_L,t) = u(x_R,t) = 0
% u(x,0) = u_0(x)

clear all
clc
%%%%%%% Defining step sizes and matrices %%%%%%

drift = 1;

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

x_points = linspace(x_L, x_R, N + 2);
t_points = linspace(0, T, M);

%%%%% Initial Conditions %%%%%%%
u0 = @(x) 0.2*sin(2*pi*x / (x_R - x_L));


%u0 = @(x) 0.1*sin(2*pi*x/(x_R-x_L));
%U_White_Noise(:, 1) = u0(x_points);
%W_Coloured_Noise(:,1) = u0(x_points);

%%%%%%%% Init. Noise Field %%%%%%%%%%%

Z = normrnd(zeros(N, M), 1);
Z_flat = Z(:); % For coloured

%%%%% White Noise
W = sqrt(dt*dx)*Z;

%%%%% Coloured Noise
% Riesz Kernel Covariance 
F = zeros(1, N * M);

gam = 3/4;
dim = 1;
c_riesz = 2^(dim-gam)*pi^(dim/2)*gamma((dim-gam)/2)/gamma(gam/2); % Constant
Cov_disjoint = @(a,b,c,d) c_riesz * (gam * (gam + 1))^(-1)  ... 
    * (abs(c-b)^(gam+1) - abs(d-b)^(gam + 1) - abs(c-a)^(gam + 1) ... 
    + abs(d-a)^(gam +1));
abs_moment = @(v, s) s^v * 2^(v/2) * gamma(v/2 + 1/2) / sqrt(pi);

K = zeros(N,N);
for i = 1:N
    for j = 1:N
            K(i,j) = Cov_disjoint(x_points(i), x_points(i + 1), x_points(j), x_points(j + 1));
    end
end
K = dt * K;
R = chol(K);

for m = 1:M
    F((m-1)*N + 1: m*N) = R*Z_flat((m-1)*N + 1: m*N);
end
F = reshape(F, [N,M]);

% Solving the systems using one-step \Theta finite differences
Theta = 0.5; % Finite-Diff theta

U_White_Noise = zeros(N + 2,M);
U_Coloured_Noise = zeros(N + 2, M);
U_White_Noise(:, 1) = u0(x_points);
U_Coloured_Noise(:, 1) = u0(x_points);

r_1 = drift * dt * Theta /(dx^2);
r_2 = drift * dt * (1 - Theta) / (dx^2);

A1 = diag((1+2*r_1)*ones(1,N)) + diag(-r_1*ones(1,N-1),1) + ... 
    diag(-r_1*ones(1,N-1),-1);

A2 = diag((1-2*r_2)*ones(1,N)) + diag(r_2*ones(1,N-1),1) + ... 
    diag(r_2*ones(1,N-1),-1);

for m = 1:M-1
    if mod(m,100) == 0
        disp("Itter: "+ num2str(m))
    end
    b_White_Noise = A2*U_White_Noise(2:end-1, m) + W(:, m)/(dx);
    b_Coloured_Noise = A2*U_Coloured_Noise(2:end-1, m) + F(:, m)/(dx);
    
    U_White_Noise(2:end-1, m + 1) = A1\b_White_Noise;
    U_Coloured_Noise(2:end-1, m + 1) = A1\b_Coloured_Noise;
end
%%
close all
figure('units','normalized','outerposition',[0 0 1 1])
h1 = surf(x_points, t_points, U_Coloured_Noise(:,:,1)');
set(h1, 'LineStyle', 'none');
grid on
set(h1, 'LineStyle', 'none')
colormap(spring);
title('Temperature Distribution $u(t,x)$ - Coloured Noise', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Space', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Time', 'FontSize', 20, 'FontWeight', 'bold');
zlabel('$u(t,x)$', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold');
view(35, 20);
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on', 'TickDir', 'out');
ax1 = gca;
ax1.XColor = 'black';
ax1.YColor = 'black';
ax1.ZColor = 'black';
lighting phong;

figure('units','normalized','outerposition',[0 0 1 1])
h2 = surf(x_points, t_points, U_White_Noise(:,:,1)');
set(h2, 'LineStyle', 'none');
grid on
set(h2, 'LineStyle', 'none')
colormap(spring);
title('Temperature Distribution $u(t,x)$ - White Noise', 'Interpreter', 'latex', 'FontSize', 24, 'FontWeight', 'bold');
xlabel('Space', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Time', 'FontSize', 20, 'FontWeight', 'bold');
zlabel('$u(t,x)$', 'Interpreter', 'latex', 'FontSize', 20, 'FontWeight', 'bold');
view(35, 20);
set(gca, 'FontSize', 20, 'LineWidth', 2, 'Box', 'on', 'TickDir', 'out');
ax2 = gca;
ax2.XColor = 'black';
ax2.YColor = 'black';
ax2.ZColor = 'black';
lighting phong;

set(gcf, 'Color', 'w');


%% col variations (drift = 1 enl Tudor Bok från variation white coloured)
% s 159 prop 8.8
HK = 1/2*(1 - (dim-gam)/2); % H * K
f = @(xi) exp(-abs(xi).^2) / ... 
    (2*(-dim/2 + gam/2 + 1)) .* abs(xi).^(-gam);
C0 = sqrt(2^(1-dim/2 + gam/2)*integral(f, -Inf, Inf)); 
C = C0^(4/(2+gam-dim));
u_coloured_time = U_Coloured_Noise(round(N/2), :);
s = 0;
delta = 1;
for j = delta:length(u_coloured_time) - 1
    s = s + abs(u_coloured_time(j+1) - u_coloured_time(j)).^(1/HK);
end
s
true_var = sqrt(2^((dim-gam)/2))* C0 * T * abs_moment(1/HK,1) % * (2^(1-dim/2 + gam/2))^(1/HK) stämmer ej överens med white noise fallet!!!!
s/true_var
%% 
U_Coloured_Noise = temp;
%% Plots
close all
figure;
amplitude = (U_White_Noise(round(N/2), :) /U_Coloured_Noise(round(N/2), :));
temp = U_Coloured_Noise;
U_Coloured_Noise = amplitude .* U_Coloured_Noise;
%U_Coloured_Noise = temp;
subplot(2, 2, 1);
plot(t_points, U_White_Noise(round(N/2), :));
xlabel('Time')
title('$t \mapsto u_W(x,t)$','FontSize', 16, 'interpreter','latex')

subplot(2, 2, 2);
plot(x_points, U_White_Noise(:,M))
xlabel('Space')
title('$x \mapsto u_W(x,t)$','FontSize', 16, 'interpreter','latex')

subplot(2, 2, 3);
plot(t_points, U_Coloured_Noise(round(N/2), :));
xlabel('Time')
title('$t \mapsto u_C(x,t)$','FontSize', 16, 'interpreter','latex')

subplot(2, 2, 4);
plot(x_points, U_Coloured_Noise(:,M))
xlabel('Space')
title('$x \mapsto u_C(x,t)$','FontSize', 16, 'interpreter','latex')

sgtitle('White Noise and Coloured Noise Solutions for same Random Outcome', 'FontSize', 20);

figure;
subplot(2, 1, 1);

plot(t_points, U_White_Noise(round(N/2), :));
hold on; 
plot(t_points, U_Coloured_Noise(round(N/2), :));
hold off; 
xlabel('Time')
legend('White Noise Solution', 'Coloured Noise Solution')
title('$t \mapsto u(x,t)$ Time Domain','FontSize', 16, 'interpreter','latex')

subplot(2, 1, 2);

plot(x_points, U_White_Noise(:,round(M/2)));
hold on;
plot(x_points, U_Coloured_Noise(:,round(M/2)));
hold off; 
xlabel('Space')
legend('White Noise Solution', 'Coloured Noise Solution')
title('$x \mapsto u(x,t)$ Space Domain','FontSize', 16, 'interpreter','latex')

sgtitle('White Noise and Coloured Noise Solutions for same Random Outcome', 'FontSize', 20);

%% Variations for White Noise Sol.

sum = 0;

u_time = U_White_Noise(round(N/2), :);
for m = 1:M-1
    sum = sum + (u_time(m + 1) - u_time(m))^4;
end
quartic_variation = sum;
drift_param_real = drift
drift_param_est = 3/(pi*sum)*T
sum = 0;
for j = 1:N
    sum = sum + (U_White_Noise(j + 1, M) - U_White_Noise(j, M)).^2;
end
quadratic_variation = sum;

q2 = @(c, th) 0.5 ./ (sqrt(1+2*c*(2*th - 1)));
q4 = @(c, th, time) 3.*c.*time.*( (1-2.*th) ./ (sqrt(1+2*c*(2.*th - 1))) + 2.*th / sqrt(1+4.*c.*th)).^2;
real_quadratic_01 = q2(c, Theta);
real_quadratic = q2(c, Theta)*(x_R-x_L)/drift;

real_quartic = q4(c, Theta, T)/(drift);

quadratic_variation
real_quadratic
quartic_variation 
real_quartic
