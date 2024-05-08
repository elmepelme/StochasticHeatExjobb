%% Simulating SHE on [x_L,x_R]
% u_t - drift * u_xx = \dot{F}
% u(x_L,t) = u(x_R,t) = 0
% u(x,0) = u_0(x)

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

x_points = linspace(x_L, x_R, N + 2);
t_points = linspace(0, T, M);

%%%%% Initial Conditions %%%%%%%
u0 = @(x) 10*(x_L-x).*(x_R-x);
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


gam = 1/2;
dim = 1;
c_riesz = 2^(dim-gam)*pi^(dim/2)*gamma((dim-gam)/2)*gamma(gam/2); % Constant
Cov_disjoint = @(a,b,c,d) c_riesz * (gam * (gam + 1))^(-1)  ... 
    * (abs(c-b)^(gam+1) - abs(d-b)^(gam + 1) - abs(c-a)^(gam + 1) ... 
    + abs(d-a)^(gam +1));

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

%% Solving the systems using one-step \Theta finite differences
Theta = 0.25; % Finite-Diff theta

U_White_Noise = zeros(N + 2,M);
U_Coloured_Noise = zeros(N + 2, M);

r_1 = drift * dt * Theta /(dx^2);
r_2 = drift * dt * (1 - Theta) / (dx^2);

A1 = diag((1+2*r_1)*ones(1,N)) + diag(-r_1*ones(1,N-1),1) + ... 
    diag(-r_1*ones(1,N-1),-1);

A2 = diag((1-2*r_2)*ones(1,N)) + diag(r_2*ones(1,N-1),1) + ... 
    diag(r_2*ones(1,N-1),-1);

for m = 1:M-1
    b_White_Noise = A2*U_White_Noise(2:end-1, m) + W(:, m)/(dx);
    b_Coloured_Noise = A2*U_Coloured_Noise(2:end-1, m) + F(:, m)/(dx);
    
    U_White_Noise(2:end-1, m + 1) = A1\b_White_Noise;
    U_Coloured_Noise(2:end-1, m + 1) = A1\b_Coloured_Noise;

end

%% Plots
close all
figure;

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
drift_param_real = drift;
drift_param_est = 3/(pi*sum)*T;
sum = 0;
for j = 1:N
    sum = sum + (U_White_Noise(j + 1, M) - U_White_Noise(j, M)).^2;
end
quadratic_variation = sum;

q2 = @(c, th) 0.5 ./ (sqrt(1+2*c*(2*th - 1)));
q4 = @(c, th, time) 3.*c.*time.*( (1-2.*th) ./ (sqrt(1+2*c*(2.*th - 1))) + 2.*th / sqrt(1+4.*c.*th)).^2;
real_quadratic_01 = q2(c, Theta);
real_quadratic = q2(c, Theta)*(x_R-x_L)/drift;
real_quartic = q4(c, Theta, T);
real_quartic = q4(c, Theta, T)/(drift);


quadratic_variation
real_quadratic
quartic_variation 
real_quartic


%% 

save("Theta 0.25.mat", "U_White_Noise", "quadratic_variation", "real_quadratic", "quartic_variation", "real_quartic", "drift_param_est", "drift_param_real")

%% Loading saved files 

theta025 = load("Theta 0.25.mat");
theta1 = load("Theta 1.mat");
theta05 = load("Theta 05.mat");

%%
u025 = theta025.U_White_Noise;
u1 = theta1.U_White_Noise;
u05 = theta05.U_White_Noise;
M = 1000;
u025space = u025(:, M);
u1space = u1(:, M);
u05space = u05(:, M);
figure
plot(u1space)
hold on
plot(u05space)
hold on
plot(u025space)


%%
figure;

plot(x_points, u05space, 'LineWidth', 1.5);  % Theta = 1
%hold on;
%plot(x_points, u05space, 'LineWidth', 1.5); % Theta = 0.5
%plot(x_points, u025space, 'LineWidth', 1.5); % Theta = 0.25

%egend('\Theta = 1', '\Theta = 0.5', '\Theta = 0.25', 'FontSize', 20);

title('$x \mapsto u(x,t)$ for $t > 0$','FontSize', 20, 'interpreter','latex')
xlabel('x');
ylabel('u(x,t)');

grid on;
grid minor;
box on;
set(gca, 'FontSize', 20);
set(gca, 'LineWidth', 1.5);

%% surfs
U_White_Noise_Solution = u05; % theta = 0.5
h = surf(t_points, x_points, U_White_Noise_Solution);
grid off;
set(h, 'LineStyle', 'none');
title('$u(t,x)$ Simulation', 'Interpreter', 'latex');
xlabel('Time');
ylabel('Space');
view(45, 45);
box on;
set(gca, 'FontSize', 20);
set(gca, 'LineWidth', 2);

%%
figure;
u05time = u05(round(N/2), :);
plot(t_points, u05time, 'LineWidth', 1.5);  % Theta = 0.5
%hold on;
%plot(x_points, u05space, 'LineWidth', 1.5); % Theta = 0.5
%plot(x_points, u025space, 'LineWidth', 1.5); % Theta = 0.25

%egend('\Theta = 1', '\Theta = 0.5', '\Theta = 0.25', 'FontSize', 20);

title('$t \mapsto u(x,t)$ for $t > 0$','FontSize', 20, 'interpreter','latex')
xlabel('t');
ylabel('u(x,t)');

grid on;
grid minor;
box on;
set(gca, 'FontSize', 20);
set(gca, 'LineWidth', 1.5);
