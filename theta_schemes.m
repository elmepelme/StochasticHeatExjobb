%% Chong and Walsh Finite Difference and Estimation of Drift
% u_t - alpha * u_xx = \dot{W}, x \in [x_l,x_r], t > 0
% u(0,t) = u(1,t) = 0
clc
alpha = exp(1)*pi^2;

M = 1000; % Number of time points. Including t = 0
N = 999; % Inner space points. N + 2 with boundary 
x_l = -2*sin(2);
x_r = pi+exp(1);
dx = (x_r-x_l)/(N + 1);

Thetas = [0.5];

U = zeros(N + 2,M, length(Thetas));
c = 1/(pi - 2); %  c = alpha * dt / dx^2 = 1/(pi -2) by the paper
dt = dx^2 * c / alpha;
T = M * dt; % Stopping time is not decided by us, still sol. is self similar

x_points = linspace(x_l, x_r, N + 2);
u0 = @(x) x.*(1-x)/2;
%U(:, 1) = u0(x_points);

W = sqrt(dt*dx)*normrnd(zeros(N, M), 1);
BC = zeros(N, M);
BC(1, :) = U(1, :, 1);
BC(end, :) = U(end, :, 1);
g = @(x,t) sin(x) + sin(t);
time_points = linspace(0, T, M);
%G = g(x_points(2:end-1), time_points)';
%g = 0;


%%

for l = 1:length(Thetas)

r_1 = alpha * dt * Thetas(l) /(dx^2);
r_2 = alpha * dt * (1 - Thetas(l)) / (dx^2);

A = diag((1+2*r_1)*ones(1,N)) + diag(-r_1*ones(1,N-1),1) + ... 
    diag(-r_1*ones(1,N-1),-1);

A2 = diag((1-2*r_2)*ones(1,N)) + diag(r_2*ones(1,N-1),1) + ... 
    diag(r_2*ones(1,N-1),-1);

for m = 1:M-1
    m
    b = A2*U(2:end-1, m, l) + W(:, m)/(dx) + r_1*BC(:, m + 1) + r_2*BC(:, m);
    U(2:end-1, m + 1, l) = A\b;
end
end
%%

sum = 0;
for l = 1:length(Thetas)
u_time = U(round(N/2), :, l);
for m = 1:M-1
    sum = sum + (u_time(m + 1) - u_time(m))^4;
end
quartic_variation(l) = sum;
drift_param_real = alpha;
drift_param_est(l) = 3/(pi*sum)*T;
sum = 0;
for j = 1:N
    sum = sum + (U(j + 1, M, l) - U(j, M, l)).^2;
end
quadratic_variation(l) = sum;
end

%% 
q2 = @(c, th) 0.5 ./ (sqrt(1+2*c*(2*th - 1)));
q4 = @(c, th, time) 3.*c.*time.*( (1-2.*th) ./ (sqrt(1+2*c*(2.*th - 1))) + 2.*th / sqrt(1+4.*c.*th)).^2;
real_quadratic_01 = q2(c, Thetas);
real_quadratic_min = q2(c, Thetas)*(x_r-x_l)/alpha;
real_quartic = q4(c, Thetas, T);
real_quartic_min = q4(c, Thetas, T)/(alpha);

%%
quadratic_variation
real_quadratic_min
quartic_variation 
real_quartic_min

%% Monte carlo ?
Alfas = zeros(N-1,1);
for i = 2:N
    i
    u_time_i = U(i, :);
    sum = 0;
    for m = 1:M-1
        sum = sum + (u_time_i(m+1) - u_time_i(m))^4;
    end
    Alfas(i) = 3/(pi*sum)*T;
end
mean(Alfas)
%%
close all
figure
h = surf(time_points, x_points, U)
set(h,'LineStyle','none')

figure
plot(time_points, u_time);
%% 
close all
x_point = round(N/2);
figure
axes = zeros(1, length(Thetas))
for i = 1:length(Thetas)
    axes(i) = subplot(2,2,i)
    plot(time_points, U(x_point, :, i))
    titleString = sprintf('theta = %0.2f, alpha = %0.3f, 4-var = %0.5f', Thetas(i), drift_param_est(i), quartic_variation(i));
    title(titleString);
end
linkaxes(axes(2:length(Thetas)));

%%
plot(x_points, U(:, round(M/2), 2))
%%
close all
h = surf(time_points, x_points, U)
set(h,'LineStyle','none')

%% quadratic variation

