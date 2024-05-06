%% Coloured noise
clc
drift = 1;

M = 1000; % Number of time points. Including t = 0
N = 999; % Inner space points. N + 2 with boundary 
x_l = 0;
x_r = 1;
dx = (x_r-x_l)/(N + 1);
Thetas = [1];
U = zeros(N + 2,M, length(Thetas));
c = 1/(pi - 2); %  c = drift * dt / dx^2 = 1/(pi -2) by the paper
dt = dx^2 * c / drift;
T = M * dt; % Stopping time is not decided by us, still sol. is self similar

x_points = linspace(x_l, x_r, N + 2);
time_points = linspace(0, T, M);

aL = 1/2;
dim = 1;
gam = 2^(dim-aL)*pi^(dim/2)*gamma((dim-aL)/2)*gamma(aL/2);
%Q = integral2(f, x_points(3), x_points(5), x_points(1), x_points(2))
%integ(x_points(3), x_points(5), x_points(1), x_points(2))

alpha = 1/2;
% ANTAG a < b < c < d
Cov_disjoint = @(a,b,c,d) gam * (alpha * (alpha + 1))^(-1) * ( abs(c-b)^(alpha+1) - abs(d-b)^(alpha + 1) - abs(c-a)^(alpha + 1) + abs(d-a)^(alpha +1));
Cov_var = @(a,b) gam * 2 /(alpha*(alpha + 1)) * (b-a)^(alpha  + 1);

Q = zeros(N,N);
for i = 1:N
    for j = 1:N
        if i ~= j
            Q(i,j) = Cov_disjoint(x_points(i), x_points(i + 1), x_points(j), x_points(j + 1));
        else
            Q(i,i) = Cov_var(x_points(i), x_points(i + 1));
        end
    end
end
Q = dt * Q; % lol
Q_cells = {};
z = normrnd(0, ones(N*M,1));
%saved_randoms = load('randoms.mat');
%z = saved_randoms.randoms;
%z = z(:);
R = chol(Q); % behöver bara chol på varje block
F = zeros(1,N*M);
for m = 1:M
    F((m-1)*N + 1: m*N) = R*z((m-1)*N + 1: m*N);
end
F = reshape(F, [N,M]);
%%
for i = 1:M
    Q_cells{end + 1} = Q;
end

K = blkdiag(Q_cells{:});
%R = chol(K);
z = normrnd(0, ones(N*M,1));
F = K*z;

%%
BC = zeros(N, M);
BC(1, :) = U(1, :, 1);
BC(end, :) = U(end, :, 1);


for l = 1:length(Thetas)

r_1 = alpha * dt * Thetas(l) /(dx^2);
r_2 = alpha * dt * (1 - Thetas(l)) / (dx^2);

A = diag((1+2*r_1)*ones(1,N)) + diag(-r_1*ones(1,N-1),1) + ... 
    diag(-r_1*ones(1,N-1),-1);

A2 = diag((1-2*r_2)*ones(1,N)) + diag(r_2*ones(1,N-1),1) + ... 
    diag(r_2*ones(1,N-1),-1);

for m = 1:M-1
    m
    b = A2*U(2:end-1, m, l) + F(:, m)/(dx) + r_1*BC(:, m + 1) + r_2*BC(:, m);
    U(2:end-1, m + 1, l) = A\b;
end
end

%%
close all
figure
h = surf(time_points, x_points, U)
set(h,'LineStyle','none')

figure
plot(time_points, U(round(N/2), :));
xlabel('Time')
title('t \mapsto u(x,t)','FontSize', 16)

figure
plot(x_points, U(:,M/2))
xlabel('Space')
title('x \mapsto u(x,t)','FontSize', 16)
%% 
close all
x_point = round(N/2);plot
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
xlabel('Time');
ylabel('x');
zlabel('U');

%% quadratic variation
c = 1; % lång
mu = 4/(1+alpha);
variation = c*2^((1-alpha)/(2+alpha-1))*


