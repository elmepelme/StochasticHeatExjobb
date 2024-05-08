%% Inference SHE drift and diffussion
 % Covariance Functions for u with drift = 1, diffusion = 1
r_x = @(t, x1, x2) sqrt(t/(2*pi)) * exp((-(x1-x2)^2) / (8*t))+ ...
    ((x1-x2)/4)*erf((x1-x2)/(2*sqrt(2*t))) - abs(x1-x2)/4;
t = 1;

r_t = @(t1,t2) 1/sqrt(4*pi) * (sqrt(t1 + t2) - sqrt(abs(t1 - t2)));

N = 3000;
points = linspace(0,1, N);

% Covariance matrices
C_x = zeros(N,N);
C_t = zeros(N,N);

for i = 1:N
    for j = 1:N
        C_x(i,j) = r_x(t, points(i), points(j));
        C_t(i,j) = r_t(points(i), points(j));
    end
end
C_t = C_t(2:end, 2:end); % Remove zeros

R_x = chol(C_x);
R_t = chol(C_t);

%%
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
two_variations = sum(diff(U_x).^2);
four_variations = sum(diff(U_t).^4);

mean(two_variations)
mean(four_variations)

%% 
close all
u_x = U_x(:, randi(K));
u_t = U_t(:, randi(K));
points = linspace(0,1, N);
points = points(2:end);

figure;
subplot(2, 1, 1);
plot(points, u_t);
xlabel('Time')
title('$t \mapsto u(t,x)$','FontSize', 16, 'interpreter','latex')

subplot(2, 1, 2);
plot(points, u_x);
xlabel('Space')
title('$x \mapsto u(t,x)$','FontSize', 16, 'interpreter','latex')

sgtitle('$u(t,x)$ for $\alpha = \sigma = 1$', 'FontSize', 20, 'interpreter','latex');
