%% SHE simulations using distributions
clc
clear all
clc
drift = 1;

% FIXA X LED
gamma = 1/2; % Riesz
c = 1/sqrt(2*pi) * (2*pi/drift)^(1/4);
f = @(x) abs(x).^(-gamma).*exp(-abs(x).^2) * (1 / gamma);
Q = integral(f,-Inf,Inf);
cov_white =  @(x,t,y,s) (1/sqrt(2*pi))* (sqrt(t + s)*exp((-abs(x-y)^2)/(2*(t+s))) ...
+ (x-y)*(cdf('Normal', (x-y)/sqrt(s+t), 0, 1) - 1));

c_Coloured = (2*drift)^((gamma  - 1)/4) * ((2 / pi) * Q)^(1/2);

T = 1;
M = 1000;
N = 1000;
dt = T/M;
t_points = linspace(0, T, M);
t_points = t_points(2:end);
x_points = linspace(0, 1, N);
R = @(t,s, H, K) (1/(2^K)) * ( (t^(2*H) + s^(2*H))^K - abs(t-s)^(2*H*K));
C_R = zeros(M - 1,M - 1);
C_R_c = zeros(M - 1, M - 1);
C_R_white = zeros(N - 1, N - 1);
H_w = 1/2;
K_w = 1/2;
H_c = 1/2;
K_c = 1 + (gamma - 1)/2;
for i = 1:M-1
    for j = 1:M-1
        C_R(i, j) = c^2 * R(t_points(i), t_points(j), H_w, K_w);
        C_R_c(i, j) = c_Coloured^2 * R(t_points(i), t_points(j), H_c, K_c);
    end
end

for i = 1:N -1
    for j = 1:N -1
        C_R_white = cov_white(0.5,0.5,x_points(i),x_points(j));
    end
end
L = chol(C_R);
L_white_x = chol(C_R_white);
L_c = chol(C_R_c);
Z = normrnd(zeros(M - 1,1),1);
Z_x = normrnd(zeros(N - 1, 1), 1);
u_w = L * Z;
u_w_x = L_white_x * Z_x;
u_c = L_c * Z;

%%
sum = 0;
sum_c = 0;
sum_c_finite = 0;
sum_white_x = 0;
for i = 2:M-1
    sum = sum + (u_w(i) - u_w(i - 1))^4;
    sum_c = sum_c + (u_c(i) - u_c(i - 1))^4;
    sum_c_finite = sum_c_finite + (u_finite_c(round(N/2), i) - u_finite_c(round(N/2), i - 1))^4;
end

%%

for i = 2:N-1
   sum_white_x = sum_white_x + (u_w_x(i) - u_w_x(i - 1))^2;
end
sum
sum_c
sum_c_finite
sum_white_x
true_variat = 3*T/(pi*drift)


%%
close all
figure;
plot(t_points, u_w, 'b', 'LineWidth', 1);
hold on;
plot(t_points, u_c, 'r', 'LineWidth', 1);
xlabel('Time t', 'FontSize', 14);
ylabel('Solution u(x,t)', 'FontSize', 14);
title('Solutions to SHE as a function of time', 'FontSize', 16);
legend('White Noise Solution', 'Coloured Noise Solution', 'FontSize', 12);
grid on;
hold off;

