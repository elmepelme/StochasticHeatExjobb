%% Distribution in space
clc
% t1 = t2
R_White_x = @(t, x1, x2) sqrt(t/(2*pi)) * exp((-(x1-x2)^2) / (8*t))+ ...
    ((x1-x2)/4)*( erf( (x1-x2)/(2*sqrt(2*t)) ) - 1);
N = 1000;
x_points = linspace(0, 1, N);
R_x_matrix = zeros(N,N);

t = 1;
for i = 1:N
    for j = 1:N
        R_x_matrix(i,j) = R_White_x(1, x_points(i), x_points(j));
    end
end
R_x = chol(R_x_matrix);

z = normrnd(zeros(N,1),1);
u_x = R_x * z;
sum = 0;
for i = 2:N
    sum = sum + (u_x(i) - u_x(i-1))^2;
end
sum
%% Distribution in time
N = 1000;
t_points = linspace(0, 5, N);
R_White_t = @(t1,t2) 1/sqrt(4*pi) * (sqrt(t1 + t2) - sqrt(abs(t1 - t2)));
R_t_matrix = zeros(N-1,N-1);



for i = 2:N
    for j = 2:N
        R_t_matrix(i-1,j-1) = R_White_t(t_points(i), t_points(j));
    end
end
R_t = chol(R_t_matrix);
z = normrnd(zeros(N-1,1),1);
u_t = R_t * z;
sum = 0;
for i = 2:N-1
    sum = sum + (u_t(i) - u_t(i-1))^4;
end
sum
%%
plot(t_points(2:end), u_t)
