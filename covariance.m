clear all
clc
C_u =  @(x,t,y,s) (1/sqrt(2*pi))* (sqrt(t + s)*exp((-abs(x-y)^2)/(2*(t+s))) ...
- sqrt(abs(t - s))*exp(-(abs(x-y)^2)/(2*abs(t-s)))) ...
+ (x-y)*(cdf('Normal', (x-y)/sqrt(s+t), 0, 1) - cdf('Normal', (x-y)/sqrt(abs(t-s)), 0, 1));
C_u_var = @(x,t,y,s) 1/sqrt(2*pi) * sqrt(t + s);
C_u_nan = @(x,t,y,s) (1/sqrt(2*pi)) * sqrt(t + s)*exp((-abs(x-y)^2)/(2*(t+s)));

%
%       b * * * b m = 3 = M
%       b * * * b m = 2
%       b * * * b m = 1
%   n = 0 1 2 3 4 
%              (4 = N + 1)
N = 3;
M = 3;
a = 0;
b = 1;
T = 1;
x_points = linspace(a, b, N);
t_points = linspace(0, T, M + 1);
t_points = t_points(2:end);

% Compute the size of the covariance matrix
N = length(x_points);
M = length(t_points);

% Initialize the covariance matrix
cov_matrix = zeros(N, M, N, M);

for i = 1:N
    for j = 1:M
        for k = 1:N
            for l = 1:M
                kov = C_u(x_points(i), t_points(j), x_points(k), t_points(l));
                if isnan(kov)
                    cov_matrix(i,j,k,l) = C_u_var(x_points(i), t_points(j), x_points(k), t_points(l));
                else
                    cov_matrix(i,j,k,l) = kov;
                end
            end
        end
    end
end

cov_matrix_2D = reshape(cov_matrix, [N*M, N*M]);

disp('Covariance Matrix:')
disp(cov_matrix_2D)
