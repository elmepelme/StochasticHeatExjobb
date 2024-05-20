%% simulate a bifractional
K = 1/4;
H = 3/4;
r_bfbm = @(t,s) 2^(-K) * ((t^(2*H) + s^(2*H))^K - abs(t-s)^(2*H*K));

N = 1000;
t = 1;
epsilon = 1e-5;
t_points = linspace(epsilon,t,N);

C_bfbm = zeros(N,N);

for i = 1:N
    for j = 1:N
        C_bfbm(i,j) = r_bfbm(t_points(i), t_points(j));
    end
end

R = chol(C_bfbm)';

J = 500;
Z = normrnd(zeros(N,J), 1);
U = R*Z;

%
abs_moment = @(v) 2^(v/2) * gamma(v/2 + 1/2) / sqrt(pi); 
close all
figure
plot(U(:, randi(J)))
q = 1/(H*K);
q_variation = mean(sum(abs(diff(U)).^q))
exact_q_variaton = sqrt(2^(1-K))^(1/(H*K))*(t-epsilon) * abs_moment(1/(H*K))
