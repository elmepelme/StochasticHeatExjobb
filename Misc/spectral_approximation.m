%% Spectral method of computation
% [0,L]
L = 1;
n = 1000;
x_points = linspace(0,L,n);
theta = 1;
K = 1000;
Phi = zeros(K, n);
lambda = zeros(1,K);
factor = zeros(K,n);
t_0 = 0.5;
for k = 1:K
    lambda(k) = k^2 * pi^2 * theta / (L^2);
    for i = 1:n
        Phi(k,i) = sqrt(2/L) * sin(pi*k*i/n);
        factor(k,i) = Phi(k,i)*((1- exp(-2*lambda(k)*t_0))^(1/2)) / ((2*lambda(k))^(1/2));
    end
end
Z = normrnd(zeros(1,K), 1);
u = Z*factor;

