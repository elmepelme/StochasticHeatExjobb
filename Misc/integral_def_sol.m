function [u] = heat_sol_brute(t_points, x_points, dy, dt, D, T)
% The solution to the heat equation in point t and x
% u(t,x) = int_0^t int_R g(s, y, t, x) dW(s,y)
% Approximated by "brute force" with definition of the "Walsh-esque" 
% multivariate stochastic integral w.r.t Brownian sheet W(s,y)
% BUT NOTE DON'T THINK WORKS!
% Input:
% t points and x points to be evaluated
% dx and dt step size
% D diameter of space, "should be" = inf since integral over R for space
mu_A = dy * dt;

g = @(t, x, s, y) exp(-((x - y).^2) ./ (4*(t - s))) ./ ((4*pi*(t - s)).^(1/2));


y_nbr_points = 2*D/(dy);
y_points = linspace(-D,D,y_nbr_points);
t_nbr_points = size(t_points,2);
u = zeros(size(t_points,2), size(x_points,2));
white_noise = sqrt(mu_A)*normrnd(zeros(t_nbr_points, y_nbr_points), 1);

for i = 2:size(t_points,2)
    i
    for j = 1:size(x_points,2)
        u_xt = 0;
        for k = 2:(i-1)
            for l = 2:y_nbr_points
            s_point = t_points(k);
            y_point = y_points(l);
            u_xt = u_xt + g(t_points(i),x_points(j), s_point, y_point) * white_noise(k,l);
            end
        end
        u(i,j) = u_xt;
    end
end
end


