%% Coloured noise
clear
clc
a = 1/2;
d = 1;
gam = 2^(d-a)*pi^(d/2)*gamma((d-a)/2)*gamma(a/2);
f = @(x,y) gam * abs(x-y).^(a-d);
N = 10;
x_points = linspace(0,1, N + 2);
Q = integral2(f, x_points(1), x_points(2), x_points(3), x_points(4))
integ = @(a,b,c,d) -4/3 * ((d-b)^(3/2) + (a-d)*sqrt(d-a) + (b-c)*sqrt(c-b) + (c-a)^(3/2)); % RÄTT TROR JAG!!
gam*integ(x_points(1), x_points(2), x_points(3), x_points(4))