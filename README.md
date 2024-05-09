# Simulations and inference of drift parameter for the stochastic heat equation

The stochastic heat eq. we are considering is the formal eq. (1)

u_t - \alpha * u_xx = \sigma * \dot{W}      
u(x,0) = u_0                                        
u(t,a) = u(t,b) = 0 // boundary if on x \in [a,b]

t > 0, x \in \mathbb{R} or x \in [a,b]

W is white noise as done by Walsh. 

Simulations based on covariance as well as finite diff schemes. Focus is convergence of variations.
See:
Chong, Y., Walsh, J.B. The Roughness and Smoothness of Numerical Solutions to the Stochastic Heat Equation. Potential Anal 37, 303–332 (2012). https://doi.org/10.1007/s11118-011-9257-6 https://personal.math.ubc.ca/~walsh/yuxiang.pdf
Estimation of the drift parameter for the fractional stochastic heat equation via power variation
Zeina Mahdi Khalil, Ciprian Tudor https://www.arxiv.org/abs/1912.07917

Etc.
