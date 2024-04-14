# Simulations and inference of drift parameter for the stochastic heat equation

The stochastic heat eq. we are considering is the formal eq.

u_t - \alpha * u_xx = \dot{F}      
u(x,0) = u_0                                        (1)
u(t,a) = u(t,b) = 0 // boundary if on x \in [a,b]

t > 0, x \in \mathbb{R} or x \in [a,b]
\dot{F} is some noise process, either space-time white, or white in time and homog. coloured in space
E(W(t,A)W(s,B)) = t \land s * \lambda(A \cap B) , t,s time, A,B measurable subsets in space
E(F(t,A)F(s,B)) = t \land s * \int_A \int_B f(x-y)dxdy where f is Fourier transform of tempered measure

This rep. is some simulations of the eq. on [a,b], studying the behaviour of finite-difference \Theta schemes 
Solutions to (1) are almost 1/4 Hölder-continuous in time and almost 1/2 Hölder-cont in space and 
they therefore have nontrivial quadratic- and quartic variations in time and space respectively.
The \Theta schemes don't always give correct limiting values.
This is interesting when we want to look at inference tools for the drift param. \alpha in (*)
