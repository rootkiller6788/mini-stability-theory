# Lyapunov Stability Theory

## Lyapunov Direct Method (1892)

A Lyapunov function V(x) generalizes the concept of energy. If V decreases along trajectories, the system is stable. The method avoids solving differential equations.

## Lyapunov Equation

For x_dot = Ax, V(x)=x'Px is a Lyapunov function iff A'P + PA = -Q for some Q>0. Equivalently: A is Hurwitz iff the Lyapunov equation has a unique positive definite solution P.

## Construction Methods

- Krasovskii: V = f(x)'Pf(x) with J'P+PJ<0
- Variable gradient: assume grad_V = g(x), enforce curl(g)=0
- SOS: polynomial V via semidefinite programming

## Converse Theorems

Massera (1949) proved that asymptotic stability implies the existence of a smooth Lyapunov function, justifying the Lyapunov approach as a necessary and sufficient method.
