# Region of Attraction Estimation

## Lyapunov Sublevel Set Method
For V(x) with V'(x) < 0, any sublevel set {x: V(x) <= c}
contained in the region where V' < 0 is an invariant subset
of the ROA.

## Algorithm
1. Choose Lyapunov function V
2. Find maximum c such that V'(x) < 0 for V(x) <= c
3. Use binary search to determine max safe level
4. Monte Carlo verification for confidence

## Convergence Rate
Exponential convergence: |x(t)| <= C*exp(-alpha*t)*|x(0)|
alpha estimated from Lyapunov function decay rate.
