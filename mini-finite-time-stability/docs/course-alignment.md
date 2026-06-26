# Course Alignment -- Finite-Time Stability
## MIT 6.241J
- **Finite-time convergence**: x(t)->0 in finite time T(x0)<inf
- **Homogeneity**: Key tool for finite-time stability analysis
- Connection: fts_homogeneity_degree computes degree
## Stanford AA203
- **Terminal controllers**: Finite-time optimal control
- **Sliding mode**: Finite-time reaching phase
## Key Formulas
dV/dt <= -c*V^alpha, 0<alpha<1 => finite-time stable
T <= V(0)^{1-alpha}/(c*(1-alpha))
Homogeneous degree: f(lambda^{ri}*xi) = lambda^{k+ri}*fi(x)
## Exercises
1. Verify finite-time stability via Lyapunov function
2. Estimate settling time bound