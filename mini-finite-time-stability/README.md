# Finite-Time Stability

## Core Theory
Finite-time stability ensures convergence to equilibrium in finite time (not asymptotically). Key conditions: Lyapunov dV/dt <= -c*V^alpha with 0<alpha<1, or negative homogeneity degree.

## API (100+ functions)
fts_create, fts_step_rk4, fts_integrate, fts_lyap_create, fts_lyap_analyze, fts_lyap_settling_time_bound, fts_hom_weights_create, fts_hom_is_finite_time, fts_st_compute, fts_st_empirical, fts_smc_create, fts_smc_simulate, fts_ft_params_create, fts_ft_settling_bound

## Build
```bash
make && make test
```

## Theory Details

### Lyapunov Condition for Finite-Time Stability

A system is finite-time stable if there exists V(x) >= 0 and constants c>0, 0<alpha<1:
```
dV/dt <= -c * V(x)^alpha
```
Then the settling time is bounded by: T <= V(x0)^(1-alpha) / (c*(1-alpha)).

### Homogeneity-Based Approach

A vector field f(x) is homogeneous of degree d < 0 with respect to weights (r1,...,rn) if:
```
f_i(lambda^r1*x1, ..., lambda^rn*xn) = lambda^(d+ri) * f_i(x)
```
If the origin of x_dot = f(x) is asymptotically stable and f has negative homogeneity degree, then the origin is finite-time stable.

### Sliding Mode Control (SMC)

First-order SMC with discontinuous control drives trajectories to a sliding surface s(x)=0 in finite time, then maintains sliding motion. Key properties: robustness to matched uncertainties, chattering in implementation.

Super-twisting algorithm: second-order SMC that provides finite-time convergence with continuous control, avoiding chattering while preserving robustness.

### Fixed-Time Stability (Polyakov, 2012)

A stronger property: the settling time is bounded uniformly for ALL initial conditions:
```
dV/dt <= -(c1*V^alpha + c2*V^beta),  with 0<alpha<1<beta
```
Then T_max <= 1/(c1*(1-alpha)) + 1/(c2*(beta-1)), independent of x0.

### Comparison of Stability Types

| Stability Type | Convergence | Settling Time Bound |
|----------------|-------------|---------------------|
| Asymptotic | t->infinity | Unbounded |
| Exponential | ~exp(-lambda*t) | Depends on x0 |
| Finite-time | T_finite < infinity | T <= V0^(1-alpha)/(c*(1-alpha)) |
| Fixed-time | T_fixed < infinity | Uniform bound independent of x0 |

### Applications

- **Terminal sliding mode control**: Sliding surface with fractional power ensures finite-time convergence after reaching
- **Consensus in multi-agent systems**: Finite-time consensus protocols for rapid coordination
- **Observation**: Finite-time observers (e.g., Levant's differentiator) for robust state estimation
- **Guidance**: Impact-time-control guidance laws for simultaneous target interception

### Implementation

| File | Purpose |
|------|---------|
| `fts_core.h/c` | Finite-time stability verification, settling time bounds |
| `fts_lyapunov.h/c` | Lyapunov-based FTS detection, V^alpha condition check |
| `fts_homogeneity.h/c` | Homogeneity degree computation, weight determination |
| `fts_sliding_mode.h/c` | SMC and super-twisting algorithm simulation |
| `fts_fixed_time.h/c` | Fixed-time stability analysis (Polyakov conditions) |
| `fts_settling.h/c` | Settling time estimation via numerical and analytical methods |

### Key Functions

| Function | Description |
|----------|-------------|
| `fts_lyap_analyze` | Check V_dot <= -c*V^alpha condition |
| `fts_settling_time_bound` | Compute upper bound on settling time |
| `fts_hom_is_finite_time` | Test negative homogeneity condition |
| `fts_smc_simulate` | Simulate sliding mode controller |
| `fts_st_empirical` | Measure settling time from simulation data |

### Examples

Example files demonstrate: finite-time convergence verification, homogeneity-based detection, and sliding mode controller simulation with chattering analysis.

### References
- Bhat, S.P. & Bernstein, D.S. (2000). Finite-time stability of continuous autonomous systems. SIAM J. Control Optim., 38:751-766.
- Polyakov, A. (2012). Nonlinear feedback design for fixed-time stabilization. IEEE TAC, 57:2106-2110.
- Moulay, E. & Perruquetti, W. (2006). Finite time stability and stabilization: State of the art. Lecture Notes in Control and Info. Sci., 334:23-41.
- Levant, A. (2003). Higher-order sliding modes, differentiation and output-feedback control. Int. J. Control, 76:924-941.
- Utkin, V., Guldner, J. & Shi, J. (2009). Sliding Mode Control in Electromechanical Systems. CRC Press.


### Examples

Example programs demonstrate:
1. **example1**: Finite-time convergence verification via V^alpha condition
2. **example2**: Homogeneity-based detection with weight computation
3. **example3**: Sliding mode controller simulation with chattering analysis

### Building

```
make && make test && make examples
```
