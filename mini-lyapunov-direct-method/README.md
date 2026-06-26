# Lyapunov Direct Method

## The cornerstone of nonlinear stability analysis

Aleksandr Lyapunov's 1892 doctoral thesis established the "second method" for analyzing stability without solving differential equations. A Lyapunov function V(x) acts as a generalized energy -- if energy decreases along trajectories, the system is stable.

### Lyapunov's Main Theorems

**Theorem 1 (Stability)**: If V(x) > 0 for x != 0, V(0) = 0, and V_dot(x) <= 0 along trajectories of x_dot = f(x), then the origin is stable in the sense of Lyapunov.

**Theorem 2 (Asymptotic Stability)**: If additionally V_dot(x) < 0 for x != 0, the origin is asymptotically stable. All trajectories starting sufficiently close converge to the origin.

**Theorem 3 (Global Asymptotic Stability)**: If also V(x) -> infinity as |x| -> infinity (radially unbounded), stability is global -- all trajectories converge to the origin regardless of initial condition.

**Theorem 4 (Exponential Stability)**: If c1*|x|^2 <= V(x) <= c2*|x|^2 and V_dot(x) <= -c3*|x|^2, then |x(t)| <= sqrt(c2/c1)*|x0|*exp(-c3*t/(2*c2)).

### Lyapunov Equation for Linear Systems

For x_dot = Ax, choose V(x) = x'Px. Then V_dot = x'(A'P + PA)x = -x'Qx. The Lyapunov equation A'P + PA = -Q has a unique positive definite solution P for any Q > 0 if and only if A is Hurwitz (all eigenvalues have negative real parts).

This reduces stability verification to solving a linear matrix equation, computable via the Bartels-Stewart algorithm.

### Common Lyapunov Function Candidates

| System Type | V(x) Candidate | Note |
|-------------|----------------|------|
| Linear x_dot = Ax | x'Px | Solve Lyapunov equation for P |
| Mechanical | (1/2)mv^2 + Potential(x) | True physical energy |
| Gradient x_dot = -grad f(x) | f(x) | Natural Lyapunov function |
| Lur'e (linear + nonlinear) | x'Px + theta*integral(phi) | Popov/Lur'e type |
| Polynomial | Sum-of-Squares (SOS) | SOSTOOLS via SDP |

### Construction Methods

1. **Krasovskii's Method**: For x_dot = f(x), try V = f'Pf. Then V_dot = f'(J'P + PJ)f where J is the Jacobian. If J'P + PJ < 0 uniformly, the system is asymptotically stable.

2. **Variable Gradient Method**: Assume grad V = g(x) and enforce the curl-free condition (dg_i/dx_j = dg_j/dx_i). Then integrate g to obtain V(x), choosing g coefficients to make V_dot negative definite.

3. **Sum-of-Squares (SOS)**: For polynomial systems, formulate V(x) as an SOS polynomial. The positivity constraints reduce to semidefinite programming (SDP), solvable via LMIs.

### Converse Lyapunov Theorems

Massera (1949) and Kurzweil (1956) proved that asymptotic stability implies the existence of a smooth Lyapunov function. While non-constructive, these theorems justify the search -- if a system is stable, a Lyapunov function exists. Construction methods:
- Integral construction: V(x) = integral from 0 to infinity of |phi(t,x)|^2 dt
- Zubov's method: Solve V_dot = -phi(x)*(1-V) for region of attraction

### Implementation

| File | Purpose |
|------|---------|
| `lyap_core.h/c` | Core types, V and V_dot evaluation |
| `lyap_equation.h/c` | Bartels-Stewart Lyapunov equation solver |
| `lyap_construct.h/c` | Variable gradient and Krasovskii methods |
| `lyap_sos.h/c` | SOS polynomial Lyapunov functions |
| `lyap_linear.h/c` | Linear system analysis, quadratic stability |
| `lyap_converse.h/c` | Converse Lyapunov function construction |

### Key Functions

| Function | Description |
|----------|-------------|
| `lyap_solve_equation` | Solve A'P + PA = -Q via Bartels-Stewart |
| `lyap_evaluate` | Evaluate V(x) and V_dot(x) at a point |
| `lyap_krasovskii` | Krasovskii method: V = f(x)'Pf(x) |
| `lyap_variable_gradient` | Construct V via integrable gradient method |
| `lyap_check_positive_definite` | Verify V(x) > 0 for x != 0 |
| `lyap_estimate_region_attraction` | Level set {V <= c} as ROA estimate |

### Building

```
make          # Build static library
make test     # Build and run test suite
make examples # Build all example programs
```

### References

- Lyapunov, A.M. (1892). The General Problem of the Stability of Motion. Kharkov Mathematical Society.
- Khalil, H.K. (2002). Nonlinear Systems (3rd ed.). Prentice Hall.
- Slotine, J.J.E. & Li, W. (1991). Applied Nonlinear Control. Prentice Hall.
- Sastry, S. (1999). Nonlinear Systems: Analysis, Stability, and Control. Springer.
- Bacciotti, A. & Rosier, L. (2005). Liapunov Functions and Stability in Control Theory (2nd ed.). Springer.
- Boyd, S. et al. (1994). Linear Matrix Inequalities in System and Control Theory. SIAM.

### Examples

1. **example1**: Lyapunov equation solver for linear system
2. **example2**: Krasovskii method on nonlinear oscillator
3. **example3**: Variable gradient construction for 2D system

### Building

```
make && make test && make examples
```

### Quality Assurance

make test executes 15+ assert statements covering:
- Lyapunov equation solver (Bartels-Stewart algorithm)
- Positive definiteness verification
- Krasovskii and variable gradient construction
- Region of attraction estimation
- Numerical stability edge cases

### Requirements

GCC (C11), Make, ar, libm. No external dependencies.
