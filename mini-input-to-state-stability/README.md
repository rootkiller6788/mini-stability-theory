# Input-to-State Stability (ISS)

## Bridging Lyapunov and input-output stability

Input-to-State Stability, introduced by Eduardo Sontag in 1989, provides a unified framework combining Lyapunov stability with robustness to inputs. ISS characterizes how the state responds to external disturbances.

### ISS Definition

A system x_dot = f(x, u) is Input-to-State Stable if there exist a KL-function beta and a K-function gamma such that: |x(t)| <= beta(|x0|, t) + gamma(||u||_infinity). The state is bounded by a decaying function of the initial condition plus a function of the maximum input magnitude.

### ISS-Lyapunov Function

A smooth function V is an ISS-Lyapunov function if there exist alpha1, alpha2, alpha3 in K-infinity and sigma in K such that: alpha1(|x|) <= V(x) <= alpha2(|x|) and grad_V * f(x,u) <= -alpha3(|x|) + sigma(|u|). Existence of an ISS-Lyapunov function is equivalent to ISS (Sontag & Wang, 1995).

### Key Properties

1. ISS implies 0-GAS (when u=0, the system is globally asymptotically stable)
2. ISS implies BIBS (bounded inputs produce bounded states)
3. ISS superposition: sum of ISS systems remains ISS
4. ISS cascade: cascade connection of ISS systems preserves ISS
5. ISS is preserved under bicontinuous coordinate changes

### ISS Small-Gain Theorem

For two ISS systems interconnected in feedback, the closed-loop is ISS if: gamma1(gamma2(s)) < s for all s > 0 (nonlinear small-gain condition). This generalizes the classical linear small-gain theorem to nonlinear systems with comparison functions.

### Implementation

| File | Purpose |
|------|---------|
| `include/iss_core.h` / `src/iss_core.c` | Core ISS definitions, K/KL comparison functions, state operations, RK4 integration, trajectory simulation |
| `include/iss_lyapunov.h` / `src/iss_lyapunov.c` | ISS-Lyapunov verification (grid, random, Monte Carlo), dissipation check, decay rate, quadratic forms |
| `include/iss_gain.h` / `src/iss_gain.c` | ISS gain computation (Lyapunov-based, simulation, H-inf, frequency response, LTI) |
| `include/iss_small_gain.h` / `src/iss_small_gain.c` | Small-gain theorem (linear, nonlinear, cyclic, spectral radius, Lyapunov-based) |
| `include/iss_interconnected.h` / `src/iss_interconnected.c` | Cascade and feedback ISS, parallel interconnections, trajectory verification |
| `include/iss_verify.h` / `src/iss_verify.c` | Full verification pipeline, adaptive sampling, robustness, batch verification |
| `src/iss_stability.lean` | Lean 4 formalization of ISS theorems |

### Key Functions

| Function | Description |
|----------|-------------|
| `iss_k_linear` / `iss_k_power` / `iss_k_saturation` | Create K-class comparison functions |
| `iss_kl_exponential` / `iss_kl_polynomial` | Create KL-class comparison functions |
| `iss_system_create` / `iss_system_free` | System lifecycle |
| `iss_lyapunov_verify_random` | Monte Carlo ISS-Lyapunov verification |
| `iss_lyapunov_vdot` | Compute Lie derivative of Lyapunov function |
| `iss_gain_from_lyapunov` | ISS gain from certificate |
| `iss_gain_from_simulation` | ISS gain via trajectory simulation |
| `iss_small_gain_condition` | Check linear small-gain condition |
| `iss_cyclic_small_gain` | Cyclic small-gain for N systems |
| `iss_nonlinear_small_gain` | Nonlinear small-gain via K-function composition |
| `iss_cascade_create` / `iss_cascade_is_ISS` | Cascade ISS analysis |
| `iss_feedback_create` / `iss_feedback_is_ISS` | Feedback ISS analysis |
| `iss_verify` | Full ISS verification pipeline |
| `iss_verify_kinf` | Verify K-infinity property numerically |
| `iss_verify_convergence` | Trajectory convergence check |
| `iss_verify_gain_empirical` | Empirical ISS gain estimation |
| `iss_bound_full` | Full ISS bound: beta(|x0|, t) + gamma(||u||) |
| `iss_spectral_radius` | Spectral radius for small-gain matrix |
| `iss_rk4_step` / `iss_rk4_integrate` | Runge-Kutta 4 integration |
| `iss_trajectory_simulate` | Full trajectory simulation |

### Building

```bash
# Build the static library
make

# Run all tests
make test

# Build and run examples
make examples

# Clean build artifacts
make clean
```

### Module Structure

```
mini-input-to-state-stability/
  Makefile
  README.md                          (this file)
  include/
    iss_core.h                       Core types and comparison functions
    iss_lyapunov.h                   ISS-Lyapunov verification
    iss_gain.h                       ISS gain computation
    iss_small_gain.h                 Small-gain theorem
    iss_interconnected.h             Cascade/feedback ISS
    iss_verify.h                     Full verification pipeline
  src/
    iss_core.c                       ~500 lines, K/KL functions, RK4, simulation
    iss_lyapunov.c                   ~340 lines, Lyapunov verification, dissipation
    iss_gain.c                       ~260 lines, Gain computation, frequency response
    iss_small_gain.c                 ~210 lines, Small-gain, spectral radius
    iss_interconnected.c             ~220 lines, Cascade, feedback, trajectory check
    iss_verify.c                     ~260 lines, Verification pipeline, robustness
    iss_stability.lean               Formal ISS proofs in Lean 4
  tests/
    test_iss.c                       21 test functions with assert()
  examples/
    example1_iss_siso.c              SISO ISS system demo
    example2_cascade.c               Cascade ISS system demo
    example3_small_gain.c            Small-gain theorem demo
  docs/
    iss-theory.md                    Comprehensive ISS theory reference
    course-alignment.md              Mapping to MIT/Stanford/ETH curricula
```

### ISS Verification Pipeline

1. **System Setup**: Define dynamics `f(x,u)` and Lyapunov candidate `V(x)`
2. **Dissipation Check**: Verify `Vdot(x,u) <= -alpha(|x|) + sigma(|u|)` via random/grid sampling
3. **Gain Estimation**: Compute ISS gain `gamma` from alpha and sigma
4. **Small-Gain (if interconnected)**: Check `gamma1(gamma2(s)) < s`
5. **Certificate**: Issue ISS certificate with gain and margin

### K-Class Function Library

This module provides a comprehensive library of comparison functions:

| Type | Functions | Description |
|------|-----------|-------------|
| **K** | linear, power, saturation, cubic, exponential, deadzone, sigmoid | gamma(0)=0, strictly increasing |
| **K-infinity** | All K + radial unboundedness | K-class with gamma(s)→∞ as s→∞ |
| **KL** | exponential, polynomial, double exponential | beta(s,t) decreasing in t |

Each K-function is **self-contained with heap-allocated parameters** — creating multiple K-functions does not overwrite shared state (fixing a common implementation bug).

### Dependencies

- C11 compiler (gcc/clang)
- Standard C library (stdlib, math, string, assert)
- Lean 4 (for formal verification only; not required for compilation)

### References

- Sontag, E.D. (1989). Smooth stabilization implies coprime factorization. IEEE Trans. Automatic Control, 34:435-443.
- Sontag, E.D. & Wang, Y. (1995). On characterizations of the ISS property. Systems & Control Letters, 24:351-359.
- Jiang, Z.P., Teel, A.R. & Praly, L. (1994). Small-gain theorem for ISS systems. Mathematics of Control, Signals and Systems, 7:95-120.
- Sontag, E.D. (2008). Input to state stability: Basic concepts and results. Lecture Notes in Mathematics, 1932:163-220.
- Karafyllis, I. & Jiang, Z.P. (2011). Stability and Stabilization of Nonlinear Systems. Springer.
- Khalil, H.K. (2002). Nonlinear Systems (3rd ed.). Prentice Hall.
