# Stability of Interconnected Systems

## Core Concept
Stability analysis of interconnected subsystems using small-gain theorem,
passivity theory, composite Lyapunov functions, M-matrix conditions, and
decentralized control.

## Key Theorems
- **Small-Gain** (Zames 1966): ||G1|| * ||G2|| < 1 => feedback stable
- **Passivity** (Willems 1972): Passive + Strictly Passive => stable
- **M-Matrix** (Siljak 1978): M-matrix comparison => composite stability
- **Gershgorin**: Diagonal dominance => eigenvalue localization
- **Decentralized**: Local controllers + interaction bound check

## Build
  make          # Build library
  make test     # Run 60+ asserts
  make examples # Build all demos
  make demo     # Run all examples
  make clean    # Clean

## Quality
| Metric | Value |
|--------|-------|
| Headers | 5 |
| C sources | 5 |
| Lean files | 1 |
| Total lines | 3000+ |
| Functions | 60+ |
| Asserts | 60+ |
| Examples | 3 |
| Docs | 2 |
| Lean theorems | 14 |

## References
- Zames (1966) IEEE TAC
- Willems (1972) Arch. Rational Mech.
- Siljak (1978) Large-Scale Dynamic Systems
- Vidyasagar, M. (1980). Input-Output Analysis of Large-Scale Interconnected Systems. Springer.

## Theory Details

### Small-Gain Theorem for Interconnections

The interconnected small-gain theorem (Zames, 1966; Hill & Moylan, 1976): For N interconnected subsystems Hi, each with finite L2 gain gamma_i, the overall interconnection is L2 stable if the spectral radius rho(Gamma) < 1, where Gamma is the gain matrix with entries gamma_{ij} = gain from input of Hj to output of Hi.

### M-Matrix Conditions (Siljak, 1978)

An NxN matrix W is an M-matrix if w_{ij} <= 0 for i != j and all leading principal minors are positive. For interconnected systems with Lyapunov functions Vi and interaction bounds:

```
dV_i/dt <= -w_{ii}*|xi|^2 + sum_{j!=i} |w_{ij}|*|xi|*|xj|
```

If W is an M-matrix, the interconnection is globally exponentially stable. This provides a scalable stability test: verify each subsystem independently, then check the interconnection matrix.

### Composite Lyapunov Functions

Three main approaches to construct Lyapunov functions for interconnections:

1. **Weighted sum**: V = sum_i d_i * V_i(xi) with positive weights d_i chosen to dominate interconnection terms. The weights exist if the test matrix is an M-matrix.

2. **Vector Lyapunov functions**: V = [V1,...,VN]' with comparison principle. Each Vi must decrease unless dominated by other Vj's. The comparison system must be asymptotically stable.

3. **ISS-based**: V_i is an ISS-Lyapunov function for subsystem i. The small-gain condition on the gain matrix ensures overall stability.

### Decentralized Control

Each subsystem is controlled using only local information. Stability conditions involve:
- Local controllers that stabilize isolated subsystems
- Bounds on interconnection strengths (matching conditions or M-matrix)
- Robustness margins against neglected interactions

### Gershgorin Diagonal Dominance

A matrix A is row-diagonally dominant if |a_{ii}| > sum_{j!=i} |a_{ij}| for all i. For interconnected systems, diagonal dominance of the loop transfer matrix (after scaling) guarantees stability via the Nyquist array method (Rosenbrock, 1974).

### Connective Stability (Siljak)

The system remains stable even when subsystems are disconnected and reconnected arbitrarily -- a property crucial for power systems and communication networks where links may fail. If the disconnected subsystems are stable and the interconnection satisfies an M-matrix condition, the system is connectively stable.

### Implementation

| File | Purpose |
|------|---------|
| `intercon_core.h/c` | Core interconnection types, interaction graph |
| `small_gain_intercon.h/c` | Multi-system small-gain, spectral radius test |
| `m_matrix.h/c` | M-matrix verification, weight computation |
| `composite_lyapunov.h/c` | Weighted sum and vector Lyapunov functions |
| `decentralized.h/c` | Decentralized controller design and verification |

### Key Functions

| Function | Description |
|----------|-------------|
| `ic_small_gain_test` | Spectral radius test for interconnected gains |
| `ic_m_matrix_check` | Verify M-matrix condition for interconnection |
| `ic_composite_lyap_weights` | Compute optimal weighting vector d |
| `ic_decentralized_verify` | Verify decentralized stability |
| `ic_connective_stability` | Check connective stability property |
| `ic_gershgorin_check` | Gershgorin diagonal dominance test |
| `ic_gain_matrix_estimate` | Estimate pairwise interaction gains |

### Building

```
make && make test && make examples
```

### References
- Zames, G. (1966). On the input-output stability of time-varying nonlinear feedback systems. IEEE TAC, 11:228-238, 465-476.
- Willems, J.C. (1972). Dissipative dynamical systems. Archive for Rational Mechanics and Analysis, 45:321-393.
- Siljak, D.D. (1978). Large-Scale Dynamic Systems: Stability and Structure. North-Holland.
- Vidyasagar, M. (1980). Input-Output Analysis of Large-Scale Interconnected Systems. Springer.
- Michel, A.N. & Miller, R.K. (1977). Qualitative Analysis of Large Scale Dynamical Systems. Academic Press.
- Lunze, J. (1992). Feedback Control of Large-Scale Systems. Prentice Hall.
