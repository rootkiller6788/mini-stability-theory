# mini-absolute-stability

Absolute stability theory for Lur'e systems with sector-bounded nonlinearities.

## Overview

This module provides a comprehensive C implementation of the classical
absolute stability theory, including the **Circle Criterion**, **Popov
Criterion**, **KYP Lemma**, and **Aizerman Conjecture** analysis.

A Lur'e system consists of a linear forward path with nonlinear feedback:
```
ẋ = A·x + b·u      (linear plant)
y = cᵀ·x            (output)
u = -φ(t, y)        (sector-bounded nonlinearity)
```

**Absolute stability** means the origin is globally asymptotically stable
for ALL nonlinearities φ in a given sector [α, β].

## Mathematical Background

### Sector Condition
A nonlinearity φ belongs to sector [α, β] if:
```
α·y² ≤ y·φ(y) ≤ β·y²    for all y
```

### Circle Criterion (SISO)
For 0 < α < β, the Nyquist plot of G(jω) must not intersect the disk
D(α,β) with center = -(α+β)/(2αβ) and radius = (β-α)/(2αβ).

### Popov Criterion
For time-invariant φ ∈ [0, k]: ∃ q ≥ 0 such that
```
Re[(1 + jωq)·G(jω)] + 1/k > 0    ∀ω ≥ 0
```

### KYP Lemma
Bridges frequency-domain inequalities and linear matrix inequalities:
```
H(jω) + H(jω)ᴴ ≥ 0  ⇔  ∃P = Pᵀ > 0: LMI(P) ≤ 0
```

### Aizerman Conjecture (FALSE for n ≥ 3)
Absolute stability in [α, β] is NOT equivalent to Hurwitz stability
of A - k·b·cᵀ for all k ∈ [α, β] (counterexamples exist).

## Module Structure

```
mini-absolute-stability/
├── Makefile              # Build system (gcc -std=c11)
├── README.md             # This file
├── include/
│   ├── abs_core.h        # Core types: LureSystem, Sector, Complex, LA
│   ├── abs_lyapunov.h    # Lyapunov equation solvers (Bartels-Stewart, etc.)
│   ├── abs_circle.h      # Circle criterion: Nyquist disk test
│   ├── abs_popov.h       # Popov criterion: multiplier method
│   ├── abs_kyp.h         # KYP lemma: LMI solver, SPR/PR checks
│   └── abs_aizerman.h    # Aizerman conjecture: Hurwitz sector, counterexamples
├── src/
│   ├── abs_core.c        # LureSystem lifecycle, sector ops, linear algebra
│   ├── abs_lyapunov.c    # Bartels-Stewart, Kronecker, matrix sign solvers
│   ├── abs_circle.c      # Circle criterion implementation
│   ├── abs_popov.c       # Popov criterion implementation
│   ├── abs_kyp.c         # KYP lemma and LMI solver
│   ├── abs_aizerman.c    # Aizerman testing, counterexample construction
│   └── abs_theory.lean   # Lean 4 formalization (theorems & structures)
├── tests/
│   └── test_abs.c        # Comprehensive tests (>=25 assert() calls)
├── examples/
│   ├── example1.c        # Circle criterion demo (2nd-order system)
│   ├── example2.c        # Popov criterion demo (3rd-order, comparison)
│   └── example3.c        # Aizerman + KYP demo (criteria comparison)
└── docs/
    ├── theory.md         # Theoretical background and derivations
    └── algorithms.md     # Algorithm descriptions and complexity analysis
```

## Building

```bash
# Build static library
make

# Run tests
make test

# Clean build artifacts
make clean
```

### Requirements
- GCC (C11 standard)
- Make

## API Quick Reference

### Core (abs_core.h)
| Function | Description |
|----------|-------------|
| `abs_lure_create()` | Create a Lur'e system |
| `abs_lure_free()` | Free a Lur'e system |
| `abs_sector_make()` | Create sector [alpha,beta] |
| `abs_linalg_symeig()` | Symmetric eigenvalue decomposition (Jacobi) |
| `abs_linalg_cholesky()` | Cholesky factorization |
| `abs_linalg_freqresp()` | Frequency response G(jw) |

### Circle Criterion (abs_circle.h)
| Function | Description |
|----------|-------------|
| `abs_circle_test()` | Full circle criterion test |
| `abs_circle_quick_check()` | Fast pass/fail check |
| `abs_circle_compute_disk()` | Compute Nyquist disk |

### Popov Criterion (abs_popov.h)
| Function | Description |
|----------|-------------|
| `abs_popov_test()` | Full Popov criterion test |
| `abs_popov_find_optimal_q()` | Find optimal Popov multiplier |
| `abs_popov_max_stable_k()` | Maximum stable sector bound |

### Lyapunov (abs_lyapunov.h)
| Function | Description |
|----------|-------------|
| `abs_lyap_bartels_stewart()` | Bartels-Stewart algorithm (O(n^3)) |
| `abs_lyap_kronecker()` | Kronecker product solver (small n) |
| `abs_lyap_matrix_sign()` | Matrix sign iteration |
| `abs_lyap_for_lure()` | Lyapunov function for Lur'e system |
| `abs_lyap_separation()` | Sep(A^T, -A) estimation |

### KYP Lemma (abs_kyp.h)
| Function | Description |
|----------|-------------|
| `abs_kyp_solve_strict()` | Solve strict KYP (SPR) |
| `abs_kyp_circle_to_lmi()` | KYP from circle criterion |
| `abs_kyp_is_spr()` | Check strict positive realness |
| `abs_kyp_is_negative_semidef()` | Check negative semidefiniteness |

### Aizerman (abs_aizerman.h)
| Function | Description |
|----------|-------------|
| `abs_aizerman_test()` | Test Aizerman conjecture |
| `abs_aizerman_hurwitz_sector()` | Compute Hurwitz sector |
| `abs_aizerman_pliss_counterexample()` | Construct Pliss counterexample |
| `abs_aizerman_conservatism_gap()` | Compare Aizerman vs criteria |

## Numerical Methods

| Method | Algorithm | Complexity |
|--------|-----------|------------|
| G(jw) evaluation | Complex linear solve (2n augmented) | O(n^3) |
| Lyapunov equation | Bartels-Stewart (Schur + Sylvester) | O(n^3) |
| Lyapunov equation | Kronecker product | O(n^6) |
| Eigenvalues | Jacobi method (symmetric) | O(n^3) |
| Optimal q | Golden-section search | O(iter * n_freqs * n^3) |
| Cholesky | Standard algorithm | O(n^3/3) |
| KYP LMI | Lyapunov + eigenvalue check | O(n^3) |

## Constants

| Constant | Value | Description |
|----------|-------|-------------|
| `ABS_PI` | 3.14159... | pi |
| `ABS_EPS` | 1e-12 | Default tolerance |
| `ABS_MAX_DIM` | 16 | Maximum system dimension |
| `ABS_MAX_ITER` | 1000 | Maximum iterations |

## Structs

- **AbsLureSystem** — Lur'e problem (A, b, c, alpha, beta, n)
- **AbsSector** — Sector [alpha, beta] with center and radius
- **AbsComplex** — Complex number (re, im)
- **AbsCircleDisk** — Nyquist disk for circle criterion
- **AbsPopovResult** — Popov criterion result
- **AbsKYPResult** — KYP lemma solution
- **AbsAizermanResult** — Aizerman conjecture test result
- **AbsLyapResult** — Lyapunov equation solution
- **AbsFrechet** — Frechet derivative structure

## Exported Functions

Over 40 exported functions across 6 modules, covering:
- Lur'e system creation and management
- Sector condition verification
- Three Lyapunov equation solvers
- Circle criterion with Nyquist sweep
- Popov criterion with optimal multiplier search
- KYP lemma and LMI solver
- Aizerman conjecture testing
- Full linear algebra toolbox (matrix ops, eigenvalues, Cholesky)

## Lean 4 Formalization

The file `src/abs_theory.lean` provides formal structures and theorem
statements (with partial proofs) for:
- `SectorCondition` — formal sector definition
- `LureSystem` — state-space representation
- `LyapunovFunction` — quadratic Lyapunov functions
- `kyp_lemma_strict` / `kyp_lemma_nonstrict` — KYP equivalence
- `circle_criterion_sufficiency` — circle criterion theorem
- `popov_criterion_sufficiency` — Popov criterion theorem
- `aizerman_conjecture_false` / `aizerman_holds_second_order`
- `lure_lyapunov_valid` — Lur'e-Postnikov Lyapunov function

## References

- Khalil, H.K. (2002). *Nonlinear Systems*, 3rd ed. Prentice Hall.
- Vidyasagar, M. (1993). *Nonlinear Systems Analysis*, 2nd ed. SIAM.
- Haddad, W.M. & Chellaboina, V. (2008). *Nonlinear Dynamical Systems and Control*. Princeton.
- Lur'e, A.I. & Postnikov, V.N. (1944). On the theory of stability of control systems.
- Aizerman, M.A. (1949). On a problem concerning stability in the large of dynamical systems.
- Popov, V.M. (1961). Absolute stability of nonlinear systems of automatic control.
- Kalman, R.E. (1963). Lyapunov functions for the problem of Lur'e in automatic control.

## License

Educational/academic use. Part of the mini-complex-control-theory project.
