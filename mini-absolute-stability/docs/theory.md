# Absolute Stability Theory

## Historical Background

The absolute stability problem was formulated by A.I. Lur'e and V.N.
Postnikov in 1944, studying the stability of control systems with a
single nonlinear element. The problem gained prominence through the
work of M.A. Aizerman (1949), who conjectured that stability of the
nonlinear system could be deduced from stability of all linear systems
in a given sector. While the conjecture is false in general, it
spurred the development of rigorous frequency-domain criteria.

## Lur'e System Definition

A Lur'e system consists of a linear time-invariant forward path with
nonlinear feedback:

```
ẋ = A·x + b·u          (state equation)
y = cᵀ·x               (output equation)
u = -φ(t, y)           (nonlinear feedback)
```

where:
- **x** ∈ ℝⁿ: state vector
- **u** ∈ ℝ: scalar control input
- **y** ∈ ℝ: scalar output
- **A** ∈ ℝⁿˣⁿ: Hurwitz state matrix
- **b, c** ∈ ℝⁿ: input/output coupling vectors
- **φ**: ℝ×ℝ → ℝ: nonlinearity belonging to sector [α, β]

The transfer function of the linear part is:
```
G(s) = cᵀ·(sI - A)⁻¹·b
```

## Sector Condition

A function φ(y) belongs to sector [α, β] (denoted φ ∈ [α, β]) if:

```
α·y² ≤ y·φ(y) ≤ β·y²    for all y ∈ ℝ
```

Equivalently:
```
(φ(y) - α·y)·(β·y - φ(y)) ≥ 0    for all y
```

Common nonlinearity types:
- **Saturation**: φ(y) = sat(y/L) · L ∈ [0, 1]
- **Deadzone**: φ(y) = max(0, |y|-δ)·sgn(y) ∈ [0, 1]
- **Relay**: φ(y) = sgn(y) ∈ [0, ∞)
- **Cubic**: φ(y) = y³ ∈ [0, ∞)
- **Sinusoidal**: φ(y) = sin(y) ∈ [sinc(π), 1] on bounded intervals

## Absolute Stability

**Definition:** The Lur'e system is absolutely stable on sector [α, β] if
the origin x = 0 is globally uniformly asymptotically stable for all
nonlinearities φ ∈ [α, β].

Key point: Absolute stability guarantees stability for the *entire class*
of sector-bounded nonlinearities, not just a specific φ.

## Circle Criterion

### Theorem (SISO Circle Criterion)

For the Lur'e system with φ ∈ [α, β] where 0 < α < β:

The system is absolutely stable if the Nyquist plot of G(jω) does not
intersect or encircle the disk D(α,β) with:
```
center: -(α+β)/(2αβ)    (on real axis)
radius: (β-α)/(2αβ)
```

For α = 0, the condition reduces to:
```
Re[G(jω)] > -1/β    for all ω ≥ 0
```

### Intuition

The circle criterion can be understood as a graphical test: draw the
Nyquist plot of G(jω) and see if it avoids the critical disk. If it
does, every linearization G̃(s) = G(s)/(1+kG(s)) for k∈[α,β] is stable.

### Proof Sketch

Via the KYP lemma and loop transformation, the circle condition is
equivalent to the existence of P = Pᵀ > 0 satisfying the LMI.

## Popov Criterion

### Theorem (Popov Criterion)

For the Lur'e system with *time-invariant* φ ∈ [0, k]:

The system is absolutely stable if there exists q ≥ 0 such that:
```
Re[(1 + jωq)·G(jω)] + 1/k > 0    for all ω ≥ 0
```

### Geometrical Interpretation

Define the Popov plot:
```
X(ω) = Re[G(jω)],   Y(ω) = ω·Im[G(jω)]
```

The condition requires that all points (X(ω), Y(ω)) lie strictly to
the right of a line (the "Popov line") passing through (-1/k, 0) with
slope 1/q.

### Lyapunov Function

The Popov criterion is proved using the Lur'e-Postnikov Lyapunov function:
```
V(x) = xᵀ·P·x + q·∫₀^{cᵀx} φ(σ) dσ
```

The integral term captures the energy stored in the nonlinearity.

### Comparison with Circle Criterion

| Property             | Circle Criterion     | Popov Criterion     |
|----------------------|----------------------|----------------------|
| Nonlinearity type    | Time-varying         | Time-invariant       |
| Conservatism         | More conservative    | Less conservative    |
| Graphical test       | Nyquist + disk       | Popov line           |
| Multiplier           | Constant             | M(s) = 1 + q·s       |
| Lyapunov function    | Quadratic only       | Quadratic + integral |

The Popov criterion is strictly sharper for time-invariant φ, but both
are sufficient (not necessary) conditions.

## KYP Lemma

### Statement

Given a minimal realization (A, B, C, D) with A Hurwitz:

The frequency-domain inequality:
```
H(jω) + H(jω)ᴴ ≥ 0    for all ω ∈ ℝ
```
where H(s) = C·(sI-A)⁻¹·B + D

is equivalent to the existence of P = Pᵀ ≥ 0 such that:
```
[ Aᵀ·P + P·A    P·B - Cᵀ ]  ≤ 0
[ Bᵀ·P - C      -(D+Dᵀ)  ]
```

### Significance

The KYP lemma is the fundamental bridge between:
- **Frequency-domain** methods (Nyquist, Bode, circle/Popov criteria)
- **State-space** methods (Lyapunov functions, LMIs)

It enables converting graphical frequency-domain tests into algebraic
matrix inequalities, which can be solved numerically.

### Applications

1. **Circle criterion → LMI**: Transform sector condition to SPR condition
2. **Popov criterion → LMI**: Augment system with multiplier dynamics
3. **Bounded-real lemma**: H∞ norm condition → Riccati inequality
4. **Passivity analysis**: Energy-based stability proofs

## Aizerman Conjecture

### Statement (1949)

The Lur'e system is absolutely stable in sector [α, β] if and only if
the linear system ẋ = (A - k·b·cᵀ)·x is stable for all k ∈ [α, β].

### Status: FALSE

Counterexamples exist for n ≥ 3:
- **Pliss (1958)**: 3rd-order counterexample
- **Fitts (1966)**: 4th-order counterexample
- **Barabanov (1988)**: Systematic construction method

### Where It Holds

The conjecture is true for:
1. **First and second order systems** (n = 1, 2)
2. **Systems with symmetric A matrix** (circle criterion proves it)
3. **Positive systems** (Metzler matrices)
4. **Minimum-phase systems** with relative degree ≤ 1 (Popov applies)

### Significance

While false, the Aizerman conjecture was historically important as it:
- Motivated the search for rigorous frequency-domain criteria
- Clarified the gap between linear and nonlinear stability
- Led to the development of multiplier methods

## References

1. **Lur'e & Postnikov** (1944) — Original formulation
2. **Aizerman** (1949) — The conjecture
3. **Kalman** (1963) — KYP lemma and its implications
4. **Popov** (1961) — Frequency-domain criterion for time-invariant systems
5. **Zames** (1966) — Circle criterion and input-output stability
6. **Khalil** (2002) — "Nonlinear Systems", Chapter 7
7. **Vidyasagar** (1993) — "Nonlinear Systems Analysis", Chapter 6
8. **Haddad & Chellaboina** (2008) — LMI-based approach

## Related Modules

- `mini-circle-popov-criterion` — Detailed frequency-domain analysis
- `mini-lyapunov-direct-method` — Lyapunov function construction
- `mini-input-output-stability` — I/O stability theory
- `mini-finite-time-stability` — Finite-time convergence
