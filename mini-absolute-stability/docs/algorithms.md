# Algorithms — Absolute Stability

This document describes the key algorithms implemented in the
`mini-absolute-stability` module.

## 1. Lur'e System Representation

**Data Structure:** `AbsLureSystem`

The Lur'e system is defined by:
```
ẋ = A·x + b·u
y = cᵀ·x
u = -φ(t, y)
```
where φ belongs to sector [α, β].

**Creation:** `abs_lure_create(A, b, c, n, alpha, beta)`
- Allocates internal storage for A[n×n], b[n], c[n]
- Validates dimensions (1 ≤ n ≤ ABS_MAX_DIM)
- Swaps alpha/beta if alpha > beta
- Returns NULL on invalid input

**Validation:** `abs_lure_validate(sys)`
- Checks non-NULL pointers
- Checks finite floating-point values in A, b, c
- Verifies dimension constraints

## 2. Sector Operations

**Sector Make:** `abs_sector_make(alpha, beta)`
- Computes center = (α+β)/2, radius = (β-α)/2
- Ensures canonical ordering α ≤ β

**Sector Check:** `abs_sector_check(sec, y, phi_y)`
- Verifies: (φ-αy)·(βy-φ) ≥ 0

**Disk Center (Circle Criterion):**
```
center = -(α+β) / (2αβ)      for 0 < α < β
center = -1/β                 for α = 0
```

**Disk Radius (Circle Criterion):**
```
radius = (β-α) / (2|αβ|)      for α ≠ 0
radius = 1/(2β)               for α = 0
```

## 3. Circle Criterion

**Function:** `abs_circle_test(sys, n_freqs, w_min, w_max)`

Algorithm:
1. Compute Nyquist disk D(α,β) from sector bounds
2. Generate log-spaced frequency grid [w_min, w_max]
3. For each frequency ω:
   - Compute G(jω) = cᵀ·(jωI - A)⁻¹·b via complex linear solve
   - Compute distance from G(jω) to disk center
   - If distance ≤ radius: criterion FAILS
4. Report minimum margin and critical frequency

**Complex linear solve:** Form 2n×2n real augmented system:
```
[ Re(sI-A)  -Im(sI-A) ] [x_re]   [b]
[ Im(sI-A)   Re(sI-A) ] [x_im] = [0]
```
Solve via Gaussian elimination with partial pivoting, then G = cᵀ·x.

Complexity: O(n³) per frequency point, O(n_freqs · n³) total.

## 4. Popov Criterion

**Function:** `abs_popov_test(sys, n_freqs, w_min, w_max)`

Algorithm:
1. Map sector to [0, k] form (k = β - α)
2. Generate Popov plot data: X(ω) = Re[G(jω)], Y(ω) = ω·Im[G(jω)]
3. Find optimal Popov multiplier q* via golden-section search:
   - For each q, evaluate: min_ω { Re[G] - q·ω·Im[G] + 1/k }
   - Search q ∈ [0, 100] to maximize the minimum
4. Report pass/fail based on margin > 0

**Golden-section search:** Unimodal optimization with φ = (√5-1)/2 conjugate ratio, 40 iterations.

## 5. KYP Lemma / LMI Solver

**Function:** `abs_kyp_solve_strict(A, B, C, D, n, m, p)`

Algorithm (simplified for strictly proper square systems):
1. Form Q = Cᵀ·C + ε·I (regularization)
2. Solve Lyapunov equation Aᵀ·P + P·A = -Q via Bartels-Stewart
3. Check P > 0 via Cholesky factorization
4. Verify LMI[P] < 0 by eigenvalue check

**LMI Formation:**
```
LMI(P) = [ Aᵀ·P + P·A    P·B - Cᵀ ]
         [ Bᵀ·P - C      -(D+Dᵀ)  ]
```

## 6. Lyapunov Equation Solvers

### Bartels-Stewart Algorithm

The most numerically stable method (O(n³)):

1. Compute symmetric Schur decomposition: A ≈ U·Λ·Uᵀ
2. Transform Q̃ = Uᵀ·Q·U
3. Solve diagonal system: Λ·Y + Y·Λ = -Q̃
   For each (i,j): Yᵢⱼ = -Q̃ᵢⱼ / (λᵢ + λⱼ)
4. Recover P = U·Y·Uᵀ

### Kronecker Product Solver

Forms the n²×n² system:
```
(I ⊗ Aᵀ + Aᵀ ⊗ I) · vec(P) = -vec(Q)
```
Solves via Gaussian elimination. Only practical for n ≤ 4.

### Matrix Sign Function

Iterative method:
```
H₀ = [Aᵀ, Q; 0, -A]
H_{k+1} = (3·H_k - H_k³) / 2   (Newton-Schulz)
P = sign(H)_(1,2) block
```

## 7. Aizerman Conjecture Testing

**Function:** `abs_aizerman_test(sys)`

Algorithm:
1. Check trivial cases (n=1,2 → holds; A symmetric → holds)
2. Compute Hurwitz sector [k_min, k_max] via bisection
3. Check if [α, β] ⊆ [k_min, k_max]
4. Test for known counterexample patterns (Pliss n=3, Fitts n=4)
5. Compare with circle/Popov criteria conservatism

## 8. Eigenvalue Computation (Jacobi Method)

**Function:** `abs_linalg_symeig(n, A, eigenvalues)`

Classical Jacobi method for real symmetric matrices:
1. Find maximum off-diagonal element |A[p][q]|
2. Compute rotation angle: θ = ½·atan2(2A[p][q], A[p][p]-A[q][q])
3. Apply Jacobi rotation: A ← Jᵀ·A·J, V ← V·J
4. Repeat until max off-diagonal < tolerance
5. Diagonal entries converge to eigenvalues

Convergence: Quadratic in practice, guaranteed for symmetric matrices.

## Complexity Summary

| Algorithm        | Complexity           | Notes                    |
|------------------|---------------------|--------------------------|
| Circle criterion | O(n_freqs · n³)     | Dominated by complex solve |
| Popov criterion  | O(iter · n_freqs · n³) | Golden-section × freq sweep |
| Bartels-Stewart  | O(n³)               | One-time cost            |
| Kronecker        | O(n⁶)               | Small n only             |
| Matrix sign      | O(iter · n³)        | Newton-Schulz iteration  |
| Jacobi eig       | O(n³ · log(1/ε))    | Quadratic convergence    |
| Cholesky         | O(n³/3)             | Standard algorithm       |
| KYP solve        | O(n³)               | Dominated by Lyapunov    |
