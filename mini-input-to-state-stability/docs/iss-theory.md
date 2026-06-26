# Input-to-State Stability (ISS) — Theory Reference

## 1. Definition and Fundamentals

Input-to-State Stability (ISS), introduced by Eduardo Sontag in 1989, bridges Lyapunov stability theory with input-output robustness analysis.

### Formal Definition

A system `dx/dt = f(x, u)` is **Input-to-State Stable (ISS)** if there exist:
- A KL-class function `beta` (decaying bound on initial condition)
- A K-class function `gamma` (asymptotic input gain)

such that for all initial conditions x(0) = x0 and all essentially bounded inputs u(·):

```
|x(t)| <= beta(|x0|, t) + gamma(||u||_inf)
```

where `|·|` is the Euclidean norm and `||u||_inf = sup_{t>=0} |u(t)|`.

### Interpretation

- `beta(|x0|, t)` bounds the transient from the initial condition
- `gamma(||u||_inf)` bounds the asymptotic effect of the input
- When u ≡ 0, ISS reduces to global asymptotic stability (0-GAS)
- Bounded inputs produce bounded states (BIBS property)

## 2. ISS-Lyapunov Functions

A smooth function `V: R^n -> R_>=0` is an **ISS-Lyapunov function** for `dx/dt = f(x,u)` if there exist K-infinity functions alpha1, alpha2, alpha3 and a K-function sigma such that:

```
alpha1(|x|) <= V(x) <= alpha2(|x|)        (positive definite, radially unbounded)
nabla V · f(x,u) <= -alpha3(|x|) + sigma(|u|)  (dissipation inequality)
```

### Sontag-Wang Theorem (1995)

**Existence of an ISS-Lyapunov function is necessary and sufficient for ISS.**

This equivalence is the cornerstone of ISS analysis:
- **Sufficiency**: If V satisfies the dissipation inequality, the system is ISS with gain `gamma = alpha1^{-1} o alpha2 o alpha3^{-1} o sigma`
- **Necessity**: If the system is ISS, a smooth ISS-Lyapunov function can be constructed (converse Lyapunov theorem)

## 3. ISS Gain Computation

The ISS gain `gamma` determines how the input magnitude maps to the state bound.

### From Lyapunov Function

Given `Vdot <= -alpha3(|x|) + sigma(|u|)`:
```
gamma(s) = alpha1^{-1} o alpha2 o alpha3^{-1} o sigma(s)
```

### For Linear Systems

For `dx/dt = Ax + Bu` with Hurwitz A:
```
gamma(s) = ||(sI - A)^{-1} B||_H-infinity * s
```
The H-infinity norm is the peak of the frequency response.

### Special Cases
- **Linear ISS**: `gamma(s) = K*s` with gain K = `||int_0^inf e^{At} B dt||`
- **Quadratic Lyapunov**: `V(x) = x'Px`, gain determined by the Lyapunov equation
- **Monotone systems**: Gain composition via max-type functions

## 4. Small-Gain Theorem

### Two-System Case

For interconnected ISS systems:
```
dx1/dt = f1(x1, x2, u1)    (ISS with gain gamma1 w.r.t. x2)
dx2/dt = f2(x2, x1, u2)    (ISS with gain gamma2 w.r.t. x1)
```

**Nonlinear Small-Gain Condition**: If `gamma1(gamma2(s)) < s` for all s > 0, the interconnection is ISS.

### Linear Small-Gain

For linear gains `gamma_i(s) = K_i * s`:
```
K1 * K2 < 1
```

### Cyclic Small-Gain (N systems)

For N ISS systems connected in a cycle with gains gamma_{i,i+1}:
The interconnection is ISS if the composition `gamma1 o gamma2 o ... o gammaN(s) < s` for all s > 0.

## 5. Integral ISS (iISS)

A relaxed notion where the dissipation involves the **integral** of the input:

```
alpha1(|x|) <= V(x) <= alpha2(|x|)
nabla V · f(x,u) <= -alpha3(|x|) + gamma(|u|)
```
where alpha3 is only positive definite (not necessarily K-infinity).

iISS is useful for systems where large but brief inputs have limited effect.

## 6. ISS of Interconnections

### Cascade

If S1 is ISS w.r.t. x2 and S2 is 0-GAS, then the cascade (S1, S2) is 0-GAS.

### Parallel

The parallel connection of ISS systems is ISS with gain `max(gamma1, gamma2)`.

### Feedback

Requires the small-gain condition. The composite ISS-Lyapunov function is:
```
V(x1, x2) = V1(x1) + k*V2(x2)
```
for some `k > 0` chosen to satisfy the dissipation inequality.

## 7. ISS for Time-Delay Systems

For systems with delays `dx/dt = f(x(t), x(t-tau), u(t))`:
- Use Lyapunov-Krasovskii functionals
- ISS property provides robustness to both input disturbances and delay mismatch
- Razumikhin-type theorems characterize ISS for delay systems

## 8. Key Theorems (Summary)

| Theorem | Statement |
|---------|-----------|
| **Sontag-Wang** | ISS ⇔ ∃ ISS-Lyapunov function |
| **Small-Gain** | γ1∘γ2(s) < s ⇒ interconnection ISS |
| **Cascade ISS** | ISS+0-GAS cascade ⇒ 0-GAS |
| **Superposition** | Sum of ISS systems is ISS |
| **ISS ⇒ BIBS** | Bounded input ⇒ bounded state |
| **ISS ⇒ 0-GAS** | Zero input ⇒ asymptotic stability |

## References

- Sontag, E.D. (1989). "Smooth stabilization implies coprime factorization." IEEE TAC, 34:435-443.
- Sontag, E.D. & Wang, Y. (1995). "On characterizations of the ISS property." Systems & Control Letters, 24:351-359.
- Jiang, Z.P., Teel, A.R. & Praly, L. (1994). "Small-gain theorem for ISS systems." MCSS, 7:95-120.
- Sontag, E.D. (2008). "Input to state stability: Basic concepts and results." LNM, 1932:163-220.
- Khalil, H.K. (2002). "Nonlinear Systems" (3rd ed.). Prentice Hall.
