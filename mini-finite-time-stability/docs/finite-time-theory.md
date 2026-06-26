# Finite-Time Stability Theory

## Definition

A dynamical system `x_dot = f(x)` with equilibrium at `x=0` is **finite-time stable** if:

1. It is Lyapunov stable (for any epsilon>0, there exists delta>0 such that ||x0||<delta implies ||x(t)||<epsilon for all t>=0)
2. There exists a settling-time function `T(x0) < infinity` such that `lim_{t->T(x0)} x(t) = 0` and `x(t) = 0` for all `t >= T(x0)`.

## Lyapunov Condition (Bhat & Bernstein, 2000)

**Theorem**: If there exists a C^1 function V(x) >= 0 and constants c>0, 0<alpha<1 such that:

```
dV/dt <= -c * V(x)^alpha
```

Then the origin is finite-time stable with settling time bound:

```
T <= V(x0)^{1-alpha} / (c * (1-alpha))
```

### Proof Sketch

From `dV/dt <= -c*V^alpha`:

```
dV / V^alpha <= -c * dt
int_{V(0)}^{V(T)} V^{-alpha} dV <= -c * T
V(T)^{1-alpha}/(1-alpha) - V(0)^{1-alpha}/(1-alpha) <= -c * T
-V(0)^{1-alpha}/(1-alpha) <= -c * T  (since V(T)=0)
T = V(0)^{1-alpha} / (c*(1-alpha))
```

## Homogeneity Approach

**Theorem (Bhat & Bernstein, 2005)**: If the origin is asymptotically stable and the vector field f is homogeneous of negative degree d<0, then the origin is finite-time stable.

A function `f: R^n -> R^n` is homogeneous of degree d with respect to weights `(r1,...,rn)` if for all lambda>0:

```
f_i(lambda^{r1} x1, ..., lambda^{rn} xn) = lambda^{d+ri} f_i(x)
```

## Comparison of Stability Types

| Type | Convergence | Settling Time |
|------|-------------|---------------|
| Lyapunov | ||x(t)|| < epsilon for t>=0 | N/A |
| Asymptotic | lim_{t->inf} x(t) = 0 | Unbounded |
| Exponential | ~exp(-lambda*t) | Depends on x0 |
| Finite-time | x(t)=0 for t>=T(x0) | T <= V0^{1-alpha}/(c*(1-alpha)) |
| Fixed-time | x(t)=0 for t>=T_max | Uniform bound for all x0 |
| Prescribed-time | x(t)=0 for t>=T_p | Designer-specified T_p |

## Key Papers

1. Bhat, S.P. & Bernstein, D.S. (2000). "Finite-time stability of continuous autonomous systems." SIAM J. Control Optim., 38(3):751-766.
2. Bhat, S.P. & Bernstein, D.S. (2005). "Geometric homogeneity with applications to finite-time stability." Math. Control Signals Syst., 17:101-127.
3. Moulay, E. & Perruquetti, W. (2006). "Finite time stability and stabilization: State of the art." LNCIS, 334:23-41.
