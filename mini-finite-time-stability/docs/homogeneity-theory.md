# Homogeneity Theory for Finite-Time Stability

## Homogeneous Vector Fields

A vector field `f: R^n -> R^n` is homogeneous of degree `d in R` with respect to the dilation weights `r = (r1,...,rn)`, `ri > 0`, if for all `lambda > 0`:

```
f_i(lambda^{r1} x1, ..., lambda^{rn} xn) = lambda^{d+ri} f_i(x)
```

The dilation operator is defined as:

```
Delta_lambda(x) = (lambda^{r1} x1, ..., lambda^{rn} xn)
```

## Dilation Group Properties

Dilations form a one-parameter transformation group:

1. **Identity**: `Delta_1(x) = x`
2. **Composition**: `Delta_a(Delta_b(x)) = Delta_{ab}(x)`
3. **Inverse**: `Delta_{1/lambda} = Delta_lambda^{-1}`

## Homogeneous Norm

A homogeneous norm `||.||_r` satisfies:

```
||Delta_lambda(x)||_r = lambda * ||x||_r
```

The canonical construction is:

```
||x||_{r,p} = (sum_i |xi|^{p/ri})^{1/p}
```

where `p > 0`, typically `p = max_i ri`.

## FTS via Homogeneity

**Key Theorem**: If `x=0` is an asymptotically stable equilibrium of `x_dot = f(x)` and `f` is homogeneous of degree `d < 0`, then `x=0` is finite-time stable.

**Intuition**: Negative homogeneity degree means the vector field grows faster than linearly near the origin, which forces trajectories to reach the origin in finite time.

## Standard Weights

For a chain of integrators of length n:
```
x1_dot = x2
x2_dot = x3
...
xn_dot = u
```

The standard weights are `r = (n, n-1, ..., 1)`, giving homogeneity degree `d = -1`.

## Bi-Limit Homogeneity

If `f(x) = f_d(x) + f_{>d}(x)` where `f_d` is the d-homogeneous approximation (lower-degree terms vanish faster under dilation), then:

- Asymptotic stability of `f_d` with `d<0` implies FTS of the full system `f`.
- This allows analyzing complex systems via their homogeneous approximations.

## Key Properties

1. **Scaling of trajectories**: If `x(t)` is a trajectory from `x0`, then `Delta_lambda(x(lambda^d t))` is a trajectory from `Delta_lambda(x0)`.
2. **Settling time scaling**: `T(Delta_lambda(x0)) = lambda^{-d} * T(x0)`
3. **For d<0**: Smaller initial conditions take LONGER to converge (counter-intuitive!)

## Examples

### Scalar FTS: `dx/dt = -|x|^alpha * sign(x)`
- Weights: `r1 = 1` (any positive works; choose 1)
- Degree: `d = alpha - 1 < 0` for `0 < alpha < 1`
- FTS for alpha in (0,1)

### Double Integrator with FTS Control
- Weights: `r = (2, 1)`
- Degree: `d = -1`
- Control: `u = -k1*|x1|^{alpha1}*sign(x1) - k2*|x2|^{alpha2}*sign(x2)` with `alpha1 = 1/3, alpha2 = 1/2`

## References

- Bhat & Bernstein (2005). "Geometric homogeneity with applications to finite-time stability." MCSS 17:101-127.
- Rosier (1992). "Homogeneous Lyapunov function for homogeneous continuous vector field." SCL 19:467-473.
- Polyakov (2012). "Nonlinear feedback design for fixed-time stabilization." IEEE TAC 57:2106-2110.
