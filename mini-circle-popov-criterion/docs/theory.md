# Circle & Popov Criterion — Theory

## Absolute Stability of Lur'e Systems

A Lur'e system consists of a linear time-invariant (LTI) forward path
G(s) = C(sI-A)^(-1)B + D in negative feedback with a memoryless nonlinearity
phi(sigma) belonging to a sector [alpha, beta]:

  alpha * sigma^2  <=  sigma * phi(sigma)  <=  beta * sigma^2

The problem of **absolute stability** asks: Is the origin globally asymptotically
stable for ALL nonlinearities phi in the given sector?

## Circle Criterion (Zames, 1966)

For a possibly time-varying phi in [alpha, beta], the closed-loop system is
absolutely stable if the Nyquist plot of G(jw) does NOT intersect or encircle
the critical disk D(c, r) where:

  c = -0.5 * (1/alpha + 1/beta)
  r =  0.5 * |1/alpha - 1/beta|

### Special Cases

| Case | Condition | Interpretation |
|------|-----------|---------------|
| alpha = 0 | Re[G(jw)] > -1/beta | Half-plane test |
| alpha = beta | G(jw) avoids -1/alpha | Linear Nyquist |
| alpha = -beta | Disk centered at origin | Symmetric sector |

## Popov Criterion (Popov, 1961)

For a **time-invariant** phi in [0, k], the system is absolutely stable if
there exists eta >= 0 such that for all w:

  Re[(1 + j*w*eta) * G(jw)] + 1/k > 0

Geometrically, in the Popov plot (Re[G] vs w*Im[G]), there must exist a line
through (-1/k, 0) with slope 1/eta such that the plot lies strictly to its right.

### Comparison

| Feature | Circle | Popov |
|---------|--------|-------|
| phi type | Time-varying allowed | Time-invariant only |
| Test | Nyquist vs disk | Popov line |
| Conservatism | More conservative | Less conservative |
| Frequency domain | G(jw) | (1+j*w*eta)*G(jw) |

## KYP Lemma (Kalman-Yakubovich-Popov)

The KYP lemma bridges frequency-domain inequalities and linear matrix
inequalities (LMIs). It establishes that:

  G(jw) + G(-jw)^T + 2*epsilon*I >= 0  (frequency domain)

is equivalent to the existence of P = P^T > 0 satisfying:

  A^T*P + P*A = -L*L^T - epsilon*P
  P*B - C^T = L*W
  D + D^T = W^T*W

This enables state-space (LMI) verification of frequency-domain conditions.

## Sector Theory

### Loop Transformation

Any sector [alpha, beta] can be transformed to canonical [0, k] via:

  G'(s) = G(s) / (1 + alpha*G(s))
  k = beta - alpha

### Passivity Connection

A system is passive iff its nonlinearity lies in sector [0, inf).
Strict passivity corresponds to sector [epsilon, inf) for epsilon > 0.

## References

1. Popov, V.M. (1961). "On absolute stability of non-linear automatic
   control systems." Automatika i Telemekhanika, 22:961-979.
2. Zames, G. (1966). "On the input-output stability of nonlinear
   time-varying feedback systems." IEEE TAC, 11:228-238, 465-476.
3. Kalman, R.E. (1963). "Lyapunov functions for the problem of Lur'e
   in automatic control." PNAS, 49(2):201-205.
4. Khalil, H.K. (2002). Nonlinear Systems (3rd ed.). Prentice Hall.
5. Vidyasagar, M. (1993). Nonlinear Systems Analysis (2nd ed.). Prentice Hall.
6. Haddad, W.M. & Chellaboina, V. (2008). Nonlinear Dynamical Systems
   and Control. Princeton University Press.
7. Boyd, S., El Ghaoui, L., Feron, E., Balakrishnan, V. (1994). Linear
   Matrix Inequalities in System and Control Theory. SIAM.

## Course Alignment

- MIT 6.241J: Dynamic Systems and Control
- Stanford AA203: Optimal and Learning-based Control
- Berkeley EE222: Nonlinear Systems
- Caltech CDS110: Introductory Control Theory
- ETH 227-0216: Nonlinear Control Systems
