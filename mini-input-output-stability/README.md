# Input-Output Stability Theory

## System analysis via input-output mappings

Input-output stability theory analyzes systems through their external behavior -- how inputs map to outputs -- rather than internal state dynamics. This operator-theoretic approach, pioneered by Zames and Sandberg in the 1960s, is the foundation of robust control theory.

### Theoretical Foundation

A system is viewed as an operator H mapping input signals u to output signals y = H(u). Stability is defined in terms of signal norms:

- **L2 stability**: u in L2 implies y in L2 (finite energy in -> finite energy out)
- **L-infinity stability**: u in L-infinity implies y in L-infinity (bounded in -> bounded out)
- **Incremental stability**: small input changes produce small output changes

### The Small-Gain Theorem (Zames, 1966)

For two interconnected systems H1, H2 with L2 gains gamma1, gamma2:

```
gamma1 * gamma2 < 1  =>  feedback interconnection is L2 stable
```

This is the single most important result in input-output stability: if the loop gain is less than unity, the closed-loop system is stable regardless of the internal dynamics of either subsystem.

### Passivity Theorem

A system is passive if the inner product (u, y) >= 0 (net energy absorbed). Key results:

- **Passivity + Passivity = Stable**: Negative feedback interconnection of two passive systems yields a stable closed-loop system
- **Strict passivity**: If either system is strictly passive, the interconnection is asymptotically stable
- **Excess/shortage of passivity**: Quantified using (epsilon, delta)-dissipativity framework

### Conic Sectors and the Circle Criterion

A memoryless nonlinearity phi(.) in sector [a, b] satisfies:

```
a * sigma^2 <= sigma * phi(sigma) <= b * sigma^2
```

The feedback interconnection of a linear system G(s) with phi in sector [a, b] is stable if the Nyquist plot of G(jw) avoids the critical disk centered at -(1/a+1/b)/2 with radius (1/a-1/b)/2.

### Operator Gain Computation

| Operator | L2 Gain | Method |
|----------|---------|--------|
| Linear time-invariant | sup_omega sigma_max(H(jw)) | H-infinity norm |
| Static nonlinearity | sup |f(x)| / |x| | Lipschitz constant |
| Cascaded operators | gamma1 * gamma2 | Multiplicative |
| Parallel operators | gamma1 + gamma2 | Additive (worst case) |

### Passivity Indices

A system is Input Feedforward Passive (IFP) with index nu if:
  (u, y) >= nu * ||u||^2

A system is Output Feedback Passive (OFP) with index rho if:
  (u, y) >= rho * ||y||^2

Combined IFP(nu) + OFP(rho) interconnection is stable if nu1 + rho2 > 0.

### Implementation

| File | Purpose |
|------|---------|
| `io_stability.h/c` | Core types, gain computation, operator norms |
| `small_gain.h/c` | Small-gain theorem, loop gain analysis |
| `passivity.h/c` | Passivity verification, storage functions |
| `conic_sectors.h/c` | Circle criterion, Popov criterion, sector bounds |
| `operator_norms.h/c` | L2, L-infinity, incremental gain computation |

### Key Functions

| Function | Description |
|----------|-------------|
| `ios_gain_l2` | Compute L2 gain of an operator via simulation |
| `ios_small_gain_test` | Check small-gain condition for interconnection |
| `ios_passivity_index` | Compute passivity indices (epsilon, delta) |
| `ios_conic_sector_fit` | Fit nonlinearity to conic sector bounds |
| `ios_circle_criterion_check` | Verify circle criterion stability |
| `ios_incremental_gain` | Estimate incremental Lipschitz gain |

### References

- Zames, G. (1966). On the input-output stability of time-varying nonlinear feedback systems. IEEE Trans. Automatic Control, 11:228-238, 465-476.
- Desoer, C.A. & Vidyasagar, M. (1975). Feedback Systems: Input-Output Properties. Academic Press.
- Van der Schaft, A. (2000). L2-Gain and Passivity Techniques in Nonlinear Control. Springer.
- Khalil, H.K. (2002). Nonlinear Systems (3rd ed.). Prentice Hall.
- Safonov, M.G. (1980). Stability and Robustness of Multivariable Feedback Systems. MIT Press.
- Sandberg, I.W. (1964). On the L2-boundedness of solutions of nonlinear functional equations. Bell System Tech. J.

### Examples

1. **example1**: Small-gain verification on interconnected systems
2. **example2_small_gain**: Loop gain computation and stability test
3. **example3_passivity**: Passivity-based feedback stability verification

### Building

```
make && make test && make examples
```
