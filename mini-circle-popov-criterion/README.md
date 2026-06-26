# Circle & Popov Criterion — mini-circle-popov-criterion

## Absolute Stability of Lur'e Systems

The Circle and Popov criteria are frequency-domain tests for absolute
stability — global asymptotic stability of a Lur'e system for ALL
nonlinearities within a given sector, not just a specific one.

### Lur'e System Structure

A Lur'e system consists of a linear time-invariant (LTI) forward path
G(s) = C(sI-A)^(-1)B + D in feedback with a memoryless nonlinearity
phi(sigma) satisfying a sector condition:



### Circle Criterion (Time-Varying phi)

For a nonlinearity phi(t, sigma) in sector [alpha, beta] that may be
time-varying, the closed-loop system is absolutely stable if the
Nyquist plot of G(jw) does not intersect or encircle the critical disk
D(alpha, beta) centered at -(1/alpha + 1/beta)/2 with radius
(1/alpha - 1/beta)/2.

Special cases:
- alpha = 0: Disk becomes half-plane Re[G(jw)] > -1/beta
- alpha = beta: Disk reduces to point -1/alpha (linear gain)
- alpha = -beta: Disk is symmetric about the origin

### Popov Criterion (Time-Invariant phi)

For a time-invariant nonlinearity phi(sigma) in sector [0, k], the
system is absolutely stable if there exists eta >= 0 such that:



Geometrically: plot Re[G(jw)] vs w*Im[G(jw)] (Popov plot). If there
exists a line through -1/k with slope 1/eta that the Popov plot lies
strictly to the right of, the system is absolutely stable.

### Kalman-Yakubovich-Popov (KYP) Lemma

The KYP lemma establishes the equivalence between frequency-domain
inequalities and linear matrix inequalities (LMIs):

Frequency condition: G(jw) + G(-jw)' + 2/k*I > 0 for all w

is equivalent to the existence of P = P' > 0 and q satisfying:


This bridges classical (frequency-domain) and modern (state-space)
control theory.

### Comparison: Circle vs Popov

| Feature          | Circle Criterion           | Popov Criterion              |
|------------------|---------------------------|------------------------------|
| Nonlinearity     | Time-varying phi(t,sigma) | Time-invariant phi(sigma)    |
| Sector           | [alpha, beta]             | [0, k] (or shifted)          |
| Test             | Nyquist avoids circle     | Popov line condition         |
| Conservatism     | More conservative         | Less conservative for TI     |
| Multiplier       | None (pure circle)        | (1 + eta*s) multiplier       |

## Implementation

### File Structure

| File                      | Lines | Purpose                                    |
|---------------------------|-------|--------------------------------------------|
|  |  170 | Circle criterion API, Nyquist data types  |
|   |  137 | Popov criterion API, Popov plot types     |
|         |  160 | KYP lemma API, state-space types          |
|       |  164 | Lur'e system API, unified stability tests |
|     |  112 | Sector theory API, nonlinearity functions |
|      |  385 | Transfer function eval, Nyquist, circle   |
|       |  438 | Popov plot, line test, multiplier search  |
|             |  458 | State-space, Leverrier, Lyapunov solver   |
|           |  519 | Lur'e system, margins, loop transform     |
|         |  503 | Sector verification, nonlinearities       |
|       |  130 | Lean 4 formalization                      |
|         |  216 | 61 assert()-based test cases              |
|  |   44 | Circle criterion demo                     |
|   |   54 | Popov criterion with optimal eta          |
|  | 86 | Full comparison via state-space          |
|              |  100 | Theory reference                          |
|           |   64 | Course alignment                          |

**Total: include/ + src/ = 3046 lines, src/*.c = 2303 lines**

### Key Functions (>= 20 exported)

1.  /  — Transfer function lifecycle
2.  — Complex frequency evaluation
3.  /  — Nyquist data
4.  — Circle criterion verification
5.  /  — Critical disk
6.  — Winding number test
7.  — Popov criterion verification
8.  — Optimal multiplier search
9.  — Circle vs Popov comparison
10.  /  — State-space lifecycle
11.  — KYP LMI solver
12.  — Leverrier algorithm
13.  — Lyapunov equation solver
14.  /  — Lur'e system lifecycle
15.  — Unified stability test
16.  /  — Margins
17.  /  — Verification
18.  /  /  — Transformations
19.  /  /  — Standard NL
20.  — Maximal sector search

### Building and Testing

gcc -std=c11 -Wall -Wextra -O2 -Iinclude -o build/test_circle tests/test_circle.c -Lbuild -lcircle -lm
Running build/test_circle...
=== Circle-Popov Criterion Test Suite ===
  Test 1: Transfer function...
  Test 2: Critical disk...
  Test 3: Half-plane case...
  Test 4: Nyquist data...
  Test 5: Log-swept Nyquist...
  Test 6: Popov criterion...
  Test 7: Optimal Popov slope...
  Test 8: Popov line test...
  Test 9: Loop transformation...
  Test 10: State-space and KYP...
  Test 11: SS-to-TF conversion...
  Test 12: Lyapunov solver...
  Test 13: Passivity check...
  Test 14: Lur'e system...
  Test 15: Lur'e stability...
  Test 16: Unified stability test...
  Test 17: Gain/phase margins...
  Test 18: H-infinity norm...
  Test 19: Sector theory...
  Test 20: Sector transformations...
  Test 21: Circle vs Popov comparison...
  Test 22: Max allowed gain...

=== ALL TESTS PASSED ===
=== Test Results: 1 passed, 0 failed ===
gcc -std=c11 -Wall -Wextra -O2 -Iinclude -o build/example1_circle examples/example1_circle.c -Lbuild -lcircle -lm
gcc -std=c11 -Wall -Wextra -O2 -Iinclude -o build/example2_popov examples/example2_popov.c -Lbuild -lcircle -lm
gcc -std=c11 -Wall -Wextra -O2 -Iinclude -o build/example3_criterion examples/example3_criterion.c -Lbuild -lcircle -lm
=== Running build/example1_circle ===
=== Example 1: Circle Criterion ===
G(s) = 1/(s+1)
  Sector [0, 0.5]: STABLE (margin=2.0001)
  Sector [0, 1.0]: STABLE (margin=1.0001)
  Sector [0, 2.0]: STABLE (margin=0.5001)
  Sector [0, 5.0]: STABLE (margin=0.2001)
  Sector [0, 10.0]: STABLE (margin=0.1001)
  Sector [0, 50.0]: STABLE (margin=0.0201)
  Sector [0, 100.0]: STABLE (margin=0.0101)

Critical disk for [0.5, 3.0]:
  Center = -1.1667
  Radius = 0.8333
  Stability: ABSOLUTELY STABLE (margin=0.3335)
=== Running build/example2_popov ===
=== Example 2: Popov Criterion ===
G(s) = 1/(s^2 + s + 1)

k          Circle          Popov           PopovEta  
-----------------------------------------------
k=0.5      STABLE          STABLE          eta=0.00
k=1.0      STABLE          STABLE          eta=0.00
k=2.0      STABLE          STABLE          eta=0.00
k=3.0      STABLE          STABLE          eta=0.00
k=5.0      inconclusive    STABLE          eta=2.04
k=10.0     inconclusive    STABLE          eta=2.04

--- Comparison ---
Comparison result: Circle=OK Popov=OK (eta=0.00)
=== Running build/example3_criterion ===
=== Example 3: Full Criterion Comparison ===

State-space system (n=2):
  A = [0.00 1.00; -1.00 -0.50]
  B = [0.00; 1.00]
  C = [2.00 1.00]
  G(s) = (0.00 s^2 + 1.00 s + 2.00) / (1.00 s^2 + 0.50 s + 1.00)

-------------------------------------------
k        Circle     Popov      KYP       
-------------------------------------------
k=0.5    OK         OK         OK        
k=1.0    OK         OK         OK        
k=2.0    FAIL       OK         OK        
k=5.0    FAIL       OK         OK        
k=10.0   FAIL       OK         OK        
k=20.0   FAIL       OK         OK        

--- Max Sector via Circle Criterion ---
Max k (circle): 1.50

--- Sector Transformations ---
Shift [1,4] by -1 => [0.0, 3.0]
Scale [1,4] by 0.5 => [0.5, 2.0]

Done.
rm -rf build

### Requirements

- GCC (C11) or compatible C compiler
- GNU Make
- POSIX-compatible shell (for test runner)
- Lean 4 (optional, for formal verification)

### References

- Popov, V.M. (1961). "On absolute stability of non-linear automatic
  control systems." Automatika i Telemekhanika, 22:961-979.
- Zames, G. (1966). "On the input-output stability of nonlinear
  feedback systems." IEEE TAC, 11:228-238, 465-476.
- Kalman, R.E. (1963). "Lyapunov functions for the problem of Lur'e
  in automatic control." PNAS, 49(2):201-205.
- Khalil, H.K. (2002). Nonlinear Systems (3rd ed.). Prentice Hall.
- Vidyasagar, M. (1993). Nonlinear Systems Analysis (2nd ed.).
  Prentice Hall.
- Haddad, W.M. & Chellaboina, V. (2008). Nonlinear Dynamical Systems
  and Control. Princeton University Press.
- Boyd, S. et al. (1994). Linear Matrix Inequalities in System and
  Control Theory. SIAM.

### Examples

**Example 1 — Circle Criterion**: Tests sector [0, k] for various k on
G(s) = 1/(s+1). Shows the half-plane case.

**Example 2 — Popov Criterion**: Compares circle vs Popov on
G(s) = 1/(s^2 + s + 1). Demonstrates Popov's reduced conservatism
for time-invariant nonlinearities.

**Example 3 — Full Comparison**: State-space to transfer function
conversion, KYP lemma, sector transformations, and robust margin search.

### Course Alignment

- MIT 6.241J — Dynamic Systems and Control
- Stanford AA203 — Optimal and Learning-based Control
- Berkeley EE222 — Nonlinear Systems
- Caltech CDS110 — Introductory Control Theory
- ETH 227-0216 — Nonlinear Control Systems
