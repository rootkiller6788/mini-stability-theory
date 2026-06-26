# mini-lasalle-invariance

## LaSalle Invariance Principle

LaSalle invariance principle (1960) extends Lyapunov direct method by relaxing
the strict negative definiteness requirement. When V_dot <= 0 (negative semidefinite),
LaSalle identifies the largest invariant set where trajectories converge.

## Module Status: COMPLETE

- **L1 Definitions**: Complete (12 struct/typedef definitions, Lean definitions)
- **L2 Core Concepts**: Complete (invariant sets, omega-limit, zero-difference)
- **L3 Math Structures**: Complete (Vector, Matrix, RK4, Lyapunov equation, FTS)
- **L4 Fundamental Theorems**: Complete (LaSalle, Barbashin-Krasovskii,
  discrete LaSalle, ranking-decay, singleton invariance -- C + Lean proofs)
- **L5 Algorithms/Methods**: Complete (Poincare map, Floquet, ROA bisection,
  Monte Carlo, convergence rate, invariant iteration, Yoshizawa)
- **L6 Canonical Problems**: Complete (7 systems: pendulum, VdP, Duffing,
  Lorenz, Lotka-Volterra, bistable, FitzHugh-Nagumo)
- **L7 Applications**: Complete (robot PD+gravity LaSalle GAS; L-V ecology;
  neural excitation)
- **L8 Advanced Topics**: Partial (discrete BK, hybrid dwell-time, Yoshizawa,
  gradient detection, monodromy -- 5 of 7 topics)
- **L9 Research Frontiers**: Partial (documented; Lean formal verification bridge)

Line count (include/ + src/): 3095 lines

## Building

    make && make test && make examples

## References
- LaSalle, J.P. (1960). IRE Trans. Circuit Theory, 7:520-527.
- LaSalle & Lefschetz (1961). Stability by Liapunov Direct Method.
- Barbashin & Krasovskii (1952). Dokl. Akad. Nauk SSSR, 86:453-456.
- Khalil, H.K. (2002). Nonlinear Systems (3rd ed.). Prentice Hall.
- Guckenheimer & Holmes (1983). Nonlinear Oscillations.
- LaSalle, J.P. (1976). The Stability of Dynamical Systems. SIAM.
- Yoshizawa, T. (1966). Stability Theory.
- Chesi, G. (2011). Domain of Attraction. Springer.
