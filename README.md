# Mini Stability Theory

A collection of **from-scratch, zero-dependency C implementations** of stability theory for dynamical systems. Covers the complete spectrum from Lyapunov's direct method (1892) to modern nonlinear stability frameworks — Input-to-State Stability (Sontag 1989), absolute stability of Lur'e systems, finite-time convergence, and stability of interconnected systems. Each sub-module maps to MIT, Stanford, UC Berkeley, Cambridge, and other top-tier university control theory courses.

## Sub-Modules

| Sub-Module | Topics | Key Courses |
|------------|--------|-------------|
| [mini-lyapunov-direct-method](mini-lyapunov-direct-method/) | Lyapunov stability theorems, positive-definite functions, quadratic/converse/nonlinear Lyapunov functions, Lyapunov equation, global asymptotic stability, radial unboundedness | MIT 6.241J, Stanford AA 203, Berkeley EE 222 |
| [mini-lasalle-invariance](mini-lasalle-invariance/) | LaSalle invariance principle, Barbashin-Krasovskii theorem, omega-limit sets, invariant sets, region of attraction estimation, convergence analysis, discrete-time LaSalle | MIT 6.241J, Berkeley EE 222, Caltech CDS 140 |
| [mini-input-to-state-stability](mini-input-to-state-stability/) | ISS framework (Sontag 1989), ISS Lyapunov functions, K/KL/K∞ function classes, nonlinear small-gain theorem, ISS for interconnected systems | MIT 6.241J, Berkeley EE 222, Stanford AA 273 |
| [mini-input-output-stability](mini-input-output-stability/) | BIBO stability, L2 gain, small-gain theorem, passivity (dissipativity), conic sectors, Zames-Falb multiplier, circle criterion, Nyquist analysis, structured singular value | MIT 6.241J, Caltech CDS 140, Cambridge 4F3 |
| [mini-absolute-stability](mini-absolute-stability/) | Lur'e system theory, Aizerman conjecture and counterexamples, Kalman-Yakubovich-Popov (KYP) lemma, sector nonlinearities, Lyapunov-based absolute stability tests | MIT 6.241J, Berkeley EE 222, Cambridge 4F3 |
| [mini-circle-popov-criterion](mini-circle-popov-criterion/) | Circle criterion (time-varying nonlinearities), Popov criterion (time-invariant nonlinearities), loop transformation, KYP lemma foundations, Nyquist/Popov plots, gain/phase margin estimation | MIT 6.241J, Cambridge 4F3, ETH 227-0216 |
| [mini-finite-time-stability](mini-finite-time-stability/) | Finite-time convergence (Bhat & Bernstein 2000), fixed-time stability (Polyakov 2012), homogeneous systems, convergence time estimation, sliding surface design, trajectory analysis | MIT 6.832, Stanford AA 203, UC San Diego MAE 281A |
| [mini-stability-of-interconnected](mini-stability-of-interconnected/) | Composite system stability, decentralized Lyapunov functions, M-matrix tests, passivity-based interconnections (Hill-Moylan), small-gain for interconnections, vector Lyapunov functions | MIT 6.241J, Berkeley EE 222, ETH 227-0216 |

## Design Philosophy

- **From-scratch, zero-dependency C** — every module uses only the C standard library (`math.h`, `stdlib.h`, `stdbool.h`); no external numerical libraries, no LAPACK/BLAS
- **Self-contained sub-modules** — each `mini-*` directory is an independent project with its own `include/`, `src/`, `test/`, and build infrastructure
- **Theory-to-code mapping** — header files embed inline derivations, key theorems (Lyapunov 1892, LaSalle 1960, Sontag 1989, Popov 1961, Zames 1966), and textbook references (Khalil, Vidyasagar, Isidori)
- **Course-aligned structure** — every sub-module's API maps to standard graduate course syllabi (MIT 6.241J Dynamic Systems & Control, Berkeley EE 222 Nonlinear Systems, Stanford AA 273 Nonlinear Control), enabling direct use as companion code for lecture study

## Building

Each sub-module is standalone. Build with any C11+ compiler:

```bash
cd mini-lyapunov-direct-method
make
./build/test_runner
```

Or compile manually:

```bash
gcc -std=c11 -O2 -Iinclude -o build/test_ldm src/*.c test/*.c -lm
./build/test_ldm
```

Requires **GCC**, **Clang**, or **MSVC** with C11 support and `-lm` for math linking.

## Project Structure

```
6. mini-stability-theory/
├── mini-absolute-stability/            # Lur'e systems, Aizerman conjecture, KYP lemma, sector nonlinearities
├── mini-circle-popov-criterion/        # Circle criterion, Popov criterion, loop transformation, Nyquist analysis
├── mini-finite-time-stability/         # Finite-time convergence, fixed-time stability, homogeneous systems
├── mini-input-output-stability/        # BIBO, small-gain theorem, passivity, L2 gain, conic sectors
├── mini-input-to-state-stability/      # ISS framework, ISS Lyapunov functions, nonlinear small-gain
├── mini-lasalle-invariance/            # LaSalle principle, Barbashin-Krasovskii, invariant sets, region of attraction
├── mini-lyapunov-direct-method/        # Lyapunov stability theorems, quadratic/converse Lyapunov functions
├── mini-stability-of-interconnected/   # Composite stability, decentralized control, passivity/small-gain for interconnections
├── .gitignore                          # Build artifacts, IDE files, OS files
├── README.md                           # This file (English)
└── README-CN.md                        # Chinese version
```

## License

MIT
