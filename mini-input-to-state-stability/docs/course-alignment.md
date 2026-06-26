# Course Alignment — ISS Module

## MIT 6.241J: Dynamic Systems and Control

| Topic | ISS Module Coverage |
|-------|-------------------|
| Lyapunov stability | ISS-Lyapunov function verification (iss_lyapunov.h/c) |
| Input-output stability | ISS gain computation (iss_gain.h/c) |
| Small-gain theorem | Linear and nonlinear small-gain (iss_small_gain.h/c) |
| Interconnected systems | Cascade and feedback ISS (iss_interconnected.h/c) |
| Robustness analysis | ISS verification with bounded inputs (iss_verify.h/c) |

ISS generalizes Lyapunov's direct method to systems with external inputs, making it ideal for the robustness portion of 6.241J (Lectures 14-16).

## Stanford AA203: Optimal and Learning-Based Control

| Topic | Coverage |
|-------|----------|
| Nonlinear stability | ISS verification pipeline |
| Robustness margins | ISS margin computation |
| Verification | Numerical ISS certification |

The ISS module provides the stability verification backbone needed before applying optimal control synthesis.

## ETH 227-0216: System Modeling

| Topic | Coverage |
|-------|----------|
| Comparison functions | K, K-infinity, KL class implementations |
| Dissipation inequalities | ISS Lyapunov condition checking |
| Interconnection analysis | Small-gain and cascade theorems |

The comparison function library (iss_core.h/c) directly supports the abstract stability characterizations taught in ETH's nonlinear systems course.

## Berkeley EE221A/EE222: Nonlinear Systems

| Topic | Coverage |
|-------|----------|
| ISS definition and properties | Core ISS types and system model |
| Converse Lyapunov | ISS ⇔ ISS-Lyapunov (Lean formalization) |
| iISS | Integral ISS analysis |
| Time-delay ISS | Framework for delay robustness |

The Lean 4 formalization (iss_stability.lean) provides proof sketches of the Sontag-Wang theorem and small-gain theorem, complementing the analytical derivations in EE222.

## MIT 16.842: Fundamentals of Systems Engineering

| Topic | Coverage |
|-------|----------|
| System stability verification | Full ISS verification pipeline |
| Interconnected system analysis | Cascade and feedback ISS |
| Robustness quantification | ISS gain estimation |

## Mapping to Classic Textbooks

| Textbook | Chapter | Aligned Functionality |
|----------|---------|----------------------|
| Khalil (2002) Ch. 4 | ISS concepts | iss_core (K/KL functions), iss_lyapunov |
| Sontag (2008) | ISS survey | iss_verify (pipeline), iss_stability.lean |
| Karafyllis & Jiang (2011) | ISS & stabilization | iss_gain, iss_small_gain |
| Isidori (1999) Ch. 10 | ISS & small-gain | iss_interconnected, nonlinear small-gain |

## Implementation Alignment

This module implements a **hands-on** version of the ISS theory:
- **Symbolic**: Lean 4 theorems (iss_stability.lean)
- **Numerical**: C computation (all .c files)
- **Verification**: Random/grid sampling + dissipation checking
- **Interconnection**: Small-gain + cascade + cyclic analysis

The design follows the principle that ISS is "Lyapunov stability + input robustness" and implements both the computational and formal sides.
