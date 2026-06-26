# Stability of Interconnected Systems - Theory

## Small-Gain Theorem (Zames 1966)

For feedback interconnection of two systems with L2 gains gamma1, gamma2:
If gamma1 * gamma2 < 1, the interconnection is finite-gain L2-stable.

Generalization to N subsystems: if the spectral radius of the gain matrix is < 1.

## Passivity Theorem (Willems 1972)

Feedback interconnection of two passive systems is L2-stable if:
- Both are passive
- At least one is strictly passive (input or output)

Generalized via dissipativity framework with (Q,S,R) supply rates.

## Lyapunov Methods (Siljak 1978)

### Composite Lyapunov Functions
V(x) = sum d_i * V_i(x_i) with d_i > 0

### Vector Lyapunov Functions
V_dot <= M * V (componentwise comparison)

### M-Matrix Stability
If M is an M-matrix, the origin is asymptotically stable.

## Gershgorin Circle Theorem
All eigenvalues of A lie within union of discs centered at a_ii with radius sum_{j!=i} |a_ij|.

## Decentralized Stability
Each subsystem stabilized by local controller. Stability guaranteed if interaction bound < stability margin of isolated subsystems.