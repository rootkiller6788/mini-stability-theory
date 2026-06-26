#ifndef LDM_QUADRATIC_H
#define LDM_QUADRATIC_H
#include "ldm_core.h"

/* ==============================================================
 * ldm_quadratic.h - Quadratic Lyapunov Functions
 *
 * V(x) = x^T P x,  P = P^T > 0
 *
 * dV/dt = x^T (A^T P + P A) x  (linear systems)
 * Lyapunov equation: A^T P + P A = -Q  (Q > 0)
 *
 * Solution methods:
 *   1. Closed-form for n=1,2,3
 *   2. Bartels-Stewart for n>3 (Schur decomposition)
 *   3. Iterative (Kleinman) for large n
 *
 * References:
 *   Gajic & Qureshi (1995) Lyapunov Matrix Equation
 *   Bartels & Stewart (1972) ACM Algorithm 432
 * ============================================================== */

/* Lyapunov equation solvers */
LDM_QuadraticForm* ldm_lyap_solve_2x2(const double* A, const double* Q);
LDM_QuadraticForm* ldm_lyap_solve_3x3(const double* A, const double* Q);
LDM_QuadraticForm* ldm_lyap_solve_kleinman(const double* A, const double* Q, int n, int max_iter, double tol);
LDM_QuadraticForm* ldm_lyap_solve(const double* A, const double* Q, int n);

/* Derivative and verification */
double ldm_quadratic_derivative_linear(const double* A, const double* P, const double* x, int n);
bool ldm_quadratic_is_lyapunov(const double* A, const double* P, int n);
bool ldm_quadratic_is_symmetric(const double* P, int n);

/* Convergence and decay rates */
double ldm_quadratic_convergence_rate(const double* A, const double* P, const double* Q, int n);
double ldm_quadratic_decay_rate(const double* P, int n);
double ldm_quadratic_condition_number(const double* P, int n);

/* Residual analysis and solution bounds */
double ldm_lyap_residual_relative(const double* A, const double* P, const double* Q, int n);
void ldm_lyap_solution_bounds(const double* A, const double* Q, int n, double* p_min_est, double* p_max_est);
bool ldm_lyap_verify_bounds(const double* A, const double* P, const double* Q, int n);

#endif
