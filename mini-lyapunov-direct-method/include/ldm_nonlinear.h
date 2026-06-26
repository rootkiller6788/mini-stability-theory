#ifndef LDM_NONLINEAR_H
#define LDM_NONLINEAR_H
#include "ldm_core.h"

/* ==============================================================
 * ldm_nonlinear.h - Nonlinear System Stability
 *
 * Linearization principle: if the linearization is asymptotically
 * stable, the nonlinear system is locally asymptotically stable.
 *
 * Construction methods for Lyapunov functions:
 *   1. Krasovskii: V = f^T f (for dx/dt = f(x))
 *   2. Variable gradient: assume grad(V) = g(x), enforce curl=0
 *   3. Energy-based: V = kinetic + potential energy
 *   4. Sum-of-squares: polynomial V via SOS decomposition
 *
 * References:
 *   Khalil (2002) Ch.4 - Lyapunov Stability
 *   Parrilo (2000) Structured Semidefinite Programs
 * ============================================================== */

/* Stability analysis */
LDM_StabilityReport ldm_nonlinear_stability(const LDM_System* sys, const double* x0, int n, double t_final, double dt);
LDM_StabilityResult ldm_indirect_method(const LDM_System* sys, const double* x_eq, int n);

/* Jacobian computation */
void ldm_jacobian(const LDM_System* sys, const double* x0, double* J, int n);
void ldm_jacobian_linear(const double* A, double* J, int n);

/* Krasovskii's method */
double ldm_krasovskii_candidate(const LDM_System* sys, int n);
double ldm_krasovskii_derivative(const LDM_System* sys, int n);
bool ldm_krasovskii_check(const LDM_System* sys, int n, int n_samples, double range);

/* Variable gradient method */
bool ldm_variable_gradient_construct(const LDM_System* sys, int n, double* P_out);
bool ldm_variable_gradient_2d(const LDM_System* sys, double* P_out);

/* Energy-based */
double ldm_energy_based(const LDM_System* sys, int n);
double ldm_energy_derivative(const LDM_System* sys, int n);

/* Local/global classification */
bool ldm_is_locally_stable(const LDM_System* sys, const double* x_eq, int n);
bool ldm_is_globally_stable(const LDM_System* sys, int n);
bool ldm_lasalle_check(const LDM_System* sys, int n, double tol);

/* Decay rate */
double ldm_estimate_decay_rate(const LDM_System* sys, const double* x0, int n, double dt, int steps);

/* Comparison */
double ldm_compare_candidates(const LDM_System* sys, int n);

#endif
