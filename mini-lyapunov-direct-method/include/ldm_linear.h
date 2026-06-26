#ifndef LDM_LINEAR_H
#define LDM_LINEAR_H
#include "ldm_core.h"
#include "ldm_quadratic.h"

/* ==============================================================
 * ldm_linear.h - Stability of Linear Systems
 *
 * dx/dt = A x
 *
 * Stability via eigenvalues: Re(lambda_i) < 0 for all i
 * Lyapunov function: V = x^T P x, with P from A^T P + P A = -I
 *
 * Exponential stability: ||x(t)|| <= k * exp(-alpha*t) * ||x(0)||
 * Convergence rate alpha = lambda_min(Q) / (2 * lambda_max(P))
 * ============================================================== */

typedef struct {
    double* eigenvalues;
    int n;
    bool is_stable;
    double spectral_abscissa;
    double* P;
    double convergence_rate;
    double* eigenvectors;
} LDM_LinearAnalysis;

/* Eigenvalue analysis */
LDM_LinearAnalysis* ldm_linear_analyze(const double* A, int n);
void ldm_linear_analysis_free(LDM_LinearAnalysis* la);

/* Stability analysis */
LDM_StabilityReport ldm_linear_stability(const double* A, const double* x0, int n, double t_final, double dt);
bool ldm_linear_is_hurwitz(const double* A, int n);
double ldm_linear_spectral_abscissa(const double* A, int n);

/* Bounds and estimates */
double ldm_linear_overshoot_bound(const double* P, int n);
double ldm_linear_settling_time_estimate(const double* A, int n, double tol);
double ldm_linear_phase_margin_estimate(const double* A, int n);

/* Dynamic properties */
double ldm_linear_damping_ratio(const double* A, int n);
double ldm_linear_natural_frequency(const double* A, int n);

/* Lyapunov matrix */
LDM_QuadraticForm* ldm_linear_lyapunov_matrix(const double* A, int n);
double ldm_linear_exponential_bound(const double* A, int n, double* M, double* alpha);

/* Structural properties */
double ldm_controllability_gramian_trace(const double* A, const double* B, int n, int m);
bool ldm_linear_is_diagonalizable(const double* A, int n);
double ldm_linear_stability_margin(const double* A, int n);
double ldm_linear_peak_bound(const double* A, int n);

#endif
