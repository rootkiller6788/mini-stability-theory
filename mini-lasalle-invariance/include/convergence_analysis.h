/* convergence_analysis.h - Convergence Analysis for Dynamical Systems
 * L4: LaSalle convergence theorem, rate of convergence estimates
 * L5: Exponential convergence rate, settling time, overshoot, contraction
 * L5: Spectral convergence, monotone convergence, contractive systems
 * Ref: LaSalle (1960), Khalil (2002) Thm 4.1, 4.4, Sontag (1998)
 */

#ifndef CONVERGENCE_ANALYSIS_H
#define CONVERGENCE_ANALYSIS_H
#include "gst_core.h"
#include <stdbool.h>

typedef enum { CONV_EXPONENTIAL, CONV_ASYMPTOTIC, CONV_FINITE_TIME,
               CONV_NONE, CONV_UNKNOWN } ConvergenceType;

typedef struct {
    ConvergenceType type;
    double rate;
    double residual;
    double estimated_order;
    bool is_converged;
    double* error_history;
    int n_error_points;
    double time_constant;
    double settling_time;
    double overshoot;
    double exponential_rate;
} ConvergenceResult;

typedef ConvergenceResult ConvergenceMetrics;

ConvergenceResult* convergence_create(void);
ConvergenceResult* conv_create(void);
void convergence_free(ConvergenceResult* cr);
void conv_free(ConvergenceResult* cr);
int convergence_analyze(ConvergenceResult* cr,
    ODEFunc f, const double* x0, const double* x_eq,
    int n, double dt, double t_final);
double convergence_rate_estimate(const ConvergenceResult* cr);
int convergence_order_estimate(ConvergenceResult* cr);
bool convergence_within_tolerance(const ConvergenceResult* cr, double tol);
void convergence_print(const ConvergenceResult* cr);

int conv_compute_from_trajectory(ConvergenceResult* cr,
    const double* traj, int n_steps, int n_dims,
    const double* x_eq, double dt);
double conv_settling_time_to_tolerance(const ConvergenceResult* cr, double tol_pct);
double conv_time_constant(const ConvergenceResult* cr);
double conv_overshoot(const ConvergenceResult* cr);
bool conv_is_contractive(const ConvergenceResult* cr);
double conv_exponential_rate_est(const ConvergenceResult* cr);

int convergence_spectral_analysis(const double* A_data, int n,
    double* dominant_eigenvalue, double* damping_ratio);

int conv_monotone_convergence_check(const ConvergenceResult* cr);
double conv_decay_envelope(const ConvergenceResult* cr, double alpha);
int conv_convergence_guarantee(ODEFunc f, LyapunovFunc V,
    const double* x0, const double* x_eq, int n,
    double dt, double t_final, double V_bound, double rate_bound);

#endif