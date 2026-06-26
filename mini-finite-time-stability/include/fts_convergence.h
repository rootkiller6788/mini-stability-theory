#ifndef FTS_CONVERGENCE_H
#define FTS_CONVERGENCE_H

#include "fts_core.h"

/* ============================================================
 * Convergence Analysis and Settling Time Estimation
 *
 * Numerical methods for detecting and characterizing convergence.
 *
 * Key concepts:
 * - Settling time: T_s = inf{t>=0 : ||x(t)|| < epsilon}
 * - Practical FTS: convergence to a neighborhood
 * - Convergence rate estimation via log-linear regression
 * - Multiple convergence detection criteria
 * ============================================================ */

typedef struct {
    double* T_estimates;  /* Array of settling time estimates */
    int n;                /* Number of estimates */
    double T_min;         /* Minimum settling time */
    double T_max;         /* Maximum settling time */
    double T_mean;        /* Mean settling time */
    double T_std;         /* Standard deviation */
    bool is_finite;       /* Finite settling time detected */
    FTSConvergenceType type; /* Convergence type */
    double rate;          /* Estimated convergence rate */
    double overshoot;     /* Maximum overshoot */
    double steady_error;  /* Steady-state error */
} FTSSettlingTime;

typedef struct {
    double alpha;         /* Lyapunov exponent */
    double beta;          /* Second exponent (for FxT) */
    double T_upper;       /* Upper bound on T */
    double T_lower;       /* Lower bound on T */
    bool is_tight;        /* Whether bounds are tight */
    double confidence;    /* Confidence in bound (0-1) */
} FTSTimeBound;

/* --- Lifecycle --- */
FTSSettlingTime* fts_st_create(int n_estimates);
void fts_st_free(FTSSettlingTime* st);

/* --- Computation --- */
int  fts_st_compute(FTSSystem* sys, FTSState x0, double T,
                    double dt, double tol, FTSSettlingTime* st);
double fts_st_empirical(FTSSystem* sys, FTSState x0,
                         double dt, double tol);
int  fts_st_multi_compute(FTSSystem* sys, FTSState* x0_array,
                           int n_initial, double T, double dt,
                           double tol, FTSSettlingTime* st);

/* --- Classification --- */
FTSConvergenceType fts_st_classify(FTSSettlingTime* st);
bool fts_st_is_finite(FTSSettlingTime* st);
bool fts_st_is_fixed_time(FTSSettlingTime* st,
                           FTSState* x0_array, int n);

/* --- Rate Estimation --- */
double fts_st_convergence_rate(FTSSettlingTime* st);
double fts_st_exponential_rate(FTSSystem* sys, FTSState x0,
                                double T, double dt);

/* --- Bounds --- */
FTSTimeBound* fts_st_bound(FTSSystem* sys, FTSState x0,
                            double alpha, double c, double beta);
void fts_tb_free(FTSTimeBound* tb);
double fts_st_residual_energy(FTSSystem* sys, FTSState x0, double t);
double fts_st_residual_energy_grid(FTSSystem* sys, FTSState x0,
                                    double* t_grid, int n);

/* --- Comparison --- */
int  fts_st_compare(FTSSettlingTime* a, FTSSettlingTime* b);
bool fts_st_meets_spec(FTSSettlingTime* st, double T_spec);

/* --- Output --- */
void fts_st_print(FTSSettlingTime* st);
void fts_tb_print(FTSTimeBound* tb);

/* --- Practical Stability --- */
bool fts_st_practical_finite_time(FTSSystem* sys, FTSState x0,
                                   double epsilon, double T,
                                   double dt);
double fts_st_practical_settling_time(FTSSystem* sys, FTSState x0,
                                       double epsilon, double dt);

#endif /* FTS_CONVERGENCE_H */
