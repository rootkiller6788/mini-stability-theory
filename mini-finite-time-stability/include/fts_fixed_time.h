#ifndef FTS_FIXED_TIME_H
#define FTS_FIXED_TIME_H

#include "fts_core.h"

/* ============================================================
 * Fixed-Time Stability (Polyakov, 2012)
 *
 * Fixed-time stability: settling time is bounded uniformly
 * for ALL initial conditions.
 *
 * Lyapunov condition:
 *   dV/dt <= -(c1*V^alpha + c2*V^beta), 0<alpha<1<beta
 * Then: T_max <= 1/(c1*(1-alpha)) + 1/(c2*(beta-1))
 *
 * Compare:
 *   Asymptotic:    T depends on x0, unbounded
 *   Exponential:   T ~ -log(eps)/lambda
 *   Finite-time:   T <= V0^(1-alpha)/(c*(1-alpha))
 *   Fixed-time:    T_max uniform, independent of x0
 *
 * Ref: Polyakov, A. (2012). "Nonlinear feedback design for
 * fixed-time stabilization." IEEE TAC, 57:2106-2110.
 * ============================================================ */

typedef struct {
    double T_max;          /* Maximum settling time bound */
    double alpha;          /* FTS exponent (0<alpha<1) */
    double beta;           /* FxT exponent (beta>1) */
    double p;              /* First coefficient c1 */
    double q;              /* Second coefficient c2 */
    double gamma;          /* Weighting factor gamma>0 */
    bool is_fixed_time;    /* Fixed-time flag */
    bool is_uniform;       /* Uniform bound property */
} FTSFixedTimeParams;

typedef struct {
    FTSFixedTimeParams params;   /* Fixed-time parameters */
    double* settling_times;      /* Array of settling times */
    int n_tests;                 /* Number of tests */
    double T_observed_max;       /* Maximum observed settling time */
    double T_observed_min;       /* Minimum observed settling time */
    double T_mean;               /* Mean settling time */
    bool verified;               /* Verification result */
    double violation_margin;     /* How far from bound */
} FTSFixedTimeTest;

/* --- Lifecycle --- */
FTSFixedTimeParams* fts_ft_params_create(double alpha, double beta,
                                          double p, double q);
void fts_ft_params_free(FTSFixedTimeParams* ftp);

/* --- Fixed-Time Detection --- */
bool fts_ft_is_fixed_time(FTSFixedTimeParams* ftp);
bool fts_ft_is_uniform(FTSFixedTimeParams* ftp);

/* --- Bounds --- */
double fts_ft_settling_bound(FTSFixedTimeParams* ftp);
double fts_ft_settling_bound_general(FTSFixedTimeParams* ftp, FTSState x0);
double fts_ft_lyapunov_bound(FTSFixedTimeParams* ftp, FTSState x0);
double fts_ft_uncertainty_bound(FTSFixedTimeParams* ftp,
                                 double delta_c1, double delta_c2);

/* --- Testing --- */
FTSFixedTimeTest* fts_ft_test(FTSSystem* sys, FTSFixedTimeParams* ftp,
                               int n_tests, double T_sim, double dt);
void fts_ft_test_free(FTSFixedTimeTest* ftt);
bool fts_ft_verify(FTSFixedTimeTest* ftt);
double fts_ft_max_settling(FTSFixedTimeTest* ftt);

/* --- Comparison --- */
bool fts_ft_compare_finite_vs_fixed(FTSFixedTimeParams* ftp,
                                     double finite_T);
double fts_ft_robustness_margin(FTSFixedTimeParams* ftp,
                                 double disturbance);

/* --- Output --- */
void fts_ft_print(FTSFixedTimeTest* ftt);
void fts_ft_params_print(FTSFixedTimeParams* ftp);

/* --- Prescribed-Time Stability --- */
bool fts_ft_is_prescribed_time(FTSFixedTimeParams* ftp, double Tp);
double fts_ft_prescribed_gain(double Tp, double t);

/* --- Implicit Lyapunov Function Methods --- */
double fts_ft_implicit_lyap_bound(FTSFixedTimeParams* ftp, FTSState x0);

#endif /* FTS_FIXED_TIME_H */
