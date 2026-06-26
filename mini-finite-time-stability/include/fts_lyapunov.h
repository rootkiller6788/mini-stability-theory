#ifndef FTS_LYAPUNOV_H
#define FTS_LYAPUNOV_H

#include "fts_core.h"

/* ============================================================
 * Lyapunov-Based Finite-Time Stability
 *
 * Lyapunov condition for FTS:
 *   dV/dt <= -c * V(x)^alpha,  with c>0, 0<alpha<1
 * Then settling time bound:
 *   T <= V(x0)^(1-alpha) / (c * (1-alpha))
 *
 * For fixed-time stability:
 *   dV/dt <= -(c1*V^alpha + c2*V^beta), 0<alpha<1<beta
 *   T_max <= 1/(c1*(1-alpha)) + 1/(c2*(beta-1))
 * ============================================================ */

typedef struct {
    double (*V)(FTSState, void*);  /* Lyapunov function */
    void* params;                   /* V function parameters */
    double alpha;                   /* FTS exponent (0<alpha<1) */
    double c;                       /* FTS coefficient (c>0) */
    double beta;                    /* Fixed-time exponent (beta>1) */
    double c2;                      /* Fixed-time second coefficient */
    bool is_fts;                    /* FTS flag */
    bool is_fxt;                    /* Fixed-time flag */
} FTSLyapunovFunc;

typedef struct {
    double* values;                 /* V(t) time series */
    int n;                          /* Number of samples */
    double settling_time;           /* Numerically observed settling time */
    double convergence_rate;        /* Estimated convergence rate */
    double V0;                      /* Initial Lyapunov value */
    double alpha_estimated;         /* Estimated exponent from data */
    FTSConvergenceType type;        /* Detected convergence type */
    double T_upper_bound;           /* Analytical upper bound on T */
} FTSLyapunovAnalysis;

/* --- Lifecycle --- */
FTSLyapunovFunc* fts_lyap_create(double (*V)(FTSState, void*),
                                  void* params, double alpha,
                                  double c, double beta);
void fts_lyap_free(FTSLyapunovFunc* lf);

/* --- Condition Checking --- */
bool   fts_lyap_check_finite_time(FTSLyapunovFunc* lf, FTSSystem* sys,
                                   FTSState x0, double dt, double T);
bool   fts_lyap_check_condition(FTSLyapunovFunc* lf, FTSSystem* sys,
                                 FTSState x, FTSState dx);
int    fts_lyap_verify_condition(FTSLyapunovFunc* lf, FTSSystem* sys,
                                  FTSState x0, double dt, double T,
                                  double* violation_time);

/* --- Lyapunov Derivative --- */
double fts_lyap_derivative(FTSLyapunovFunc* lf, FTSSystem* sys, FTSState x);
double fts_lyap_derivative_numerical(FTSLyapunovFunc* lf, FTSSystem* sys,
                                      FTSState x, double h);

/* --- Analysis --- */
FTSLyapunovAnalysis* fts_lyap_analyze(FTSSystem* sys, FTSState x0,
                                       double T, double dt, int rec);
void fts_lyap_analysis_free(FTSLyapunovAnalysis* a);
bool fts_lyap_is_finite_time_stable(FTSLyapunovAnalysis* a);
bool fts_lyap_is_fixed_time_stable(FTSLyapunovAnalysis* a);

/* --- Bounds --- */
double fts_lyap_settling_time_bound(FTSLyapunovFunc* lf, FTSState x0);
double fts_lyap_fixed_time_bound(FTSLyapunovFunc* lf);
double fts_lyap_settling_time_numerical(FTSLyapunovFunc* lf,
                                         FTSSystem* sys, FTSState x0,
                                         double tol, double dt);

/* --- Utility --- */
void   fts_lyap_print(FTSLyapunovAnalysis* a);
double fts_lyap_homogeneous_degree(FTSLyapunovFunc* lf);
void   fts_lyap_set_fixed_time_params(FTSLyapunovFunc* lf,
                                       double beta, double c2);
double fts_lyap_vdot_bound(FTSLyapunovFunc* lf, FTSSystem* sys,
                            FTSState x, double* V_val);

/* --- Prescribed-Time --- */
double fts_lyap_prescribed_time_gain(FTSLyapunovFunc* lf, double Tp, double t);

/* --- Robust FTS --- */
bool fts_lyap_robust_fts_check(FTSLyapunovFunc* lf, FTSSystem* sys,
                                FTSState x0, double disturbance, double dt, double T);

#endif /* FTS_LYAPUNOV_H */
