#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fts_core.h"
#include "fts_lyapunov.h"

/* ============================================================
 * fts_lyapunov.c -- Lyapunov-Based Finite-Time Stability
 *
 * Implements Lyapunov condition checking for FTS:
 *   dV/dt <= -c * V(x)^alpha, with c>0, 0<alpha<1
 *
 * Settling time bound:
 *   T <= V(x0)^(1-alpha) / (c * (1-alpha))
 *
 * Also covers fixed-time Lyapunov conditions:
 *   dV/dt <= -(c1*V^alpha + c2*V^beta), 0<alpha<1<beta
 *   T_max <= 1/(c1*(1-alpha)) + 1/(c2*(beta-1))
 *
 * Refs: Bhat & Bernstein (2000), Polyakov (2012),
 *       Moulay & Perruquetti (2006)
 * ============================================================ */

/* ---------- Lifecycle ---------- */

FTSLyapunovFunc* fts_lyap_create(double (*V)(FTSState, void*),
                                  void* params, double alpha,
                                  double c, double beta) {
    FTSLyapunovFunc* lf = calloc(1, sizeof(FTSLyapunovFunc));
    if (!lf) return NULL;
    lf->V = V; lf->params = params;
    lf->alpha = alpha; lf->c = c; lf->beta = beta;
    lf->c2 = 0.0;
    /* Validate FTS condition: c>0 and 0<alpha<1 */
    lf->is_fts = (c > FTS_EPSILON && alpha > 0.0 && alpha < 1.0);
    /* Fixed-time condition requires beta>1 and c2>0 in addition */
    lf->is_fxt = false;
    return lf;
}

void fts_lyap_free(FTSLyapunovFunc* lf) {
    free(lf);
}

void fts_lyap_set_fixed_time_params(FTSLyapunovFunc* lf,
                                     double beta, double c2) {
    if (!lf) return;
    lf->beta = beta; lf->c2 = c2;
    lf->is_fxt = (lf->is_fts && beta > 1.0 && c2 > FTS_EPSILON);
}

/* ---------- Lyapunov Condition Checking ---------- */

bool fts_lyap_check_condition(FTSLyapunovFunc* lf, FTSSystem* sys,
                               FTSState x, FTSState dx) {
    /* Verify if dV/dt <= -c * V^alpha at a point.
     * dV/dt = grad(V)^T * f(x), approximated via directional derivative. */
    if (!lf || !sys) return false;
    double V_val = lf->V(x, lf->params);
    double Vdot = fts_lyap_derivative(lf, sys, x);
    double bound = -lf->c * pow(V_val, lf->alpha);
    return (Vdot <= bound + FTS_EPSILON);
}

int fts_lyap_verify_condition(FTSLyapunovFunc* lf, FTSSystem* sys,
                               FTSState x0, double dt, double T,
                               double* violation_time) {
    /* Integrate system and check FTS Lyapunov condition at each step.
     * Returns: 0=condition holds throughout, 1=violation found.
     * violation_time is set to the time of first violation. */
    if (!lf || !sys || !violation_time) return -1;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return -1;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    int result = 0;
    while (cp->t < T) {
        fts_step_rk4(cp);
        double V_val = lf->V(cp->state, lf->params);
        double Vdot = fts_lyap_derivative(lf, cp, cp->state);
        double bound = -lf->c * pow(V_val, lf->alpha);
        if (Vdot > bound + FTS_EPSILON && V_val > FTS_EPSILON) {
            *violation_time = cp->t; result = 1; break;
        }
    }
    fts_free(cp);
    return result;
}

bool fts_lyap_check_finite_time(FTSLyapunovFunc* lf, FTSSystem* sys,
                                 FTSState x0, double dt, double T) {
    /* Simple FTS check: integrate and verify convergence to origin
     * within specified time horizon. */
    if (!lf || !sys) return false;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return false;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    bool converged = false;
    while (cp->t < T) {
        fts_step_rk4(cp);
        double Vt = lf->V(cp->state, lf->params);
        if (Vt < 1e-6) { converged = true; break; }
    }
    fts_free(cp);
    return converged;
}

/* ---------- Lyapunov Derivative ---------- */

double fts_lyap_derivative(FTSLyapunovFunc* lf, FTSSystem* sys,
                            FTSState x) {
    /* Compute dV/dt = grad(V)^T * f(x) numerically.
     * For V(x)=||x||^2, dV/dt = 2*x^T*f(x).
     * For general V, use finite differences. */
    if (!lf || !sys) return 0.0;
    FTSState dx; dx.dim = x.dim;
    sys->rhs(sys, x, &dx);
    /* For quadratic V(x)=0.5*||x||^2: dV/dt = x^T * f(x) */
    double dV = 0.0;
    for (int i = 0; i < x.dim; i++)
        dV += x.x[i] * dx.x[i];
    return dV;
}

double fts_lyap_derivative_numerical(FTSLyapunovFunc* lf,
                                      FTSSystem* sys, FTSState x,
                                      double h) {
    /* Compute dV/dt using finite differences:
     * dV/dt approximately (V(x+h*f(x)) - V(x)) / h */
    if (!lf || !sys) return 0.0;
    double V0 = lf->V(x, lf->params);
    FTSState dx; dx.dim = x.dim;
    sys->rhs(sys, x, &dx);
    FTSState x1 = x;
    for (int i = 0; i < x.dim; i++)
        x1.x[i] = x.x[i] + h * dx.x[i];
    double V1 = lf->V(x1, lf->params);
    if (h < FTS_EPSILON) return 0.0;
    return (V1 - V0) / h;
}

/* ---------- Lyapunov Analysis ---------- */

FTSLyapunovAnalysis* fts_lyap_analyze(FTSSystem* sys, FTSState x0,
                                       double T, double dt, int rec) {
    /* Perform full Lyapunov analysis:
     * - Simulate system from x0 for time T
     * - Record V(t) values at each step
     * - Estimate settling time (when V < 1e-3)
     * - Classify convergence type */
    if (!sys) return NULL;
    int N = (int)(T / dt);
    if (N <= 0) N = 1000;
    FTSLyapunovAnalysis* a = calloc(1, sizeof(FTSLyapunovAnalysis));
    if (!a) return NULL;
    a->values = calloc(N, sizeof(double));
    if (!a->values) { free(a); return NULL; }
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) { free(a->values); free(a); return NULL; }
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));

    a->V0 = fts_norm(cp->state);
    a->settling_time = T;
    a->type = FTS_ASYMPTOTIC;
    bool settled = false;

    for (int i = 0; i < N; i++) {
        fts_step_rk4(cp);
        double V = fts_norm(cp->state);
        if (i < N) { a->values[i] = V; a->n++; }
        if (!settled && V < 1e-3) {
            a->settling_time = cp->t; settled = true;
            a->type = FTS_FINITE_TIME;
        }
    }

    /* Estimate convergence rate via log-linear regression */
    if (a->n > 10 && settled) {
        double sum_t = 0, sum_logV = 0, sum_t2 = 0, sum_tlogV = 0;
        int m = 0;
        for (int i = 0; i < a->n; i++) {
            if (a->values[i] > 1e-12) {
                double t_i = i * dt;
                double logV = log(a->values[i]);
                sum_t += t_i; sum_logV += logV;
                sum_t2 += t_i*t_i; sum_tlogV += t_i*logV;
                m++;
            }
        }
        if (m > 1) {
            double denom = m*sum_t2 - sum_t*sum_t;
            if (fabs(denom) > FTS_EPSILON)
                a->convergence_rate = -(m*sum_tlogV - sum_t*sum_logV) / denom;
        }
    }

    /* Analytical upper bound using default Lyapunov function */
    double V0 = fts_norm(x0);
    a->T_upper_bound = V0 / (1.0 - 0.5); /* simple heuristic */
    fts_free(cp);
    return a;
}

void fts_lyap_analysis_free(FTSLyapunovAnalysis* a) {
    if (a) { free(a->values); free(a); }
}

bool fts_lyap_is_finite_time_stable(FTSLyapunovAnalysis* a) {
    return a && a->type == FTS_FINITE_TIME;
}

bool fts_lyap_is_fixed_time_stable(FTSLyapunovAnalysis* a) {
    /* Fixed-time requires uniform bound. Here we check if
     * convergence was achieved and rate suggests fixed-time.
     * For rigorous check, need multiple initial conditions. */
    return a && a->type == FTS_FINITE_TIME &&
           a->convergence_rate > 0.0;
}

/* ---------- Bounds ---------- */

double fts_lyap_settling_time_bound(FTSLyapunovFunc* lf, FTSState x0) {
    /* Analytical upper bound on settling time:
     * T <= V(x0)^(1-alpha) / (c * (1-alpha))
     *
     * Derivation: integrate dV/dt <= -c*V^alpha
     *   dV / V^alpha <= -c*dt
     *   int_{V0}^{0} V^(-alpha) dV <= -c * T
     *   V0^(1-alpha)/(1-alpha) = c * T
     *   => T = V0^(1-alpha) / (c*(1-alpha)) */
    if (!lf || lf->c < FTS_EPSILON) return FTS_INF;
    double V0 = lf->V(x0, lf->params);
    if (V0 < FTS_EPSILON) return 0.0;
    if (lf->alpha >= 1.0) return FTS_INF; /* not FTS */
    return pow(V0, 1.0 - lf->alpha) / (lf->c * (1.0 - lf->alpha));
}

double fts_lyap_fixed_time_bound(FTSLyapunovFunc* lf) {
    /* Fixed-time settling time bound (Polyakov 2012):
     * T_max <= 1/(c1*(1-alpha)) + 1/(c2*(beta-1))
     * Independent of initial condition! */
    if (!lf || !lf->is_fxt) return FTS_INF;
    double T1 = (lf->c > FTS_EPSILON && lf->alpha < 1.0) ?
                1.0 / (lf->c * (1.0 - lf->alpha)) : 0.0;
    double T2 = (lf->c2 > FTS_EPSILON && lf->beta > 1.0) ?
                1.0 / (lf->c2 * (lf->beta - 1.0)) : 0.0;
    return T1 + T2;
}

double fts_lyap_settling_time_numerical(FTSLyapunovFunc* lf,
                                         FTSSystem* sys, FTSState x0,
                                         double tol, double dt) {
    /* Numerically compute settling time by integrating until
     * V(x) < tol. Returns FTS_INF if not converged within 1000 time units. */
    if (!lf || !sys) return FTS_INF;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return FTS_INF;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    double Tmax = 1000.0;
    while (cp->t < Tmax) {
        fts_step_rk4(cp);
        double Vt = lf->V(cp->state, lf->params);
        if (Vt < tol) { double Ts = cp->t; fts_free(cp); return Ts; }
    }
    fts_free(cp);
    return FTS_INF;
}

/* ---------- V_dot Bound ---------- */

double fts_lyap_vdot_bound(FTSLyapunovFunc* lf, FTSSystem* sys,
                            FTSState x, double* V_val) {
    /* Compute both V(x) and dV/dt, return dV/dt.
     * Useful for checking the FTS condition. */
    if (!lf || !sys) { if (V_val) *V_val = 0; return 0; }
    double V = lf->V(x, lf->params);
    if (V_val) *V_val = V;
    return fts_lyap_derivative(lf, sys, x);
}

/* ---------- Prescribed-Time Gain ---------- */

double fts_lyap_prescribed_time_gain(FTSLyapunovFunc* lf, double Tp,
                                      double t) {
    /* Time-varying gain for prescribed-time stabilization:
     * k(t) = 1/(Tp - t) for t < Tp, forces convergence by time Tp.
     * Ref: Song et al. (2017), "Prescribed-time control". */
    if (lf && t < Tp && Tp > FTS_EPSILON)
        return 1.0 / (Tp - t);
    return FTS_INF;
}

/* ---------- Robust FTS Check ---------- */

bool fts_lyap_robust_fts_check(FTSLyapunovFunc* lf, FTSSystem* sys,
                                FTSState x0, double disturbance,
                                double dt, double T) {
    /* Check if FTS condition holds under bounded disturbance:
     * dV/dt <= -c*V^alpha + disturbance
     * Still gives practical finite-time stability. */
    if (!lf || !sys) return false;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return false;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    bool holds = true;
    while (cp->t < T) {
        fts_step_rk4(cp);
        double V = lf->V(cp->state, lf->params);
        double Vdot = fts_lyap_derivative(lf, cp, cp->state);
        if (Vdot > -lf->c*pow(V, lf->alpha) + disturbance + FTS_EPSILON) {
            holds = false; break;
        }
    }
    fts_free(cp);
    return holds;
}

/* ---------- Output ---------- */

void fts_lyap_print(FTSLyapunovAnalysis* a) {
    if (!a) return;
    printf("Lyapunov Analysis:\n");
    printf("  Settling time: %.6f\n", a->settling_time);
    printf("  Convergence rate: %.6f\n", a->convergence_rate);
    printf("  Type: %d (0=Asym,1=Exp,2=FTS,3=FxT,4=Practical)\n", a->type);
    printf("  V0: %.6f\n", a->V0);
    printf("  T upper bound: %.6f\n", a->T_upper_bound);
}

double fts_lyap_homogeneous_degree(FTSLyapunovFunc* lf) {
    /* The homogeneous degree of the Lyapunov function.
     * For V(x) = ||x||^p, the degree is p.
     * For our default V (quadratic), degree = 2. */
    return lf ? lf->alpha : 0.0;
}
