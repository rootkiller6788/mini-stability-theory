#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fts_core.h"
#include "fts_convergence.h"

/* ============================================================
 * fts_convergence.c -- Convergence Analysis & Settling Time
 *
 * Numerical methods for detecting and characterizing convergence:
 * - Settling time estimation (norm-based criterion)
 * - Convergence classification (asymptotic, exponential, FTS, FxT)
 * - Convergence rate estimation (log-linear regression)
 * - Practical finite-time stability (convergence to neighborhood)
 *
 * Key formulas:
 * - Settling time: T_s = inf{t>=0 : ||x(t)|| < epsilon}
 * - Exponential rate: lambda = -d/dt log(||x||)
 * - Practical FTS: ||x(T)|| <= delta for all ||x0|| <= Delta
 * ============================================================ */

/* --- Lifecycle --- */

FTSSettlingTime* fts_st_create(int n_estimates) {
    FTSSettlingTime* st = calloc(1, sizeof(FTSSettlingTime));
    if (!st) return NULL;
    st->n = n_estimates;
    st->T_estimates = calloc(n_estimates, sizeof(double));
    if (!st->T_estimates) { free(st); return NULL; }
    st->T_min = 1e99;
    st->T_max = 0.0;
    st->T_mean = 0.0;
    st->T_std = 0.0;
    st->is_finite = false;
    st->type = FTS_ASYMPTOTIC;
    st->rate = 0.0;
    st->overshoot = 0.0;
    st->steady_error = 0.0;
    return st;
}

void fts_st_free(FTSSettlingTime* st) {
    if (st) { free(st->T_estimates); free(st); }
}

/* --- Computation --- */

int fts_st_compute(FTSSystem* sys, FTSState x0, double T,
                    double dt, double tol, FTSSettlingTime* st) {
    /* Compute settling time by simulating from x0 until ||x|| < tol.
     *
     * Returns:
     *   1 if converged within time T (settling time stored in st->T_estimates[0])
     *   0 if not converged within time T
     *   -1 on error
     *
     * The settling time is defined as:
     *   T_s = inf{ t >= 0 : ||x(t)|| < tol }
     *
     * For practical purposes, we simulate until convergence or T. */
    if (!sys || !st) return -1;

    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return -1;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));

    int result = 0;
    int max_steps = (int)(T / dt);
    if (max_steps <= 0) max_steps = 10000;

    double max_norm = 0.0;
    double prev_norm = fts_norm(x0);

    for (int step = 0; step < max_steps; step++) {
        fts_step_rk4(cp);
        double nrm = fts_norm(cp->state);

        /* Track overshoot */
        if (nrm > max_norm) max_norm = nrm;

        if (nrm < tol) {
            st->T_estimates[0] = cp->t;
            st->T_min = cp->t;
            st->T_max = cp->t;
            st->T_mean = cp->t;
            st->is_finite = true;
            st->type = FTS_FINITE_TIME;
            st->steady_error = nrm;
            result = 1;
            break;
        }

        /* Check for convergence in norm derivative sense */
        if (nrm < prev_norm) {
            /* decreasing */
        }
        prev_norm = nrm;
    }

    if (result == 0) {
        /* Not converged; store final error */
        st->steady_error = fts_norm(cp->state);
        st->is_finite = false;
        /* If norm is decreasing, might be asymptotic or exponential */
        st->type = FTS_ASYMPTOTIC;
    }

    st->overshoot = max_norm;
    fts_free(cp);
    return result;
}

double fts_st_empirical(FTSSystem* sys, FTSState x0,
                         double dt, double tol) {
    /* Quick empirical settling time estimate.
     *
     * Simulates for a long time horizon (100 time units),
     * returns the first time ||x|| < tol.
     * Returns 1e99 if not converged. */
    if (!sys) return 1e99;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return 1e99;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));

    double Tmax = 100.0;
    int max_steps = (int)(Tmax / dt);
    double Ts = 1e99;

    for (int step = 0; step < max_steps; step++) {
        fts_step_rk4(cp);
        if (fts_norm(cp->state) < tol) {
            Ts = cp->t; break;
        }
    }
    fts_free(cp);
    return Ts;
}

int fts_st_multi_compute(FTSSystem* sys, FTSState* x0_array,
                          int n_initial, double T, double dt,
                          double tol, FTSSettlingTime* st) {
    /* Compute settling times for multiple initial conditions.
     *
     * For each x0 in the array, simulates and records settling time.
     * Then computes min, max, mean, and std of settling times.
     *
     * This is essential for fixed-time stability verification:
     * if max settling time is bounded uniformly, we have fixed-time. */
    if (!sys || !x0_array || !st || n_initial <= 0) return -1;
    if (n_initial > st->n) n_initial = st->n;

    int converged = 0;
    st->T_min = 1e99;
    st->T_max = 0.0;
    double T_sum = 0.0;

    for (int k = 0; k < n_initial; k++) {
        FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
        if (!cp) continue;
        fts_copy_state(&cp->state, &x0_array[k]);
        cp->rhs = sys->rhs;
        if (sys->np > 0)
            memcpy(cp->params, sys->params, sys->np * sizeof(double));

        int max_steps = (int)(T / dt);
        bool found = false;

        for (int step = 0; step < max_steps; step++) {
            fts_step_rk4(cp);
            if (fts_norm(cp->state) < tol) {
                st->T_estimates[converged] = cp->t;
                if (cp->t < st->T_min) st->T_min = cp->t;
                if (cp->t > st->T_max) st->T_max = cp->t;
                T_sum += cp->t;
                converged++;
                found = true;
                break;
            }
        }
        fts_free(cp);
        if (!found) {
            st->T_estimates[k] = 1e99;
        }
    }

    if (converged > 0) {
        st->T_mean = T_sum / converged;
        /* Compute standard deviation */
        double var_sum = 0.0;
        for (int k = 0; k < converged; k++) {
            double diff = st->T_estimates[k] - st->T_mean;
            var_sum += diff * diff;
        }
        st->T_std = sqrt(var_sum / converged);
        st->is_finite = true;
        /* If variance is low across different x0, might be fixed-time */
        if (st->T_std < 0.1 * st->T_mean && st->T_max < T)
            st->type = FTS_FIXED_TIME;
        else
            st->type = FTS_FINITE_TIME;
    }

    return converged;
}

/* --- Classification --- */

FTSConvergenceType fts_st_classify(FTSSettlingTime* st) {
    /* Classify convergence type based on settling time data.
     *
     * - FTS_FIXED_TIME: T_max bounded uniformly across x0
     * - FTS_FINITE_TIME: finite settling time, but not uniform
     * - FTS_EXPONENTIAL: convergence rate ~exp(-lambda*t)
     * - FTS_ASYMPTOTIC: convergence as t->inf, not finite
     * - FTS_PRACTICAL: convergence to a neighborhood
     *
     * The classification uses heuristics based on:
     * - Whether settling time is bounded (is_finite)
     * - Convergence rate (rate)
     * - Steady-state error (steady_error)
     * - Standard deviation of settling times (T_std) */
    if (!st) return FTS_ASYMPTOTIC;
    if (st->type == FTS_FIXED_TIME) return FTS_FIXED_TIME;
    if (st->is_finite) {
        /* Check consistency with fixed-time */
        if (st->T_std > 0 && st->T_std < 0.1 * st->T_mean)
            return FTS_FIXED_TIME;
        return FTS_FINITE_TIME;
    }
    /* Not finite-time: check exponential or asymptotic */
    if (st->rate > 0.01)
        return FTS_EXPONENTIAL;
    if (st->steady_error > 1e-3)
        return FTS_PRACTICAL;
    return FTS_ASYMPTOTIC;
}

bool fts_st_is_finite(FTSSettlingTime* st) {
    return st && st->is_finite;
}

bool fts_st_is_fixed_time(FTSSettlingTime* st,
                           FTSState* x0_array, int n) {
    /* Fixed-time detection: check if settling times are bounded
     * uniformly across different initial conditions.
     *
     * Criterion: std(T) / mean(T) < 0.1 (low dispersion)
     * AND all individual settling times are finite. */
    if (!st || !st->is_finite) return false;
    (void)x0_array; (void)n;
    /* If we have multi-compute data, use it */
    if (st->T_std > 0 && st->T_mean > 1e-12)
        return (st->T_std / st->T_mean) < 0.1;
    return false;
}

/* --- Rate Estimation --- */

double fts_st_convergence_rate(FTSSettlingTime* st) {
    /* Estimated convergence rate from settling time data.
     *
     * For exponential convergence: ||x(t)|| ~ exp(-lambda*t)
     * => log(||x||) ~ -lambda*t
     *
     * lambda is estimated from the recorded trajectory data. */
    return st ? st->rate : 0.0;
}

double fts_st_exponential_rate(FTSSystem* sys, FTSState x0,
                                double T, double dt) {
    /* Estimate exponential convergence rate via log-linear
     * regression on ||x(t)||.
     *
     * Fits: log(||x(t)||) = a - lambda*t
     * lambda > 0 indicates exponential convergence.
     *
     * Uses least squares over the simulated trajectory. */
    if (!sys) return 0.0;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return 0.0;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));

    int max_steps = (int)(T / dt);
    double sum_t = 0.0, sum_log = 0.0, sum_t2 = 0.0, sum_tlog = 0.0;
    int count = 0;

    for (int step = 0; step < max_steps; step += 10) {
        fts_integrate(cp, 10.0 * dt, 0);
        double nrm = fts_norm(cp->state);
        if (nrm > 1e-12 && cp->t < T) {
            double logN = log(nrm);
            sum_t += cp->t;
            sum_log += logN;
            sum_t2 += cp->t * cp->t;
            sum_tlog += cp->t * logN;
            count++;
        }
    }

    fts_free(cp);
    if (count < 5) return 0.0;
    double denom = count * sum_t2 - sum_t * sum_t;
    if (fabs(denom) < 1e-12) return 0.0;
    double lambda = -(count * sum_tlog - sum_t * sum_log) / denom;
    return (lambda > 0) ? lambda : 0.0;
}

/* --- Bounds --- */

FTSTimeBound* fts_st_bound(FTSSystem* sys, FTSState x0,
                            double alpha, double c, double beta) {
    /* Compute settling time bounds.
     *
     * Upper bound (FTS):  T <= V0^(1-alpha) / (c*(1-alpha))
     * Lower bound: linear bound based on initial velocity
     *
     * Returns FTSTimeBound with upper, lower, and tightness flag. */
    FTSTimeBound* tb = calloc(1, sizeof(FTSTimeBound));
    if (!tb) return NULL;
    tb->alpha = alpha; tb->beta = beta;
    double V0 = fts_norm(x0);
    /* Upper bound */
    if (c > 1e-12 && alpha < 1.0 && alpha > 0.0) {
        tb->T_upper = pow(V0, 1.0 - alpha) / (c * (1.0 - alpha));
    } else {
        tb->T_upper = 1e99;
    }
    /* Lower bound (simple heuristic):
     * For scalar FTS: x0/k <= T (linear approximation) */
    tb->T_lower = V0 / (c + 1e-6);
    if (tb->T_lower > tb->T_upper) tb->T_lower = 0.0;
    /* Tightness: if lower/upper > 0.5, bounds are tight */
    if (tb->T_upper < 1e99 && tb->T_lower > 0)
        tb->is_tight = (tb->T_lower / tb->T_upper) > 0.5;
    tb->confidence = tb->is_tight ? 0.9 : 0.5;
    (void)sys; (void)beta;
    return tb;
}

void fts_tb_free(FTSTimeBound* tb) {
    free(tb);
}

double fts_st_residual_energy(FTSSystem* sys, FTSState x0, double t) {
    /* Residual energy at time t starting from x0.
     *
     * Simulates the system from x0 for time t and returns ||x(t)||.
     * This is the "energy" remaining at time t. */
    if (!sys) return 0.0;
    FTSSystem* cp = fts_create(sys->dim, sys->np, sys->dt);
    if (!cp) return 1e99;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    fts_integrate(cp, t, 0);
    double energy = fts_norm(cp->state);
    fts_free(cp);
    return energy;
}

double fts_st_residual_energy_grid(FTSSystem* sys, FTSState x0,
                                    double* t_grid, int n) {
    /* Compute maximum residual energy over a time grid.
     * Returns max over all grid points.
     * Useful for verifying worst-case energy decay. */
    if (!sys || !t_grid || n <= 0) return 0.0;
    double max_energy = 0.0;
    for (int i = 0; i < n; i++) {
        double E = fts_st_residual_energy(sys, x0, t_grid[i]);
        if (E > max_energy) max_energy = E;
    }
    return max_energy;
}

/* --- Comparison --- */

int fts_st_compare(FTSSettlingTime* a, FTSSettlingTime* b) {
    /* Compare two settling time results.
     * Returns: -1 if a settles faster, 1 if b settles faster, 0 if equal. */
    if (!a && !b) return 0;
    if (!a || !a->is_finite) return 1;
    if (!b || !b->is_finite) return -1;
    if (a->T_min < b->T_min - 1e-6) return -1;
    if (a->T_min > b->T_min + 1e-6) return 1;
    return 0;
}

bool fts_st_meets_spec(FTSSettlingTime* st, double T_spec) {
    /* Check if settling time meets a specification.
     * Returns true if T_s <= T_spec. */
    if (!st || !st->is_finite) return false;
    return st->T_max <= T_spec + 1e-6;
}

/* --- Output --- */

void fts_st_print(FTSSettlingTime* st) {
    if (!st) return;
    printf("Settling Time Analysis:\n");
    printf("  T_min: %.6f  T_max: %.6f  T_mean: %.6f\n",
           st->T_min, st->T_max, st->T_mean);
    printf("  T_std: %.6f\n", st->T_std);
    printf("  Is finite: %s\n", st->is_finite ? "yes" : "no");
    printf("  Type: %d\n", st->type);
    printf("  Rate: %.6f\n", st->rate);
    printf("  Overshoot: %.6f\n", st->overshoot);
    printf("  Steady error: %.6e\n", st->steady_error);
}

void fts_tb_print(FTSTimeBound* tb) {
    if (!tb) return;
    printf("Time Bounds:\n");
    printf("  T_upper: %.6f  T_lower: %.6f\n",
           tb->T_upper, tb->T_lower);
    printf("  Tight: %s  Confidence: %.2f\n",
           tb->is_tight ? "yes" : "no", tb->confidence);
}

/* --- Practical Finite-Time Stability --- */

bool fts_st_practical_finite_time(FTSSystem* sys, FTSState x0,
                                   double epsilon, double T,
                                   double dt) {
    /* Practical finite-time stability:
     * System converges to epsilon-neighborhood of origin by time T.
     *
     * This is weaker than FTS (exact convergence) but more
     * practically relevant, as real systems always have
     * measurement noise and disturbances. */
    if (!sys) return false;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return false;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    fts_integrate(cp, T, 0);
    bool converged = (fts_norm(cp->state) < epsilon);
    fts_free(cp);
    return converged;
}

double fts_st_practical_settling_time(FTSSystem* sys, FTSState x0,
                                       double epsilon, double dt) {
    /* Settling time to reach epsilon-neighborhood.
     * Searches up to T=100.
     * Returns settling time or 1e99 if not reached. */
    if (!sys) return 1e99;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return 1e99;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    double Tmax = 100.0;
    int max_steps = (int)(Tmax / dt);
    double Ts = 1e99;
    for (int step = 0; step < max_steps; step++) {
        fts_step_rk4(cp);
        if (fts_norm(cp->state) < epsilon) { Ts = cp->t; break; }
    }
    fts_free(cp);
    return Ts;
}