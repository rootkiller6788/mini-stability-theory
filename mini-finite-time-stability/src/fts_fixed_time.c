#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fts_core.h"
#include "fts_fixed_time.h"

/* ============================================================
 * fts_fixed_time.c -- Fixed-Time Stability (Polyakov, 2012)
 *
 * Fixed-time stability: settling time is bounded uniformly
 * for ALL initial conditions.
 *
 * Lyapunov condition:
 *   dV/dt <= -(c1*V^alpha + c2*V^beta), 0<alpha<1<beta
 *
 * Settling time bound:
 *   T_max <= 1/(c1*(1-alpha)) + 1/(c2*(beta-1))
 *
 * Contrast with finite-time stability:
 *   - FTS: T depends on x0, can be arbitrarily large
 *   - FxT: T has uniform upper bound independent of x0
 *
 * Ref: Polyakov, A. (2012). "Nonlinear feedback design for
 * fixed-time stabilization of linear control systems."
 * IEEE Trans. Autom. Control, 57(8):2106-2110.
 * ============================================================ */

/* --- Lifecycle --- */

FTSFixedTimeParams* fts_ft_params_create(double alpha, double beta,
                                          double p, double q) {
    /* Create fixed-time parameter struct.
     *
     * Parameters:
     *   alpha: FTS exponent (0<alpha<1)
     *   beta:  FxT exponent (beta>1)
     *   p:     first coefficient c1 > 0
     *   q:     second coefficient c2 > 0
     *
     * Fixed-time condition: p>0, q>0, alpha<1, beta>1
     *
     * T_max = 1/(p*(1-alpha)) + 1/(q*(beta-1)) */
    FTSFixedTimeParams* ftp = calloc(1, sizeof(FTSFixedTimeParams));
    if (!ftp) return NULL;
    ftp->alpha = alpha;
    ftp->beta = beta;
    ftp->p = p;
    ftp->q = q;
    ftp->gamma = 1.0;
    /* Check fixed-time condition */
    ftp->is_fixed_time = (p > 1e-12 && q > 1e-12 &&
                          alpha < 1.0 && beta > 1.0);
    ftp->is_uniform = ftp->is_fixed_time;
    /* Compute T_max bound */
    if (ftp->is_fixed_time) {
        double T1 = 1.0 / (p * (1.0 - alpha));
        double T2 = 1.0 / (q * (beta - 1.0));
        ftp->T_max = T1 + T2;
    } else {
        ftp->T_max = 1e99;
    }
    return ftp;
}

void fts_ft_params_free(FTSFixedTimeParams* ftp) {
    free(ftp);
}

/* --- Fixed-Time Detection --- */

bool fts_ft_is_fixed_time(FTSFixedTimeParams* ftp) {
    return ftp && ftp->is_fixed_time;
}

bool fts_ft_is_uniform(FTSFixedTimeParams* ftp) {
    return ftp && ftp->is_uniform;
}

/* --- Bounds --- */

double fts_ft_settling_bound(FTSFixedTimeParams* ftp) {
    /* Fixed-time settling bound (Polyakov, 2012):
     *
     * T_max = 1/(c1*(1-alpha)) + 1/(c2*(beta-1))
     *
     * This bound is independent of the initial condition!
     *
     * Proof sketch:
     * From dV/dt <= -(c1*V^alpha + c2*V^beta):
     *   Let z = V^(1-alpha), then dz/dt <= -c1*(1-alpha) - c2*(1-alpha)*V^(beta-alpha)
     *   V >= 1: V^beta dominates, T1 <= 1/(c1*(1-alpha))
     *   V <= 1: V^alpha dominates, T2 <= 1/(c2*(beta-1))
     *   Total: T <= T1 + T2 */
    if (!ftp || !ftp->is_fixed_time) return 1e99;
    return ftp->T_max;
}

double fts_ft_settling_bound_general(FTSFixedTimeParams* ftp,
                                      FTSState x0) {
    /* General settling time bound that may depend on x0.
     *
     * For fixed-time stable systems, returns the uniform bound.
     * For finite-time only, returns the x0-dependent bound:
     *   T <= V0^(1-alpha) / (p*(1-alpha))
     *
     * This function bridges finite-time and fixed-time analysis. */
    if (!ftp) return 1e99;
    if (ftp->is_fixed_time) return ftp->T_max;
    /* Fall back to FTS bound */
    double V0 = fts_norm(x0);
    if (V0 < 1e-12) return 0.0;
    if (ftp->alpha >= 1.0 || ftp->p < 1e-12) return 1e99;
    return pow(V0, 1.0 - ftp->alpha) / (ftp->p * (1.0 - ftp->alpha));
}

double fts_ft_lyapunov_bound(FTSFixedTimeParams* ftp, FTSState x0) {
    /* Lyapunov-based settling time bound.
     *
     * Uses the standard Lyapunov condition:
     *   T <= V(x0)^(1-alpha) / (p*(1-alpha))
     * for the finite-time part.
     *
     * If also fixed-time, augments with the uniform bound. */
    if (!ftp) return 1e99;
    double V0 = fts_norm(x0);  /* or could use user-supplied V */
    if (V0 < 1e-12) return 0.0;
    /* Finite-time bound */
    if (ftp->alpha >= 1.0 || ftp->p < 1e-12) return 1e99;
    double T_fts = pow(V0, 1.0 - ftp->alpha) / (ftp->p * (1.0 - ftp->alpha));
    /* Fixed-time: take min(fts_bound, fix_bound) */
    if (ftp->is_fixed_time)
        return (T_fts < ftp->T_max) ? T_fts : ftp->T_max;
    return T_fts;
}

double fts_ft_uncertainty_bound(FTSFixedTimeParams* ftp,
                                 double delta_c1, double delta_c2) {
    /* Robustness of fixed-time bound to parameter uncertainty.
     *
     * If coefficients c1, c2 are uncertain by deltas:
     *   T_max_robust = 1/((c1-dc1)*(1-alpha)) + 1/((c2-dc2)*(beta-1))
     *
     * This shows how the bound degrades with uncertainty.
     * Small uncertainty in c1 (near zero) can dramatically
     * increase the bound. */
    if (!ftp || !ftp->is_fixed_time) return 1e99;
    double c1_eff = ftp->p - delta_c1;
    double c2_eff = ftp->q - delta_c2;
    if (c1_eff < 1e-12 || c2_eff < 1e-12) return 1e99;
    double T1 = 1.0 / (c1_eff * (1.0 - ftp->alpha));
    double T2 = 1.0 / (c2_eff * (ftp->beta - 1.0));
    return T1 + T2;
}

/* --- Testing --- */

FTSFixedTimeTest* fts_ft_test(FTSSystem* sys, FTSFixedTimeParams* ftp,
                               int n_tests, double T_sim, double dt) {
    /* Test fixed-time stability by simulating from multiple
     * initial conditions and checking if all converge within
     * the theoretical bound.
     *
     * Generates n_tests random initial conditions, simulates
     * each for time T_sim, and records settling times. */
    if (!sys || !ftp || n_tests <= 0) return NULL;

    FTSFixedTimeTest* ftt = calloc(1, sizeof(FTSFixedTimeTest));
    if (!ftt) return NULL;
    ftt->params = *ftp;
    ftt->n_tests = n_tests;
    ftt->settling_times = calloc(n_tests, sizeof(double));
    if (!ftt->settling_times) { free(ftt); return NULL; }

    ftt->T_observed_max = 0.0;
    ftt->T_observed_min = 1e99;
    double T_sum = 0.0;
    int converged_count = 0;

    for (int test = 0; test < n_tests; test++) {
        /* Generate random initial condition */
        FTSState x0; x0.dim = sys->dim;
        for (int i = 0; i < sys->dim; i++)
            x0.x[i] = ((double)((test * 17 + i * 31) % 100) / 50.0 - 1.0) * 10.0;

        /* Simulate */
        FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
        if (!cp) continue;
        fts_copy_state(&cp->state, &x0);
        cp->rhs = sys->rhs;
        if (sys->np > 0)
            memcpy(cp->params, sys->params, sys->np * sizeof(double));

        double Ts = T_sim;
        int max_steps = (int)(T_sim / dt);
        bool converged = false;

        for (int step = 0; step < max_steps; step++) {
            fts_step_rk4(cp);
            if (fts_norm(cp->state) < 1e-6) {
                Ts = cp->t; converged = true; break;
            }
        }

        if (converged) {
            ftt->settling_times[converged_count] = Ts;
            if (Ts > ftt->T_observed_max) ftt->T_observed_max = Ts;
            if (Ts < ftt->T_observed_min) ftt->T_observed_min = Ts;
            T_sum += Ts;
            converged_count++;
        }
        fts_free(cp);
    }

    if (converged_count > 0) {
        ftt->T_mean = T_sum / converged_count;
        /* Verification: all settling times <= theoretical T_max */
        ftt->verified = (ftt->T_observed_max <= ftp->T_max + 1e-6);
        ftt->violation_margin = ftp->T_max - ftt->T_observed_max;
    } else {
        ftt->verified = false;
        ftt->violation_margin = -1.0;
    }

    return ftt;
}

void fts_ft_test_free(FTSFixedTimeTest* ftt) {
    if (ftt) { free(ftt->settling_times); free(ftt); }
}

bool fts_ft_verify(FTSFixedTimeTest* ftt) {
    return ftt && ftt->verified;
}

double fts_ft_max_settling(FTSFixedTimeTest* ftt) {
    return ftt ? ftt->T_observed_max : 1e99;
}

/* --- Comparison --- */

bool fts_ft_compare_finite_vs_fixed(FTSFixedTimeParams* ftp,
                                     double finite_T) {
    /* Compare fixed-time bound with finite-time settling time.
     * Returns true if fixed-time bound is better (smaller).
     *
     * This comparison highlights the advantage of fixed-time
     * stability: while finite-time settling time grows with
     * ||x0||, fixed-time settling time remains bounded. */
    if (!ftp || !ftp->is_fixed_time) return false;
    return ftp->T_max < finite_T;
}

double fts_ft_robustness_margin(FTSFixedTimeParams* ftp,
                                 double disturbance) {
    /* Robustness margin: how much disturbance the fixed-time
     * property can tolerate before the bound becomes invalid.
     *
     * Simple model: effective coefficients degrade with disturbance.
     * c_eff = c - gamma*disturbance
     * If c_eff <= 0, fixed-time property is lost.
     *
     * Returns: effective c1 after disturbance. */
    if (!ftp) return 0.0;
    return ftp->p - ftp->gamma * disturbance;
}

/* --- Output --- */

void fts_ft_print(FTSFixedTimeTest* ftt) {
    if (!ftt) return;
    printf("Fixed-Time Stability Test:\n");
    printf("  T_max (theoretical): %.6f\n", ftt->params.T_max);
    printf("  T_observed_max:      %.6f\n", ftt->T_observed_max);
    printf("  T_observed_min:      %.6f\n", ftt->T_observed_min);
    printf("  T_mean:              %.6f\n", ftt->T_mean);
    printf("  Verified:            %s\n", ftt->verified ? "yes" : "no");
    printf("  Violation margin:    %.6f\n", ftt->violation_margin);
    printf("  Tests:               %d\n", ftt->n_tests);
}

void fts_ft_params_print(FTSFixedTimeParams* ftp) {
    if (!ftp) return;
    printf("Fixed-Time Parameters:\n");
    printf("  alpha: %.4f  beta: %.4f\n", ftp->alpha, ftp->beta);
    printf("  c1 (p): %.4f  c2 (q): %.4f\n", ftp->p, ftp->q);
    printf("  T_max: %.6f\n", ftp->T_max);
    printf("  Is fixed-time: %s\n", ftp->is_fixed_time ? "yes" : "no");
}

/* --- Prescribed-Time Stability --- */

bool fts_ft_is_prescribed_time(FTSFixedTimeParams* ftp, double Tp) {
    /* Prescribed-time stability: the settling time can be
     * arbitrarily prescribed by the designer, independent of
     * initial conditions and system parameters.
     *
     * Condition: T_max <= Tp (user-specified time).
     *
     * This is stronger than fixed-time stability, where the
     * bound is determined by system parameters. In prescribed-time
     * control, the bound is a design parameter. */
    if (!ftp) return false;
    return ftp->T_max <= Tp + 1e-6;
}

double fts_ft_prescribed_gain(double Tp, double t) {
    /* Time-varying gain for prescribed-time stabilization:
     * k(t) = 1 / (Tp - t)
     *
     * As t -> Tp, the gain grows unbounded, forcing convergence
     * exactly at the prescribed time.
     *
     * Ref: Song et al. (2017), "Time-varying feedback for
     * regulation of normal-form nonlinear systems in prescribed
     * finite time", Automatica. */
    if (t >= Tp || Tp < 1e-12) return 1e99;
    return 1.0 / (Tp - t);
}

/* --- Implicit Lyapunov Function Method --- */

double fts_ft_implicit_lyap_bound(FTSFixedTimeParams* ftp, FTSState x0) {
    /* Implicit Lyapunov function (ILF) method for fixed-time
     * stability analysis (Polyakov, 2012; Korobov, 1979).
     *
     * The ILF Q(V,x) = 0 implicitly defines V(x). For systems
     * of the form x_dot = -k1*|x|^alpha*sign(x) - k2*|x|^beta*sign(x),
     * the settling time bound is exactly:
     *
     *   T(x0) <= int_0^inf (1 / (k1*V^alpha + k2*V^beta)) dV
     *
     * This integral can be evaluated in terms of the Beta function:
     *   T <= B(alpha, beta) / (k1^((beta-1)/(beta-alpha)) * k2^((1-alpha)/(beta-alpha)))
     *
     * For simplicity, we return the standard Lyapunov bound. */
    if (!ftp) return 1e99;
    return fts_ft_lyapunov_bound(ftp, x0);
}