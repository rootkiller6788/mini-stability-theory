#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "popov_criterion.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * Popov Data Management
 *
 * Popov plot maps: x = Re[G(jw)], y = w * Im[G(jw)]
 * This is a transformed version of the Nyquist plot that enables
 * a simple graphical test for absolute stability.
 * ================================================================ */

PopovData* popov_data_create(int n) {
    PopovData* pd = (PopovData*)calloc(1, sizeof(PopovData));
    if (!pd) return NULL;
    pd->n = n;
    pd->omega_re = (double*)calloc((size_t)n, sizeof(double));
    pd->omega_im = (double*)calloc((size_t)n, sizeof(double));
    if (!pd->omega_re || !pd->omega_im) {
        free(pd->omega_re); free(pd->omega_im); free(pd); return NULL;
    }
    return pd;
}

void popov_data_free(PopovData* pd) {
    if (!pd) return;
    free(pd->omega_re);
    free(pd->omega_im);
    free(pd);
}

/* ================================================================
 * Popov Plot Computation
 *
 * For each frequency w, compute:
 *   x = Re[G(jw)]
 *   y = w * Im[G(jw)]
 *
 * The Popov plot emphasizes the high-frequency behavior through
 * the w-multiplication of the imaginary part.
 * ================================================================ */

void popov_compute(const TransferFunction* tf, PopovData* pd,
                   double wmin, double wmax) {
    if (!tf || !pd || pd->n < 2) return;
    double dw = (wmax - wmin) / (double)(pd->n - 1);
    for (int i = 0; i < pd->n; i++) {
        double w = wmin + (double)i * dw;
        double re, im;
        tf_eval_complex(tf, w, &re, &im);
        pd->omega_re[i] = re;
        pd->omega_im[i] = w * im;
    }
}

/* Log-spaced Popov plot computation */
void popov_compute_log(const TransferFunction* tf, PopovData* pd,
                        double wmin, double wmax) {
    if (!tf || !pd || pd->n < 2) return;
    if (wmin < 1e-8) wmin = 1e-8;
    double lmin = log10(wmin), lmax = log10(wmax);
    double dl = (lmax - lmin) / (double)(pd->n - 1);
    for (int i = 0; i < pd->n; i++) {
        double w = pow(10.0, lmin + (double)i * dl);
        double re, im;
        tf_eval_complex(tf, w, &re, &im);
        pd->omega_re[i] = re;
        pd->omega_im[i] = w * im;
    }
}

/* ================================================================
 * Sector Loop Transformation
 *
 * Given G(s) with phi in sector [alpha, beta], apply the loop
 * transformation to obtain an equivalent system with phi' in [0, k]
 * where k = beta - alpha.
 *
 * Transformed linear part:
 *   G'(s) = G(s) / (1 + alpha * G(s))
 *
 * This is derived from:
 *   y   = G(s) * u
 *   u   = -phi(y)
 *   Let v = u + alpha*y  =>  y = G(s)/(1+alpha*G(s)) * v
 *   Then v = -phi'(y) where phi' in [0, beta-alpha]
 * ================================================================ */

TransferFunction* popov_loop_transform(const TransferFunction* tf,
                                        SectorBounds sb, double *k) {
    *k = sb.beta - sb.alpha;
    if (!tf || *k <= 0.0) return NULL;

    /* For small alpha, G'(s) ~ G(s). We return a copy.
     * For general alpha, we need polynomial algebra:
     * G'(s) = num(s) / (den(s) + alpha*num(s))
     *
     * New numerator = num(s)
     * New denominator = den(s) + alpha*num(s)
     * The degree is max(n_den, n_num) */

    int n_num = tf->n_num;
    int n_den_new = tf->n_den;
    if (n_num > n_den_new) n_den_new = n_num;

    TransferFunction* tf_new = tf_create(n_num, n_den_new);
    if (!tf_new) return NULL;

    /* Copy numerator */
    for (int i = 0; i < n_num; i++)
        tf_new->num[i] = tf->num[i];

    /* New denominator = den + alpha * num */
    for (int i = 0; i < tf->n_den; i++)
        tf_new->den[i] = tf->den[i];
    for (int i = 0; i < n_num; i++)
        tf_new->den[i] += sb.alpha * tf->num[i];

    return tf_new;
}

/* ================================================================
 * Popov Line Test
 *
 * The Popov criterion: there exists eta >= 0 such that
 *   Re[(1 + j*w*eta) * G(jw)] + 1/k > 0  for all w >= 0
 *
 * In the Popov plot (x = Re[G], y = w*Im[G]):
 *   y < (1/eta) * (x + 1/k)
 *
 * i.e., the Popov plot lies strictly to the RIGHT of the line
 * passing through (-1/k, 0) with slope 1/eta.
 *
 * This function tests the line condition:
 *   y < slope * x + intercept
 *
 * For the Popov condition, we need:
 *   slope = 1/eta, intercept = 1/(eta*k)
 *   or equivalently: y < slope * (x + 1/k)
 *
 * Returns true if all points satisfy the condition.
 * *violation = max(y - (slope*x + intercept)), >0 means violation.
 * ================================================================ */

bool popov_line_test(const PopovData* pd, double slope, double intercept,
                     double *violation) {
    *violation = 0.0;
    if (!pd || pd->n < 2) return false;

    double max_violation = -1e30;
    for (int i = 0; i < pd->n; i++) {
        double line_y = slope * pd->omega_re[i] + intercept;
        double v = pd->omega_im[i] - line_y;
        if (v > max_violation) max_violation = v;
    }

    *violation = max_violation;
    return max_violation <= 1e-9;
}

/* ================================================================
 * Popov Criterion Check (with specified eta)
 *
 * Algorithm:
 *   1. Apply loop transformation to shift sector to [0, k].
 *   2. Generate Popov plot of the transformed G'(s).
 *   3. Test the Popov line condition with specified eta.
 *   4. Return margin = line_distance (positive => stable).
 * ================================================================ */

int popov_criterion_check(const TransferFunction* tf, SectorBounds sb,
                           double eta, double *margin) {
    *margin = 0.0;
    if (!tf || eta < 0.0) return 0;

    /* Loop-transform to canonical sector [0, k] */
    double k;
    TransferFunction* tf_trans = popov_loop_transform(tf, sb, &k);
    if (!tf_trans) {
        /* If alpha=0, G' = G directly */
        tf_trans = tf_create(tf->n_num, tf->n_den);
        if (!tf_trans) return 0;
        for (int i = 0; i < tf->n_num; i++) tf_trans->num[i] = tf->num[i];
        for (int i = 0; i < tf->n_den; i++) tf_trans->den[i] = tf->den[i];
        k = sb.beta;
    }

    if (k <= 0.0) { tf_free(tf_trans); return 0; }

    /* Generate Popov plot */
    PopovData* pd = popov_data_create(POPOV_N_POINTS);
    if (!pd) { tf_free(tf_trans); return 0; }
    popov_compute_log(tf_trans, pd, 0.001, POPOV_MAX_FREQ);

    /* Popov line: passes through (-1/k, 0) with slope 1/eta
     * Line equation: y = (1/eta) * (x + 1/k)
     *               = slope * x + slope/k
     * slope = 1/eta, intercept = slope/k = 1/(eta*k) */
    double slope, intercept;
    if (eta > 1e-10) {
        slope = 1.0 / eta;
        intercept = slope / k;  /* = 1/(eta*k) */
    } else {
        /* eta = 0 => vertical line: x > -1/k
         * We use a very large slope to approximate */
        slope = 1e10;
        intercept = slope / k;
    }

    double violation;
    bool pass = popov_line_test(pd, slope, intercept, &violation);

    /* Margin: negative violation means stable (all points below line) */
    *margin = -violation;

    popov_data_free(pd);
    tf_free(tf_trans);
    return pass ? 1 : 0;
}

/* ================================================================
 * Optimal Popov Multiplier Search
 *
 * Searches for the eta >= 0 that maximizes the stability margin.
 * Uses a uniform grid search over [eta_min, eta_max].
 *
 * The Popov multiplier eta exploits the time-invariance of the
 * nonlinearity. Larger eta means more "phase lead" compensation,
 * but there is a limit beyond which the condition becomes
 * conservative again.
 * ================================================================ */

void popov_find_optimal_slope(const TransferFunction* tf,
                               SectorBounds sb,
                               double *eta_opt, double *margin_opt) {
    *eta_opt = 0.0;
    *margin_opt = -1e30;

    for (int i = 0; i < POPOV_ETA_STEPS; i++) {
        double eta = POPOV_ETA_MAX * (double)i / (double)(POPOV_ETA_STEPS - 1);
        double margin;
        if (popov_criterion_check(tf, sb, eta, &margin)) {
            if (margin > *margin_opt) {
                *margin_opt = margin;
                *eta_opt = eta;
            }
        }
    }
}

/* ================================================================
 * Grid Search with Custom Range
 *
 * More flexible search allowing user to specify the eta range.
 * Returns the best eta found. If no feasible eta exists,
 * *best_margin is set to -1e30.
 * ================================================================ */

double popov_grid_search_eta(const TransferFunction* tf,
                              SectorBounds sb,
                              double eta_min, double eta_max,
                              int n_steps, double *best_margin) {
    double best_eta = 0.0;
    *best_margin = -1e30;

    if (n_steps < 2) n_steps = 2;
    double deta = (eta_max - eta_min) / (double)(n_steps - 1);

    for (int i = 0; i < n_steps; i++) {
        double eta = eta_min + (double)i * deta;
        if (eta < 0.0) eta = 0.0;
        double margin;
        if (popov_criterion_check(tf, sb, eta, &margin)) {
            if (margin > *best_margin) {
                *best_margin = margin;
                best_eta = eta;
            }
        }
    }

    return best_eta;
}

/* ================================================================
 * Popov Plot Data Generator (for external use)
 *
 * Fills arrays with the Popov plot data for plotting or
 * further analysis.
 * ================================================================ */

void popov_plot_data(const TransferFunction* tf, SectorBounds sb,
                     double *x, double *y, int n) {
    (void)sb;
    if (!tf || !x || !y || n < 2) return;
    PopovData* pd = popov_data_create(n);
    if (!pd) return;
    popov_compute(tf, pd, 0.01, POPOV_MAX_FREQ);
    for (int i = 0; i < n; i++) {
        x[i] = pd->omega_re[i];
        y[i] = pd->omega_im[i];
    }
    popov_data_free(pd);
}

/* ================================================================
 * Compare Circle vs Popov Criteria
 *
 * Runs both tests on the same system and reports which criterion
 * succeeds. Since Popov exploits time-invariance, it can sometimes
 * certify stability where the (more conservative) circle criterion
 * cannot.
 *
 * Returns bitmask:
 *   bit 0 (value 1) = Circle criterion passes
 *   bit 1 (value 2) = Popov criterion passes
 *   => 0=neither, 1=circle only, 2=popov only, 3=both
 * ================================================================ */

int popov_compare_criteria(const TransferFunction* tf, SectorBounds sb,
                            double *circle_margin, double *popov_margin,
                            double *eta_found) {
    *circle_margin = 0.0;
    *popov_margin  = 0.0;
    *eta_found     = 0.0;

    int result = 0;

    /* Run circle criterion */
    if (circle_criterion_check(tf, sb, circle_margin))
        result |= 1;

    /* Run Popov criterion with optimal eta search */
    double eta_opt, margin_opt;
    popov_find_optimal_slope(tf, sb, &eta_opt, &margin_opt);

    if (margin_opt > 0.0) {
        result |= 2;
        *popov_margin = margin_opt;
        *eta_found = eta_opt;
    } else {
        /* Try with eta=0.0 specifically */
        if (popov_criterion_check(tf, sb, 0.0, popov_margin))
            result |= 2;
    }

    return result;
}

/* ================================================================
 * Self-Test Routine
 * ================================================================ */

int popov_criterion_self_test(void) {
    /* G(s) = 1/(s+1), phi in sector [0, 2] */
    TransferFunction* tf = tf_create(1, 2);
    assert(tf != NULL);
    tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;

    /* Popov data creation */
    PopovData* pd = popov_data_create(100);
    assert(pd != NULL);
    popov_compute(tf, pd, 0.01, 10.0);
    assert(pd->n == 100);
    /* For 1/(s+1): Re[G] = 1/(1+w^2), Im[G] = -w/(1+w^2)
     * At w=0: (Re, w*Im) = (1, 0)
     * At w=1: (Re, w*Im) = (0.5, -0.5) */
    assert(fabs(pd->omega_re[0] - 1.0) < 0.05);
    popov_data_free(pd);

    /* Log-spaced Popov */
    PopovData* pd2 = popov_data_create(200);
    assert(pd2 != NULL);
    popov_compute_log(tf, pd2, 0.01, 100.0);
    assert(pd2->n == 200);
    assert(fabs(pd2->omega_re[0] - 1.0) < 0.1);
    popov_data_free(pd2);

    /* Popov criterion check: sector [0, 2], eta=1.0 */
    SectorBounds sb = {0.0, 2.0, 0, 0, 0.0, 0.0};
    double margin;
    int result = popov_criterion_check(tf, sb, 1.0, &margin);
    assert(result == 1);  /* should be stable */
    assert(margin > 0.0);

    /* Popov criterion with eta=0 should also pass */
    result = popov_criterion_check(tf, sb, 0.0, &margin);
    assert(result == 1);

    /* Find optimal eta */
    double eta_opt, margin_opt;
    popov_find_optimal_slope(tf, sb, &eta_opt, &margin_opt);
    assert(eta_opt >= 0.0);
    assert(margin_opt > 0.0);

    /* Grid search */
    double best_m;
    double best_eta = popov_grid_search_eta(tf, sb, 0.0, 10.0, 20, &best_m);
    assert(best_eta >= 0.0);
    assert(best_m > 0.0);

    /* Popov line test */
    PopovData* pd3 = popov_data_create(50);
    assert(pd3 != NULL);
    popov_compute(tf, pd3, 0.01, 10.0);
    double violation;
    bool line_ok = popov_line_test(pd3, 1.0, 1.0, &violation);
    (void)line_ok; (void)violation;
    popov_data_free(pd3);

    /* Compare criteria */
    SectorBounds sb2 = {0.0, 3.0, 0, 0, 0.0, 0.0};
    double cm, pm, ef;
    int comp = popov_compare_criteria(tf, sb2, &cm, &pm, &ef);
    (void)comp; (void)cm; (void)pm; (void)ef;

    /* Loop transform: sector [0.5, 2.5] => [0, 2.0] */
    SectorBounds sb3 = {0.5, 2.5, 0, 0, 0.0, 0.0};
    double k;
    TransferFunction* tf2 = popov_loop_transform(tf, sb3, &k);
    assert(tf2 != NULL);
    assert(fabs(k - 2.0) < 0.01);
    /* G'(s) = G(s)/(1+0.5*G(s)) = (1/(s+1)) / (1+0.5/(s+1)) = 1/(s+1.5) */
    /* So G'(j0) = 1/1.5 ~ 0.667 */
    double re2, im2;
    tf_eval_complex(tf2, 0.0, &re2, &im2);
    assert(fabs(re2 - 0.6667) < 0.01);
    tf_free(tf2);

    tf_free(tf);
    return 1;
}
