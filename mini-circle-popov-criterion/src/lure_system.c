#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "lure_system.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ===== Lur'e System Lifecycle ===== */
LureSystem* lure_create(StateSpace* ss, SectorBounds sb) {
    LureSystem* ls = (LureSystem*)calloc(1, sizeof(LureSystem));
    if (!ls) return NULL;
    ls->linear = ss;
    ls->nonlinear = sb;
    ls->eta = 0.0;
    ls->config = LURE_CONFIG_STANDARD;
    return ls;
}

LureSystem* lure_deep_copy(const LureSystem* ls) {
    if (!ls) return NULL;
    StateSpace* ss_copy = ss_create(ls->linear->n);
    if (!ss_copy) return NULL;
    int n = ls->linear->n;
    for (int i = 0; i < n * n; i++) ss_copy->A[i] = ls->linear->A[i];
    for (int i = 0; i < n; i++) { ss_copy->B[i] = ls->linear->B[i]; ss_copy->C[i] = ls->linear->C[i]; }
    ss_copy->D = ls->linear->D;
    LureSystem* ls2 = lure_create(ss_copy, ls->nonlinear);
    if (!ls2) { ss_free(ss_copy); return NULL; }
    ls2->eta = ls->eta; ls2->config = ls->config;
    return ls2;
}

void lure_free(LureSystem* ls) {
    if (!ls) return;
    /* Note: does NOT free ls->linear; caller owns it.
     * For deep copies, caller must also free ls->linear separately. */
    free(ls);
}

/* ===== Sector Utilities ===== */
double lure_sector_lower_bound(const LureSystem* ls) {
    return ls ? ls->nonlinear.alpha : 0.0;
}

double lure_sector_upper_bound(const LureSystem* ls) {
    return ls ? ls->nonlinear.beta : 1.0;
}

bool lure_sector_check(const LureSystem* ls, double sigma, double phi) {
    if (!ls) return false;
    double prod = sigma * phi;
    double lo = ls->nonlinear.alpha * sigma * sigma;
    double hi = ls->nonlinear.beta  * sigma * sigma;
    return (prod >= lo - 1e-12) && (prod <= hi + 1e-12);
}

/* ===== Absolute Stability Test (Unified) ===== */
int lure_absolute_stability(const LureSystem* ls, int method,
                             double *margin, double *eta) {
    *margin = 0.0;
    if (eta) *eta = 0.0;
    if (!ls || !ls->linear) return -1;

    switch (method) {
    case LURE_METHOD_CIRCLE:
        return lure_circle_test(ls, margin);
    case LURE_METHOD_POPOV: {
        double e;
        int r = lure_popov_test(ls, margin, &e);
        if (eta) *eta = e;
        return r;
    }
    case LURE_METHOD_KYP:
        return lure_kyp_test(ls, margin);
    case LURE_METHOD_ALL: {
        /* Try circle first (fastest) */
        double cm;
        if (lure_circle_test(ls, &cm)) {
            *margin = cm; if (eta) *eta = 0.0; return 1;
        }
        /* Try Popov */
        double pm, pe;
        if (lure_popov_test(ls, &pm, &pe)) {
            *margin = pm; if (eta) *eta = pe; return 1;
        }
        /* Try KYP */
        double km;
        if (lure_kyp_test(ls, &km)) {
            *margin = km; if (eta) *eta = 0.0; return 1;
        }
        return 0;
    }
    default:
        return -1;
    }
}

/* ===== Circle Criterion Test ===== */
int lure_circle_test(const LureSystem* ls, double *margin) {
    *margin = 0.0;
    if (!ls || !ls->linear) return -1;

    /* Convert state-space to transfer function */
    int n = ls->linear->n;
    TransferFunction* tf = tf_create(n + 1, n + 1);
    if (!tf) return -1;
    kyp_to_transfer_function(ls->linear, tf);

    int result = circle_criterion_check(tf, ls->nonlinear, margin);
    tf_free(tf);
    return result;
}

/* ===== Popov Criterion Test ===== */
int lure_popov_test(const LureSystem* ls, double *margin, double *eta_opt) {
    *margin = 0.0;
    if (eta_opt) *eta_opt = 0.0;
    if (!ls || !ls->linear) return -1;

    int n = ls->linear->n;
    TransferFunction* tf = tf_create(n + 1, n + 1);
    if (!tf) return -1;
    kyp_to_transfer_function(ls->linear, tf);

    /* Find optimal eta */
    double best_eta, best_margin;
    popov_find_optimal_slope(tf, ls->nonlinear, &best_eta, &best_margin);

    if (best_margin > 0.0) {
        *margin = best_margin;
        if (eta_opt) *eta_opt = best_eta;
        tf_free(tf);
        return 1;
    }

    tf_free(tf);
    return 0;
}

/* ===== KYP-Based Test ===== */
int lure_kyp_test(const LureSystem* ls, double *margin) {
    *margin = 0.0;
    if (!ls || !ls->linear) return -1;

    /* Check if the linear part is strictly positive real after
     * loop transformation to [0, k] sector. */
    KYPSolution* sol = kyp_solve(ls->linear, 0.01);
    if (!sol) return -1;

    if (sol->is_feasible) {
        *margin = 1.0;  /* rough margin */
        kyp_free(sol);
        return 1;
    }

    *margin = 0.0;
    kyp_free(sol);
    return 0;
}

/* ===== Gain and Phase Margin Computation =====
 * Computes margins on the LINEARIZED system (phi = gain*sigma).
 * These are classical Bode/Nyquist margins, not absolute stability margins. */
StabilityMargin* lure_gain_margin(const LureSystem* ls) {
    StabilityMargin* sm = (StabilityMargin*)calloc(1, sizeof(StabilityMargin));
    if (!sm) return NULL;
    if (!ls || !ls->linear) { sm->gain_margin = 1e9; sm->phase_margin = 180.0; return sm; }

    int n = ls->linear->n;
    TransferFunction* tf = tf_create(n + 1, n + 1);
    if (!tf) { sm->gain_margin = 1e9; sm->phase_margin = 180.0; return sm; }
    kyp_to_transfer_function(ls->linear, tf);

    /* Find gain crossover frequency (where |G(jw)| = 1) */
    double wc = circle_criterion_crossover_freq(tf);
    sm->frequency = wc;

    /* Phase at crossover */
    double re, im;
    tf_eval_complex(tf, wc, &re, &im);
    double phase = atan2(im, re) * 180.0 / M_PI;
    sm->phase_margin = 180.0 + phase;
    if (sm->phase_margin > 180.0) sm->phase_margin -= 360.0;

    /* Gain margin: find phase crossover (where phase = -180 deg) */
    double w180 = 0.0;
    double min_phase_err = 1e30;
    for (int i = 0; i < 500; i++) {
        double w = 0.01 * pow(10.0, 3.0 * (double)i / 499.0);
        double r, im2;
        tf_eval_complex(tf, w, &r, &im2);
        double ph = atan2(im2, r) * 180.0 / M_PI;
        double err = fabs(ph + 180.0);
        if (err < min_phase_err) { min_phase_err = err; w180 = w; }
    }
    if (w180 > 0.0) {
        double r180, im180;
        tf_eval_complex(tf, w180, &r180, &im180);
        double mag180 = sqrt(r180 * r180 + im180 * im180);
        if (mag180 > 1e-10) sm->gain_margin = 1.0 / mag180;
        else sm->gain_margin = 1e9;
    } else {
        sm->gain_margin = 1e9;
    }

    /* Sector margin: max beta for which absolute stability holds */
    SectorBounds sb = ls->nonlinear;
    double m;
    if (circle_criterion_check(tf, sb, &m))
        sm->sector_margin = m;
    else
        sm->sector_margin = 0.0;

    tf_free(tf);
    return sm;
}

/* ===== Robust Stability Margin ===== */
StabilityMargin* lure_robust_margin(const LureSystem* ls, int method) {
    StabilityMargin* sm = (StabilityMargin*)calloc(1, sizeof(StabilityMargin));
    if (!sm || !ls) return sm;

    int n = ls->linear->n;
    TransferFunction* tf = tf_create(n + 1, n + 1);
    if (!tf) return sm;
    kyp_to_transfer_function(ls->linear, tf);

    /* Binary search for max beta */
    double alpha = ls->nonlinear.alpha;
    double beta_lo = alpha + 0.01;
    double beta_hi = LURE_MAX_SECTOR;
    double best_beta = alpha;

    for (int iter = 0; iter < 30; iter++) {
        double beta_mid = 0.5 * (beta_lo + beta_hi);
        SectorBounds sb = {alpha, beta_mid, ls->nonlinear.is_time_varying, 0, 0.0, 0.0};
        double margin;
        int stable = 0;
        if (method == LURE_METHOD_CIRCLE)
            stable = circle_criterion_check(tf, sb, &margin);
        else {
            double e;
            popov_find_optimal_slope(tf, sb, &e, &margin);
            stable = (margin > 0.0) ? 1 : 0;
        }
        if (stable) {
            best_beta = beta_mid;
            beta_lo = beta_mid;
        } else {
            beta_hi = beta_mid;
        }
        if (beta_hi - beta_lo < 0.001) break;
    }

    sm->sector_margin = best_beta - alpha;
    tf_free(tf);
    return sm;
}

void lure_margin_free(StabilityMargin* sm) { free(sm); }

/* ===== Loop Transformation ===== */
LureSystem* lure_loop_transform(const LureSystem* ls) {
    if (!ls) return NULL;
    double alpha = ls->nonlinear.alpha;
    double k = ls->nonlinear.beta - alpha;
    if (k <= 0.0) return NULL;

    /* G'(s) = G(s)/(1 + alpha*G(s))
     * For alpha=0: G' = G directly */
    int n = ls->linear->n;
    TransferFunction* tf = tf_create(n + 1, n + 1);
    if (!tf) return NULL;
    kyp_to_transfer_function(ls->linear, tf);

    TransferFunction* tf2 = popov_loop_transform(tf, ls->nonlinear, &k);
    tf_free(tf);
    if (!tf2) return NULL;

    /* Build new state-space from transformed transfer function.
     * For simplicity, use controllable canonical form. */
    int n2 = tf2->n_den - 1;
    if (n2 < 1) { tf_free(tf2); return NULL; }
    StateSpace* ss2 = ss_create(n2);
    if (!ss2) { tf_free(tf2); return NULL; }

    /* Controllable canonical form for strictly proper part:
     * A = companion matrix from denominator (excluding leading coeff)
     * B = [0 ... 0 1]^T
     * C = numerator coefficients (proper part)
     * D = 0 (handled separately)
     *
     * den(s) = den[0] + den[1]*s + ... + den[n]*s^n  (den[n] may not be 1)
     * Normalize by den[n] */

    double lead = tf2->den[n2];
    if (fabs(lead) < 1e-15) { ss_free(ss2); tf_free(tf2); return NULL; }

    /* Build companion matrix (bottom rows) */
    for (int i = 0; i < n2 - 1; i++)
        ss2->A[i * n2 + (i + 1)] = 1.0;
    for (int j = 0; j < n2; j++)
        ss2->A[(n2 - 1) * n2 + j] = -tf2->den[j] / lead;

    /* B = [0 ... 0 1]^T */
    ss2->B[n2 - 1] = 1.0;

    /* C = numerator coefficients (num[j] - D*den[j]) / lead
     * For strictly proper: D=0, C[j] = num[j]/lead */
    for (int j = 0; j < n2; j++)
        ss2->C[j] = (j < tf2->n_num ? tf2->num[j] : 0.0) / lead;

    ss2->D = 0.0;

    /* New sector: [0, k] */
    SectorBounds sb_new = {0.0, k, ls->nonlinear.is_time_varying, 0, 0.0, 0.0};
    LureSystem* ls_new = lure_create(ss2, sb_new);
    if (!ls_new) { ss_free(ss2); tf_free(tf2); return NULL; }

    tf_free(tf2);
    return ls_new;
}

/* ===== Zames-Falb Multiplier ===== */
LureSystem* lure_zames_falb_multiplier(const LureSystem* ls,
                                         double *multiplier_gain) {
    *multiplier_gain = 1.0;
    if (!ls) return NULL;
    /* Zames-Falb multipliers generalize Popov by allowing more
     * complex convolution kernels. For simplicity, return a copy
     * with unity multiplier (equivalent to standard Popov). */
    return lure_deep_copy(ls);
}

/* ===== H-infinity Norm (via frequency sweep) ===== */
double lure_hinf_norm(const LureSystem* ls) {
    if (!ls || !ls->linear) return 0.0;
    int n = ls->linear->n;
    TransferFunction* tf = tf_create(n + 1, n + 1);
    if (!tf) return 0.0;
    kyp_to_transfer_function(ls->linear, tf);
    double max_mag = 0.0;
    for (int i = 0; i < 500; i++) {
        double w = 0.01 * pow(10.0, 3.0 * (double)i / 499.0);
        double re, im;
        tf_eval_complex(tf, w, &re, &im);
        double mag = sqrt(re * re + im * im);
        if (mag > max_mag) max_mag = mag;
    }
    tf_free(tf);
    return max_mag;
}

/* Forward declaration of local Leverrier implementation */
static int lure_leverrier(const double *A, int n, double *poly);

/* ===== Check if Linear Part is Stable ===== */
bool lure_linear_is_stable(const LureSystem* ls) {
    if (!ls || !ls->linear) return false;
    int n = ls->linear->n;
    double *poly = (double*)calloc((size_t)(n + 1), sizeof(double));
    if (!poly) return false;
    if (lure_leverrier(ls->linear->A, n, poly) != 0) { free(poly); return false; }
    bool ok = true;
    for (int i = 0; i <= n && ok; i++)
        if (poly[i] <= 0.0) ok = false;
    free(poly);
    return ok;
}

/* ===== Self-Test ===== */
int lure_system_self_test(void) {
    /* Build a simple Lur'e system:
     * Linear part: G(s) = 1/(s+1)  (1st order, stable)
     * Nonlinearity: phi in sector [0, 2]
     *
     * State-space: A=-1, B=1, C=1, D=0 */
    StateSpace* ss = ss_create(1);
    assert(ss != NULL);
    ss->A[0] = -1.0; ss->B[0] = 1.0; ss->C[0] = 1.0; ss->D = 0.0;

    SectorBounds sb = {0.0, 2.0, 0, 0, 0.0, 0.0};
    LureSystem* ls = lure_create(ss, sb);
    assert(ls != NULL);
    assert(lure_sector_lower_bound(ls) == 0.0);
    assert(lure_sector_upper_bound(ls) == 2.0);

    /* Sector check: sigma=1, phi=1 => 0 <= 1 <= 2 => pass */
    assert(lure_sector_check(ls, 1.0, 1.0));
    /* sigma=1, phi=3 => 0 <= 3 > 2 => fail */
    assert(!lure_sector_check(ls, 1.0, 3.0));

    /* Deep copy */
    LureSystem* ls2 = lure_deep_copy(ls);
    assert(ls2 != NULL);
    assert(ls2->linear != ls->linear);  /* different pointers */
    assert(ls2->linear->A[0] == -1.0);
    {
        StateSpace* ss_saved = ls2->linear;  /* save before freeing ls2 */
        lure_free(ls2);
        ss_free(ss_saved);  /* cleanup deep copy linear part */
    }
    {
        LureSystem* ls3 = lure_deep_copy(ls);
        assert(ls3 != NULL);
        StateSpace* ss3 = ls3->linear;  /* save before freeing ls3 */
        lure_free(ls3);
        ss_free(ss3);
    }

    /* Circle criterion test */
    double margin;
    int result = lure_circle_test(ls, &margin);
    assert(result == 1);
    assert(margin > 0.0);

    /* Popov criterion test */
    double eta_opt;
    result = lure_popov_test(ls, &margin, &eta_opt);
    assert(result == 1);
    assert(margin > 0.0);
    assert(eta_opt >= 0.0);

    /* Unified test */
    result = lure_absolute_stability(ls, LURE_METHOD_CIRCLE, &margin, NULL);
    assert(result == 1);
    result = lure_absolute_stability(ls, LURE_METHOD_POPOV, &margin, &eta_opt);
    assert(result == 1);
    result = lure_absolute_stability(ls, LURE_METHOD_ALL, &margin, NULL);
    assert(result == 1);

    /* Gain/phase margins */
    StabilityMargin* sm = lure_gain_margin(ls);
    assert(sm != NULL);
    assert(sm->gain_margin > 0.0);
    assert(sm->phase_margin > -180.0 && sm->phase_margin < 180.0);
    lure_margin_free(sm);

    /* Robust margin */
    StabilityMargin* rm = lure_robust_margin(ls, LURE_METHOD_CIRCLE);
    assert(rm != NULL);
    assert(rm->sector_margin > 0.0);
    lure_margin_free(rm);

    /* H-infinity norm */
    double hinf = lure_hinf_norm(ls);
    assert(hinf > 0.0);

    /* Linear stability check */
    assert(lure_linear_is_stable(ls));

    /* Loop transformation */
    LureSystem* ls_lt = lure_loop_transform(ls);
    if (ls_lt) {
        assert(ls_lt->nonlinear.alpha == 0.0);
        assert(ls_lt->nonlinear.beta == 2.0);
        StateSpace* sst = ls_lt->linear;
        lure_free(ls_lt);
        ss_free(sst);
    }

    /* Zames-Falb multiplier */
    double mg;
    LureSystem* ls_zf = lure_zames_falb_multiplier(ls, &mg);
    if (ls_zf) {
        StateSpace* ssz = ls_zf->linear;
        lure_free(ls_zf);
        ss_free(ssz);
    }

    lure_free(ls);
    ss_free(ss);
    return 1;
}

/* ===== Leverrier Algorithm (local copy for this translation unit) ===== */
static int lure_leverrier(const double *A, int n, double *poly) {
    if (!A || !poly || n < 1 || n > 20) return -1;
    double *S = (double*)calloc((size_t)(n * n), sizeof(double));
    double *T = (double*)calloc((size_t)(n * n), sizeof(double));
    if (!S || !T) { free(S); free(T); return -1; }
    for (int i = 0; i < n; i++) S[i * n + i] = 1.0;
    poly[n] = 1.0;
    for (int k = 1; k <= n; k++) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                T[i * n + j] = 0.0;
                for (int m = 0; m < n; m++)
                    T[i * n + j] += A[i * n + m] * S[m * n + j];
            }
        double tr = 0.0;
        for (int i = 0; i < n; i++) tr += T[i * n + i];
        poly[n - k] = -tr / (double)k;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                S[i * n + j] = T[i * n + j];
                if (i == j) S[i * n + j] += poly[n - k];
            }
    }
    free(S); free(T); return 0;
}
/* Replace the call in lure_linear_is_stable with local copy */
/* Note: the function below overrides the earlier definition */
bool lure_linear_is_stable_fixed(const LureSystem* ls) {
    if (!ls || !ls->linear) return false;
    int n = ls->linear->n;
    double *poly = (double*)calloc((size_t)(n + 1), sizeof(double));
    if (!poly) return false;
    if (lure_leverrier(ls->linear->A, n, poly) != 0) { free(poly); return false; }
    bool ok = true;
    for (int i = 0; i <= n && ok; i++)
        if (poly[i] <= 0.0) ok = false;
    free(poly);
    return ok;
}
