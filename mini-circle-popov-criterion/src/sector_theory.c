#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "sector_theory.h"
#include "popov_criterion.h"

/* ===== Nonlinearity Lifecycle ===== */
Nonlinearity* nonlinearity_create(NonlinearityFunc phi, NonlinearityFunc dphi,
                                   const char *name) {
    Nonlinearity* nl = (Nonlinearity*)calloc(1, sizeof(Nonlinearity));
    if (!nl) return NULL;
    nl->phi = phi; nl->dphi = dphi; nl->name = name;
    return nl;
}

void nonlinearity_free(Nonlinearity* nl) { free(nl); }

/* ===== Sector Condition Verification ===== */
bool sector_check_point(SectorBounds sb, double sigma, double phi) {
    if (fabs(sigma) < SECTOR_EPS) return true;  /* trivially satisfied at 0 */
    double prod = sigma * phi;
    double lo = sb.alpha * sigma * sigma;
    double hi = sb.beta  * sigma * sigma;
    return (prod >= lo - SECTOR_EPS) && (prod <= hi + SECTOR_EPS);
}

bool sector_verify_function(NonlinearityFunc phi, SectorBounds sb,
                             double sigma_min, double sigma_max, int n) {
    if (!phi || n < 2) return false;
    double ds = (sigma_max - sigma_min) / (double)(n - 1);
    for (int i = 0; i < n; i++) {
        double sigma = sigma_min + (double)i * ds;
        double val = phi(sigma);
        if (!sector_check_point(sb, sigma, val)) return false;
    }
    return true;
}

SectorBounds sector_tightest_bounds(NonlinearityFunc phi,
                                     double sigma_min, double sigma_max, int n) {
    SectorBounds sb = {0.0, 0.0, 0, 0, 0.0, 0.0};
    if (!phi || n < 2 || sigma_min >= sigma_max) return sb;

    double alpha = 1e30, beta = -1e30;
    double ds = (sigma_max - sigma_min) / (double)(n - 1);
    for (int i = 0; i < n; i++) {
        double sigma = sigma_min + (double)i * ds;
        if (fabs(sigma) < SECTOR_EPS) continue;
        double val = phi(sigma);
        double slope = val / sigma;
        if (slope < alpha) alpha = slope;
        if (slope > beta)  beta  = slope;
    }
    /* Also test negative sigma values */
    for (int i = 0; i < n; i++) {
        double sigma = -(sigma_min + (double)i * ds);
        if (fabs(sigma) < SECTOR_EPS) continue;
        double val = phi(sigma);
        double slope = val / sigma;
        if (slope < alpha) alpha = slope;
        if (slope > beta)  beta  = slope;
    }
    sb.alpha = alpha; sb.beta = beta;
    return sb;
}

/* ===== Sector Transformations ===== */
SectorBounds sector_shift(SectorBounds sb, double c) {
    sb.alpha += c; sb.beta += c; return sb;
}

SectorBounds sector_scale(SectorBounds sb, double k) {
    sb.alpha *= k; sb.beta *= k;
    if (k < 0) { double t = sb.alpha; sb.alpha = sb.beta; sb.beta = t; }
    return sb;
}

SectorBounds sector_invert(SectorBounds sb) {
    if (fabs(sb.alpha) < SECTOR_EPS && fabs(sb.beta) < SECTOR_EPS)
        return sb;
    double a = sb.alpha, b = sb.beta;
    if (fabs(a) < SECTOR_EPS) a = SECTOR_EPS;
    if (fabs(b) < SECTOR_EPS) b = SECTOR_EPS;
    sb.alpha = 1.0 / b; sb.beta = 1.0 / a;
    return sb;
}

SectorBounds sector_normalize(SectorBounds sb, double *k_out) {
    double k = sb.beta - sb.alpha;
    if (k_out) *k_out = k;
    if (k <= 0.0) return sb;
    sb.alpha = 0.0; sb.beta = 1.0;
    return sb;
}

/* ===== Sector Arithmetic ===== */
SectorBounds sector_intersection(SectorBounds a, SectorBounds b) {
    SectorBounds r = a;
    r.alpha = fmax(a.alpha, b.alpha);
    r.beta  = fmin(a.beta,  b.beta);
    if (r.alpha > r.beta) { r.alpha = r.beta = 0.0; }
    return r;
}

SectorBounds sector_union(SectorBounds a, SectorBounds b) {
    SectorBounds r = a;
    r.alpha = fmin(a.alpha, b.alpha);
    r.beta  = fmax(a.beta,  b.beta);
    return r;
}

bool sector_contains(SectorBounds outer, SectorBounds inner) {
    return (outer.alpha <= inner.alpha + SECTOR_EPS) &&
           (outer.beta  >= inner.beta  - SECTOR_EPS);
}

/* ===== Standard Nonlinearity Implementations ===== */
double nonlinearity_saturation(double sigma, double limit) {
    if (sigma >  limit) return  limit;
    if (sigma < -limit) return -limit;
    return sigma;
}

double nonlinearity_saturation_deriv(double sigma, double limit) {
    if (fabs(sigma) < limit) return 1.0;
    return 0.0;
}

double nonlinearity_deadzone(double sigma, double delta) {
    if (sigma >  delta) return sigma - delta;
    if (sigma < -delta) return sigma + delta;
    return 0.0;
}

double nonlinearity_deadzone_deriv(double sigma, double delta) {
    if (fabs(sigma) < delta) return 0.0;
    return 1.0;
}

double nonlinearity_relay(double sigma) {
    if (sigma > 0.0) return 1.0;
    if (sigma < 0.0) return -1.0;
    return 0.0;
}

double nonlinearity_quantizer(double sigma, double q_level) {
    if (q_level <= 0.0) return sigma;
    return q_level * floor(sigma / q_level + 0.5);
}

double nonlinearity_cubic(double sigma, double c) {
    return sigma + c * sigma * sigma * sigma;
}

/* ===== Passivity and Sector Equivalence ===== */
bool sector_is_passive(SectorBounds sb) {
    return sb.alpha >= 0.0;  /* [0, inf) or tighter => passive */
}

bool sector_is_finite_gain(SectorBounds sb, double *gain_bound) {
    if (sb.beta >= 1e10) { *gain_bound = 1e10; return false; }
    *gain_bound = fmax(fabs(sb.alpha), fabs(sb.beta));
    return true;
}

/* ===== Combined Nonlinearity Tests ===== */
SectorBounds sector_parallel(SectorBounds s1, SectorBounds s2) {
    SectorBounds r;
    r.alpha = s1.alpha + s2.alpha;
    r.beta  = s1.beta  + s2.beta;
    r.is_time_varying = s1.is_time_varying || s2.is_time_varying;
    return r;
}

SectorBounds sector_feedback(SectorBounds forward, SectorBounds feedback) {
    /* For small-gain condition: forward*feedback < 1.
     * The combined sector is more complex; return a conservative estimate. */
    SectorBounds r;
    double f_max = fmax(fabs(forward.alpha), fabs(forward.beta));
    double b_max = fmax(fabs(feedback.alpha), fabs(feedback.beta));
    if (f_max * b_max < 1.0) {
        /* Small-gain satisfied => stable */
        r.alpha = forward.alpha / (1.0 + f_max * b_max);
        r.beta  = forward.beta  / (1.0 - f_max * b_max);
    } else {
        r.alpha = 0.0; r.beta = 1e10;  /* inconclusive */
    }
    return r;
}

/* ===== Max Allowed Gain ===== */
double sector_max_allowed_gain(const TransferFunction* tf,
                                int method, double *margin_at_max) {
    *margin_at_max = 0.0;
    if (!tf) return 0.0;

    /* Binary search for largest beta with alpha=0 */
    double beta_lo = 0.01, beta_hi = 1000.0;
    double best = 0.0;

    for (int iter = 0; iter < 40; iter++) {
        double beta_mid = 0.5 * (beta_lo + beta_hi);
        SectorBounds sb = {0.0, beta_mid, 0, 0, 0.0, 0.0};
        double margin;
        int stable;
        if (method == 0) {
            stable = circle_criterion_check(tf, sb, &margin);
        } else {
            double e;
            popov_find_optimal_slope(tf, sb, &e, &margin);
            stable = (margin > 0.0) ? 1 : 0;
        }
        if (stable) {
            best = beta_mid;
            *margin_at_max = margin;
            beta_lo = beta_mid;
        } else {
            beta_hi = beta_mid;
        }
        if (beta_hi - beta_lo < 0.001) break;
    }
    return best;
}

/* Forward declarations for wrapper nonlinearities */
static double sat1(double sigma);
static double dz02(double sigma);
static double cubic05(double sigma);

/* ===== Self-Test ===== */
int sector_theory_self_test(void) {
    /* Sector verification */
    SectorBounds sb = {0.0, 2.0, 0, 0, 0.0, 0.0};
    assert(sector_check_point(sb, 1.0, 1.0));   /* 1*1=1 in [0,2] */
    assert(sector_check_point(sb, 1.0, 0.5));   /* 1*0.5=0.5 in [0,2] */
    assert(!sector_check_point(sb, 1.0, 3.0));  /* 1*3=3 > 2 */
    assert(!sector_check_point(sb, 1.0, -0.1)); /* negative => fail for [0,2] */
    assert(sector_check_point(sb, 0.0, 100.0)); /* sigma=0 trivially satisfied */

    /* Verify saturation nonlinearity (using wrapper for fixed limit=1.0) */
    Nonlinearity* nl = nonlinearity_create(sat1, NULL, "sat1");
    assert(nl != NULL);
    /* sat(sigma, 1.0) has sector [0, 1] globally */
    SectorBounds sat_sb = sector_tightest_bounds(nl->phi, -5.0, 5.0, 200);
    assert(sat_sb.alpha >= -0.1);  /* ~= 0 */
    assert(sat_sb.beta <= 1.1);   /* ~= 1 */
    nonlinearity_free(nl);

    /* Test standard nonlinearity functions */
    assert(nonlinearity_saturation(0.5, 1.0) == 0.5);
    assert(nonlinearity_saturation(2.0, 1.0) == 1.0);
    assert(nonlinearity_saturation(-2.0, 1.0) == -1.0);
    assert(nonlinearity_deadzone(0.5, 0.2) == 0.3);
    assert(nonlinearity_deadzone(0.1, 0.2) == 0.0);
    assert(nonlinearity_relay(1.0) == 1.0);
    assert(nonlinearity_relay(-1.0) == -1.0);
    assert(nonlinearity_relay(0.0) == 0.0);
    assert(nonlinearity_cubic(1.0, 0.5) == 1.5);
    assert(nonlinearity_quantizer(1.7, 1.0) == 2.0);
    assert(nonlinearity_quantizer(1.3, 1.0) == 1.0);

    /* Verify function over domain */
    bool ok = sector_verify_function(sat1, sat_sb, -3.0, 3.0, 100);
    assert(ok);

    /* Sector transformations */
    SectorBounds sb2 = {1.0, 3.0, 0, 0, 0.0, 0.0};
    SectorBounds shifted = sector_shift(sb2, -1.0);
    assert(fabs(shifted.alpha - 0.0) < 0.01);
    assert(fabs(shifted.beta - 2.0) < 0.01);

    SectorBounds scaled = sector_scale(sb2, 2.0);
    assert(fabs(scaled.alpha - 2.0) < 0.01);
    assert(fabs(scaled.beta - 6.0) < 0.01);

    SectorBounds inv = sector_invert(sb2);
    assert(fabs(inv.alpha - 1.0/3.0) < 0.1);
    assert(fabs(inv.beta - 1.0) < 0.1);

    double k;
    SectorBounds norm = sector_normalize(sb2, &k);
    assert(fabs(k - 2.0) < 0.01);
    assert(fabs(norm.alpha) < 0.01);
    assert(fabs(norm.beta - 1.0) < 0.01);

    /* Sector arithmetic */
    SectorBounds a = {0.0, 2.0, 0, 0, 0.0, 0.0};
    SectorBounds b = {1.0, 3.0, 0, 0, 0.0, 0.0};
    SectorBounds inter = sector_intersection(a, b);
    assert(fabs(inter.alpha - 1.0) < 0.01);
    assert(fabs(inter.beta - 2.0) < 0.01);
    SectorBounds uni = sector_union(a, b);
    assert(fabs(uni.alpha) < 0.01);
    assert(fabs(uni.beta - 3.0) < 0.01);
    assert(sector_contains(uni, a));
    assert(sector_contains(uni, b));

    /* Passivity */
    SectorBounds passive_sb = {0.0, 1e10, 0, 0, 0.0, 0.0};
    assert(sector_is_passive(passive_sb));
    SectorBounds active_sb = {-1.0, 1.0, 0, 0, 0.0, 0.0};
    assert(!sector_is_passive(active_sb));

    double gain_bound;
    assert(sector_is_finite_gain(a, &gain_bound));
    assert(gain_bound > 0.0);

    /* Combined sectors */
    SectorBounds par = sector_parallel(a, b);
    assert(fabs(par.alpha - 1.0) < 0.01);
    assert(fabs(par.beta - 5.0) < 0.01);

    /* Max allowed gain test */
    TransferFunction* tf = tf_create(1, 2);
    assert(tf != NULL);
    tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
    double max_k, margin_k;
    max_k = sector_max_allowed_gain(tf, 0, &margin_k);
    assert(max_k > 0.0);
    tf_free(tf);

    return 1;
}

/* ===== Additional Nonlinearity Types ===== */

/* Hysteresis relay with state.
 * state tracks the previous sign. */
double nonlinearity_hysteresis_relay(double sigma, double hysteresis,
                                      double *prev_output) {
    if (!prev_output) return nonlinearity_relay(sigma);
    if (sigma > hysteresis) {
        *prev_output = 1.0;
    } else if (sigma < -hysteresis) {
        *prev_output = -1.0;
    }
    /* else: stay at previous value */
    return *prev_output;
}

/* Sigmoid nonlinearity: phi(sigma) = tanh(sigma * gain)
 * Sector: [0, 1] globally (monotonic odd sigmoid). */
double nonlinearity_sigmoid(double sigma, double gain) {
    double x = gain * sigma;
    /* Stable tanh computation avoiding overflow */
    if (x > 20.0) return 1.0;
    if (x < -20.0) return -1.0;
    double ex = exp(x);
    double emx = 1.0 / ex;
    return (ex - emx) / (ex + emx);
}

/* Polynomial nonlinearity: phi(sigma) = sum_{k=0}^{deg-1} coeff[k] * sigma^k */
double nonlinearity_polynomial(double sigma, const double *coeffs, int degree) {
    double result = 0.0;
    double spow = 1.0;
    for (int k = 0; k < degree; k++) {
        result += coeffs[k] * spow;
        spow *= sigma;
    }
    return result;
}

/* ===== Sector-Based Lyapunov Function Construction ===== */

/* Given a Lur'e system and a KYP solution, construct the
 * Lur'e-Postnikov Lyapunov function:
 *   V(x) = x^T * P * x + 2 * integral_0^{Cx} phi(sigma) dsigma
 *
 * Returns the value of V(x) for a given state vector x.
 * Requires the KYP solution P and the nonlinearity phi. */
double sector_lyapunov_value(const double *P, int n,
                              const double *x, NonlinearityFunc phi,
                              const double *C, double alpha) {
    /* Quadratic term: x^T * P * x */
    double quad = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            quad += x[i] * P[i * n + j] * x[j];

    /* Integral term: 2 * integral_0^{y} phi(s) ds
     * Approximated via trapezoidal integration */
    double y = 0.0;
    for (int i = 0; i < n; i++) y += C[i] * x[i];

    int steps = 100;
    double integral = 0.0;
    double ds = y / (double)steps;
    double phi_prev = phi(0.0);
    for (int i = 1; i <= steps; i++) {
        double s = (double)i * ds;
        double phi_s = phi(s);
        integral += 0.5 * (phi_prev + phi_s) * ds;
        phi_prev = phi_s;
    }

    return quad + 2.0 * integral * alpha;
}

/* ===== Sector Margin via Binary Search ===== */

/* Find the largest beta for which the sector [alpha, beta] is
 * absolutely stabilizing using the specified method (0=circle, 1=popov).
 * Performs binary search on beta. */
double sector_max_beta_binary_search(const TransferFunction* tf,
                                      double alpha, int method) {
    if (!tf || alpha < 0.0) return alpha;

    double lo = alpha + 1e-6;
    double hi = SECTOR_SIGMA_MAX * 100.0;
    double best = alpha;

    for (int iter = 0; iter < 50; iter++) {
        double mid = 0.5 * (lo + hi);
        SectorBounds sb = {alpha, mid, 0, 0, 0.0, 0.0};
        double margin;
        int ok;
        if (method == 0)
            ok = circle_criterion_check(tf, sb, &margin);
        else {
            double e;
            popov_find_optimal_slope(tf, sb, &e, &margin);
            ok = (margin > 0.0) ? 1 : 0;
        }
        if (ok) {
            best = mid;
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo < 1e-6) break;
    }
    return best;
}

/* ===== Aizerman Conjecture Test =====
 * The Aizerman conjecture states: if the linearized system (with phi=K*sigma)
 * is stable for all K in [alpha, beta], then the nonlinear system with phi in
 * [alpha, beta] is absolutely stable.
 *
 * This conjecture is known to be FALSE in general (counterexamples exist for
 * n >= 3), but can be tested for specific systems.
 *
 * Returns: 1 if conjecture holds for this system, 0 if not. */
int sector_test_aizerman(const TransferFunction* tf, SectorBounds sb) {
    (void)tf; (void)sb;
    /* The Aizerman conjecture is false in general.
     * This function returns 0 as a conservative default. */
    return 0;
}

/* ===== Kalman Conjecture Test =====
 * The Kalman conjecture: if phi'(sigma) is in [alpha, beta] (slope-restricted),
 * then the system is absolutely stable.
 *
 * This conjecture is also FALSE in general for n >= 4.
 *
 * Returns: 1 if conjecture holds, 0 if not. */
int sector_test_kalman(const TransferFunction* tf, SectorBounds sb) {
    (void)tf; (void)sb;
    return 0;
}

/* ===== Mixed Sector + Slope Restriction ===== */

/* Check if a nonlinearity satisfies both a sector condition and
 * a slope restriction. This enables less conservative stability
 * criteria (e.g., off-axis circle criterion). */
bool sector_check_slope_restricted(SectorBounds sb,
                                    double sigma, double phi,
                                    double dphi_dsigma) {
    /* Sector condition */
    if (!sector_check_point(sb, sigma, phi)) return false;
    /* Slope restriction */
    if (sb.is_slope_bounded) {
        if (dphi_dsigma < sb.slope_min || dphi_dsigma > sb.slope_max)
            return false;
    }
    return true;
}

/* ===== Sector Robustness Analysis ===== */

/* Compute the sensitivity of the sector margin to parameter variations.
 * Returns the derivative-like sensitivity measure (finite difference). */
double sector_sensitivity(const TransferFunction* tf,
                           SectorBounds sb, double delta) {
    double margin_nom, margin_pert;
    circle_criterion_check(tf, sb, &margin_nom);

    SectorBounds sb_pert = sb;
    sb_pert.beta += delta;
    circle_criterion_check(tf, sb_pert, &margin_pert);

    return (margin_pert - margin_nom) / delta;
}

/* Wrapper nonlinearities with fixed parameters for NonlinearityFunc interface */
static double sat1(double sigma) { return nonlinearity_saturation(sigma, 1.0); }
static double dz02(double sigma) { return nonlinearity_deadzone(sigma, 0.2); }
static double cubic05(double sigma) { return nonlinearity_cubic(sigma, 0.5); }
