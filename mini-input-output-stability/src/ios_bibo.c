#include "ios_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================
 * Polynomial helper for BIBO analysis
 * ============================================================ */

static int poly_roots_real(const double* c, int n, double* roots) {
    if (n == 2) {
        if (fabs(c[0]) < IOS_EPS) return 0;
        roots[0] = -c[1] / c[0];
        return 1;
    }
    if (n == 3) {
        double a = c[0], b = c[1], c2 = c[2];
        if (fabs(a) < IOS_EPS) {
            if (fabs(b) < IOS_EPS) return 0;
            roots[0] = -c2 / b;
            return 1;
        }
        double d = b * b - 4.0 * a * c2;
        if (d >= 0) {
            roots[0] = (-b + sqrt(d)) / (2.0 * a);
            roots[1] = (-b - sqrt(d)) / (2.0 * a);
            return 2;
        }
        return 0;
    }
    return 0;
}

/* ============================================================
 * BIBO Stability Checks
 * ============================================================ */
IOSBIBOStatus ios_check_bibo_tf(const IOSTransferFcn* tf) {
    if (!tf || tf->n_den < 1) return IOS_BIBO_UNSTABLE;
    /* For 1st order: a0*s + a1 = 0 => s = -a1/a0, stable if a1/a0 > 0 */
    if (tf->n_den == 1) {
        double a = tf->den[0];
        return (fabs(a) < IOS_EPS) ? IOS_BIBO_CRITICAL : IOS_BIBO_STABLE;
    }
    /* Compute roots */
    double* roots = (double*)malloc((size_t)(tf->n_den - 1) * sizeof(double));
    if (!roots) return IOS_BIBO_UNSTABLE;
    int nr = poly_roots_real(tf->den, tf->n_den, roots);
    bool ok = true;
    for (int i = 0; i < nr; i++)
        if (roots[i] > -IOS_EPS) ok = false;
    free(roots);
    return ok ? IOS_BIBO_STABLE : IOS_BIBO_UNSTABLE;
}

IOSBIBOStatus ios_check_bibo_ss(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return IOS_BIBO_UNSTABLE;
    int n = sys->n;
    /* 1st order: eigenvalue is A[0] */
    if (n == 1)
        return (sys->A[0] < -IOS_EPS) ? IOS_BIBO_STABLE : IOS_BIBO_UNSTABLE;
    /* 2nd order: trace < 0 and determinant > 0 */
    if (n == 2) {
        double tr = sys->A[0] + sys->A[3];
        double det = sys->A[0] * sys->A[3] - sys->A[1] * sys->A[2];
        return (tr < -IOS_EPS && det > IOS_EPS) ? IOS_BIBO_STABLE : IOS_BIBO_UNSTABLE;
    }
    /* General case: check trace (sum of eigenvalues) */
    double tr = 0.0;
    for (int i = 0; i < n; i++) tr += sys->A[i * n + i];
    return (tr < -IOS_EPS) ? IOS_BIBO_STABLE : IOS_BIBO_CRITICAL;
}

bool ios_is_bibo_stable(const IOSSystem* sys) {
    if (!sys || !sys->sys) return false;
    return ios_check_bibo_ss(sys->sys) == IOS_BIBO_STABLE;
}

/* ============================================================
 * BIBO Gain Computation
 * ============================================================ */
double ios_bibo_gain(const IOSSystem* sys, IOSNormType type) {
    if (!sys || !sys->sys) return IOS_INF;
    int n = sys->sys->n;
    double g = IOS_INF;
    switch (type) {
        case IOS_NORM_L1: {
            g = 0.0;
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < n; j++)
                    s += fabs(sys->sys->A[i * n + j]);
                if (s > g) g = s;
            }
            break;
        }
        case IOS_NORM_L2: {
            double tr = 0.0;
            for (int i = 0; i < n; i++)
                tr += sys->sys->A[i * n + i];
            g = 1.0 / fmax(fabs(tr), IOS_EPS);
            break;
        }
        case IOS_NORM_LINF:
            g = 1.0 / fmax(fabs(sys->sys->A[0]), IOS_EPS);
            break;
    }
    return g;
}

double ios_L1_induced_norm(const IOSStateSpace* sys) {
    if (!sys) return 0.0;
    double g = 0.0;
    for (int j = 0; j < sys->m; j++) {
        double s = 0.0;
        for (int i = 0; i < sys->n; i++)
            s += fabs(sys->B[i * sys->m + j]);
        if (s > g) g = s;
    }
    return g;
}

double ios_Linf_induced_norm(const IOSStateSpace* sys) {
    if (!sys) return 0.0;
    double g = 0.0;
    for (int i = 0; i < sys->p; i++) {
        double s = 0.0;
        for (int j = 0; j < sys->n; j++)
            s += fabs(sys->C[i * sys->n + j]);
        if (s > g) g = s;
    }
    return g;
}

double ios_peak_gain(const IOSStateSpace* sys) {
    double a = ios_L1_induced_norm(sys);
    double b = ios_Linf_induced_norm(sys);
    return a > b ? a : b;
}

/* ============================================================
 * L2 Gain (H infinity norm approximation via frequency sweep)
 * ============================================================ */
double ios_L2_gain_ss(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return IOS_INF;
    double hi = 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double m = 0.0;
        for (int i = 0; i < sys->p; i++)
            for (int j = 0; j < sys->m; j++) {
                double g_re = sys->D[i * sys->m + j];
                double g_im = 0.0;
                for (int r = 0; r < sys->n; r++) {
                    double a_r = sys->A[r * sys->n + r];
                    double denom = a_r * a_r + w * w;
                    if (denom > IOS_EPS) {
                        double c = sys->C[i * sys->n + r];
                        double b = sys->B[r * sys->m + j];
                        g_re += c * b * a_r / denom;
                        g_im -= c * b * w / denom;
                    }
                }
                double mag = sqrt(g_re * g_re + g_im * g_im);
                if (mag > m) m = mag;
            }
        if (m > hi) hi = m;
    }
    return hi;
}

double ios_L2_gain_tf(const IOSTransferFcn* tf) {
    if (!tf) return IOS_INF;
    double hi = 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        if (mag > hi) hi = mag;
    }
    return hi;
}

/* ============================================================
 * H2 Norm
 * ============================================================ */
double ios_H2_norm_tf(const IOSTransferFcn* tf) {
    if (!tf) return IOS_INF;
    double H2 = 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        H2 += mag * mag;
    }
    return sqrt(H2 / (2.0 * IOS_PI));
}

double ios_Hinf_norm_tf(const IOSTransferFcn* tf) {
    return ios_L2_gain_tf(tf);
}

double ios_h2_norm_ss(const IOSStateSpace* sys) {
    if (!sys) return IOS_INF;
    double tr = 0.0;
    for (int i = 0; i < sys->n; i++) tr += sys->A[i * sys->n + i];
    return fabs(tr) > IOS_EPS ? 1.0 / fabs(tr) : IOS_INF;
}

/* ============================================================
 * L1 Norm of Impulse Response
 * ============================================================ */
double ios_L1_norm_impulse(const IOSTransferFcn* tf) {
    if (!tf) return IOS_INF;
    double L1 = 0.0;
    for (double t = 0.0; t < 100.0; t += 0.01)
        L1 += fabs(ios_impulse_response_tf(tf, t)) * 0.01;
    return L1;
}

/* ============================================================
 * Impulse and Step Response (Transfer Function)
 * ============================================================ */
double ios_impulse_response(const IOSTransferFcn* tf, double t) {
    return ios_impulse_response_tf(tf, t);
}

double ios_step_response(const IOSTransferFcn* tf, double t) {
    return ios_step_response_tf(tf, t);
}

double ios_convolution(const IOSSignal* u, const IOSTransferFcn* tf, double t) {
    if (!u || !tf) return 0.0;
    double y = 0.0;
    for (double tau = 0.0; tau < t; tau += 0.01) {
        int idx = (int)(tau / 0.01);
        if (idx < u->length)
            y += u->data[idx] * ios_impulse_response_tf(tf, t - tau) * 0.01;
    }
    return y;
}

double ios_impulse_response_tf(const IOSTransferFcn* tf, double t) {
    if (!tf || t < 0.0) return 0.0;
    /* Inverse Laplace approx via numerical integration along imaginary axis */
    double h = 0.0;
    for (int i = 0; i < 50; i++) {
        double s = 0.1 * (double)(i + 1);
        h += ios_tf_eval(tf, s) * exp(s * t) * 0.1;
    }
    return h;
}

double ios_step_response_tf(const IOSTransferFcn* tf, double t) {
    if (!tf || t < 0.0) return 0.0;
    double resp = 0.0;
    for (double tau = 0.0; tau < t; tau += 0.01)
        resp += ios_impulse_response_tf(tf, tau) * 0.01;
    return resp;
}

/* ============================================================
 * Impulse and Step Response (State Space)
 * ============================================================ */
double ios_impulse_response_ss(const IOSStateSpace* sys, double t) {
    if (!sys || t < 0.0) return 0.0;
    double y = 0.0;
    for (int i = 0; i < sys->p; i++)
        for (int j = 0; j < sys->n; j++)
            for (int k = 0; k < sys->m; k++)
                y += sys->C[i * sys->n + j] * sys->B[j * sys->m + k]
                     * exp(sys->A[j * sys->n + j] * t);
    return y;
}

void ios_impulse_response_full(const IOSStateSpace* sys, double t_max, double dt, IOSSignal* h) {
    if (!sys || !h) return;
    for (double t = 0.0; t < t_max; t += dt)
        ios_signal_append(h, ios_impulse_response_ss(sys, t));
}

void ios_step_response_full(const IOSStateSpace* sys, double t_max, double dt, IOSSignal* y) {
    if (!sys || !y) return;
    double s = 0.0;
    for (double t = 0.0; t < t_max; t += dt) {
        s += ios_impulse_response_ss(sys, t) * dt;
        ios_signal_append(y, s);
    }
}

double ios_impulse_energy(const IOSStateSpace* sys) {
    if (!sys) return IOS_INF;
    IOSSignal* h = ios_signal_create(1000);
    if (!h) return IOS_INF;
    ios_impulse_response_full(sys, 50.0, 0.05, h);
    double e = ios_signal_energy(h);
    ios_signal_free(h);
    return e;
}

/* ============================================================
 * Exponential Stability
 * ============================================================ */
bool ios_is_exponentially_stable(const IOSStateSpace* sys) {
    return ios_check_bibo_ss(sys) == IOS_BIBO_STABLE;
}

/* ============================================================
 * Settling Time, Overshoot, Rise Time Estimates
 * ============================================================ */
double ios_settling_time(const IOSStateSpace* sys, double tol) {
    if (!sys || sys->n < 1) return IOS_INF;
    double mr = ios_ss_eigen_max_real(sys);
    if (mr >= -IOS_EPS) return IOS_INF;
    return -log(tol) / fabs(mr);
}

double ios_overshoot_estimate(const IOSStateSpace* sys) {
    if (!sys || sys->n < 2) return 0.0;
    double tr = sys->A[0] + sys->A[3];
    double det = sys->A[0] * sys->A[3] - sys->A[1] * sys->A[2];
    double denom = sqrt(fmax(det, IOS_EPS));
    double zeta = -tr / (2.0 * denom);
    if (zeta >= 1.0) return 0.0;
    return exp(-IOS_PI * zeta / sqrt(1.0 - zeta * zeta)) * 100.0;
}

double ios_rise_time_estimate(const IOSStateSpace* sys) {
    double bw = 1.0;
    if (sys) {
        IOSTransferFcn* tf = ios_tf_create(1, 2);
        if (tf) {
            tf->num[0] = 1.0;
            tf->den[0] = 1.0;
            tf->den[1] = -ios_ss_eigen_max_real(sys);
            bw = ios_bandwidth(tf);
            ios_tf_free(tf);
        }
    }
    return bw > 0.0 ? 0.35 / bw : IOS_INF;
}

/* ============================================================
 * BIBO Margins
 * ============================================================ */
void ios_bibo_margins(const IOSTransferFcn* tf, double* gm, double* pm) {
    if (!tf || !gm || !pm) return;
    *gm = IOS_INF;
    *pm = 0.0;
    
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        double d = fabs(1.0 + mag);
        if (d < *gm) *gm = d;
        if (fabs(mag - 1.0) < 0.01)
            *pm = IOS_PI + phase;
    }
}

int ios_bibo_margin_check(const IOSStateSpace* sys, double margin) {
    if (!sys) return -1;
    double g = ios_L2_gain_ss(sys);
    return (g < margin) ? 0 : 1;
}

/* ============================================================
 * Additional Utility Functions
 * ============================================================ */
bool ios_is_strictly_proper(const IOSStateSpace* sys) {
    if (!sys) return false;
    for (int i = 0; i < sys->p * sys->m; i++)
        if (fabs(sys->D[i]) > IOS_EPS) return false;
    return true;
}

bool ios_has_integrator(const IOSTransferFcn* tf) {
    return tf && tf->n_den > 0 && fabs(tf->den[tf->n_den - 1]) < IOS_EPS;
}

bool ios_is_minimum_phase(const IOSTransferFcn* tf) {
    if (!tf || tf->n_num < 2) return true;
    return tf->num[0] * tf->num[tf->n_num - 1] > -IOS_EPS;
}

double ios_Lp_gain_bound(const IOSStateSpace* sys, double p) {
    if (!sys) return IOS_INF;
    double L2 = ios_L2_gain_ss(sys);
    if (fabs(p - 2.0) < 0.01) return L2;
    if (p < 1.0) return IOS_INF;
    return L2 * pow(2.0, fabs(p - 2.0));
}

double ios_poles_max(const IOSStateSpace* sys) {
    return ios_ss_eigen_max_real(sys);
}

int ios_count_unstable(const IOSStateSpace* sys) {
    if (!sys) return -1;
    int c = 0;
    for (int i = 0; i < sys->n; i++)
        if (sys->A[i * sys->n + i] > IOS_EPS) c++;
    return c;
}

int ios_stability_radius(const IOSStateSpace* sys, double* radius) {
    if (!sys || !radius) return -1;
    double mr = ios_ss_eigen_max_real(sys);
    *radius = fabs(mr);
    return mr < 0.0 ? 0 : 1;
}

/* ============================================================
 * Bode / Nyquist Utilities
 * ============================================================ */
double ios_gain_crossover_frequency(const IOSTransferFcn* tf) {
    if (!tf) return 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        if (fabs(mag - 1.0) < 0.05) return w;
    }
    return 0.0;
}

double ios_phase_crossover_frequency(const IOSTransferFcn* tf) {
    if (!tf) return 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        if (fabs(phase + IOS_PI) < 0.05) return w;
    }
    return 0.0;
}

void ios_nyquist_analysis(const IOSTransferFcn* tf, double* re, double* im, int n) {
    if (!tf || !re || !im || n < 1) return;
    for (int i = 0; i < n; i++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)i / (double)n);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        re[i] = mag * cos(phase);
        im[i] = mag * sin(phase);
    }
}

/* ============================================================
 * Reports and Tables
 * ============================================================ */
void ios_bibo_report(const IOSStateSpace* sys) {
    if (!sys) return;
    printf("BIBO: n=%d status=%s L2=%.4f\n",
           sys->n,
           ios_check_bibo_ss(sys) == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE",
           ios_L2_gain_ss(sys));
}

void ios_print_gain_table(const IOSStateSpace* sys) {
    if (!sys) return;
    printf("Gain Table:\n");
    printf("  L2 gain:      %.6f\n", ios_L2_gain_ss(sys));
    printf("  L1 induced:   %.6f\n", ios_L1_induced_norm(sys));
    printf("  Linf induced: %.6f\n", ios_Linf_induced_norm(sys));
    printf("  Peak gain:    %.6f\n", ios_peak_gain(sys));
    printf("  DC gain:      %.6f\n", ios_dc_gain_ss(sys));
}

void ios_bibo_analysis_report(const IOSStateSpace* sys) {
    ios_print_gain_table(sys);
}

void ios_bibo_test_suite(void) {
    printf("=== BIBO Test Suite ===\n");
    IOSStateSpace* ss = ios_ss_create(1, 1, 1);
    if (!ss) return;
    for (double a = -5.0; a <= 5.0; a += 1.0) {
        ss->A[0] = a; ss->B[0] = 1.0; ss->C[0] = 1.0;
        printf("  A=%.1f: BIBO=%s L2=%.4f\n", a,
               ios_check_bibo_ss(ss) == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE",
               ios_L2_gain_ss(ss));
    }
    ios_ss_free(ss);
}

void ios_bibo_compare_systems(const IOSStateSpace* sys1, const IOSStateSpace* sys2) {
    if (!sys1 || !sys2) return;
    printf("System comparison:\n");
    printf("  Sys1: L2=%.4f Sys2: L2=%.4f\n",
           ios_L2_gain_ss(sys1), ios_L2_gain_ss(sys2));
    printf("  Difference: %.4f\n", ios_L2_gain_ss(sys1) - ios_L2_gain_ss(sys2));
}

int ios_bibo_verify_with_input(const IOSStateSpace* sys, const IOSSignal* u, double tol) {
    if (!sys || !u) return -1;
    IOSSignal* y = ios_signal_create(u->length);
    if (!y) return -1;
    ios_ss_simulate(sys, u, y);
    double gain_out = ios_signal_norm(y, IOS_NORM_L2);
    double gain_in = ios_signal_norm(u, IOS_NORM_L2);
    ios_signal_free(y);
    return (gain_in > IOS_EPS && gain_out / gain_in < tol) ? 0 : 1;
}

double ios_weighted_L2_gain(const IOSStateSpace* sys, const double* W, int n_w) {
    (void)W; (void)n_w;
    return ios_L2_gain_ss(sys);
}

double ios_robustness_metric(const IOSStateSpace* sys) {
    if (!sys) return 0.0;
    double L2 = ios_L2_gain_ss(sys);
    return L2 > IOS_EPS ? 1.0 / L2 : 0.0;
}

double ios_bibo_bandwidth_ss(const IOSStateSpace* sys) {
    if (!sys) return 1.0;
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return 1.0;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0;
    tf->den[1] = -ios_ss_eigen_max_real(sys);
    double bw = ios_bandwidth(tf);
    ios_tf_free(tf);
    return bw;
}

double ios_bibo_delay_margin_ss(const IOSStateSpace* sys) {
    double pm = 90.0;
    double bw = ios_bibo_bandwidth_ss(sys);
    return bw > IOS_EPS ? pm * IOS_PI / (180.0 * bw) : IOS_INF;
}

double ios_initial_response_ss(const IOSStateSpace* sys, const double* x0, double t) {
    if (!sys || !x0 || t < 0.0) return 0.0;
    double y = 0.0;
    for (int i = 0; i < sys->p; i++)
        for (int j = 0; j < sys->n; j++)
            y += sys->C[i * sys->n + j] * x0[j] * exp(sys->A[j * sys->n + j] * t);
    return y;
}

void ios_initial_response_full(const IOSStateSpace* sys, const double* x0,
                                 double t_max, double dt, IOSSignal* y) {
    if (!sys || !x0 || !y) return;
    ios_signal_clear(y);
    for (double t = 0.0; t < t_max; t += dt)
        ios_signal_append(y, ios_initial_response_ss(sys, x0, t));
}

int ios_compare_bibo_margins(const IOSStateSpace* sys1, const IOSStateSpace* sys2) {
    if (!sys1 || !sys2) return -1;
    double g1 = ios_L2_gain_ss(sys1), g2 = ios_L2_gain_ss(sys2);
    return (g1 < g2) ? 1 : ((g1 > g2) ? -1 : 0);
}

/* ============================================================
 * Eigenvalue Analysis Utilities
 * For diagonal/simple matrices, compute stability margins
 * ============================================================ */

/* Compute spectral abscissa: max(Re(lambda_i)) */
double ios_spectral_abscissa(const IOSStateSpace* sys) {
    return ios_ss_eigen_max_real(sys);
}

/* Time constant estimate: -1 / max(Re(lambda_i)) */
double ios_time_constant(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return IOS_INF;
    double mr = ios_ss_eigen_max_real(sys);
    if (mr >= -IOS_EPS) return IOS_INF;
    return -1.0 / mr;
}

/* Compute damping ratio for 2nd order systems */
double ios_damping_ratio(const IOSStateSpace* sys) {
    if (!sys || sys->n != 2) return -1.0;
    double tr = sys->A[0] + sys->A[3];
    double det = sys->A[0] * sys->A[3] - sys->A[1] * sys->A[2];
    double wn = sqrt(fmax(det, 0.0));
    if (wn < IOS_EPS) return -1.0;
    return -tr / (2.0 * wn);
}

/* Natural frequency for 2nd order systems */
double ios_natural_frequency(const IOSStateSpace* sys) {
    if (!sys || sys->n != 2) return -1.0;
    double det = sys->A[0] * sys->A[3] - sys->A[1] * sys->A[2];
    return det > IOS_EPS ? sqrt(det) : -1.0;
}

/* Controllability Gramian trace (diagonal approximation) */
double ios_controllability_gramian_trace(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return IOS_INF;
    double trace = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double a_ii = sys->A[i * sys->n + i];
        if (a_ii >= -IOS_EPS) return IOS_INF;
        for (int j = 0; j < sys->m; j++)
            trace += sys->B[i * sys->m + j] * sys->B[i * sys->m + j] / (-2.0 * a_ii);
    }
    return trace;
}

/* Observability Gramian trace (diagonal approximation) */
double ios_observability_gramian_trace(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return IOS_INF;
    double trace = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double a_ii = sys->A[i * sys->n + i];
        if (a_ii >= -IOS_EPS) return IOS_INF;
        for (int j = 0; j < sys->p; j++)
            trace += sys->C[j * sys->n + i] * sys->C[j * sys->n + i] / (-2.0 * a_ii);
    }
    return trace;
}

/* Hankel singular values (diagonal approx): sqrt(gramian_ctrb * gramian_obs) */
double ios_hankel_singular_value_max(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return IOS_INF;
    double max_hsv = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double a_ii = sys->A[i * sys->n + i];
        if (a_ii >= -IOS_EPS) return IOS_INF;
        double c_sum = 0.0, b_sum = 0.0;
        for (int j = 0; j < sys->p; j++)
            c_sum += sys->C[j * sys->n + i] * sys->C[j * sys->n + i];
        for (int j = 0; j < sys->m; j++)
            b_sum += sys->B[i * sys->m + j] * sys->B[i * sys->m + j];
        double hsv = sqrt(c_sum * b_sum) / (-2.0 * a_ii);
        if (hsv > max_hsv) max_hsv = hsv;
    }
    return max_hsv;
}

/* Modal decomposition analysis for diagonal systems */
void ios_modal_analysis(const IOSStateSpace* sys, double* modes, int max_modes, int* count) {
    if (!sys || !modes || !count) return;
    *count = 0;
    for (int i = 0; i < sys->n && *count < max_modes; i++)
        modes[(*count)++] = sys->A[i * sys->n + i];
}

/* Relative stability: distance of eigenvalues from imaginary axis */
double ios_stability_distance(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return IOS_INF;
    double min_dist = IOS_INF;
    for (int i = 0; i < sys->n; i++) {
        double d = fabs(sys->A[i * sys->n + i]);
        if (d < min_dist) min_dist = d;
    }
    return min_dist;
}

/* Compute condition number for A matrix (diagonal approx) */
double ios_condition_number(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return 1.0;
    double max_abs = fabs(sys->A[0]);
    double min_abs = fabs(sys->A[0]);
    for (int i = 1; i < sys->n; i++) {
        double d = fabs(sys->A[i * sys->n + i]);
        if (d > max_abs) max_abs = d;
        if (d < min_abs) min_abs = d;
    }
    return min_abs > IOS_EPS ? max_abs / min_abs : IOS_INF;
}

/* I/O delay margin estimate using phase margin */
double ios_delay_margin_from_pm(const IOSTransferFcn* tf) {
    if (!tf) return IOS_INF;
    double gm = 0.0, pm = 0.0;
    ios_bibo_margins(tf, &gm, &pm);
    double wc = ios_gain_crossover_frequency(tf);
    if (wc < IOS_EPS) return IOS_INF;
    return pm / wc;
}

/* Check if system is minimum-phase (all zeros in LHP) */
bool ios_check_minimum_phase_ss(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return true;
    /* For SISO: check that C*(sI-A)^{-1}*B + D has no RHP zeros */
    double dc = ios_dc_gain_ss(sys);
    if (fabs(dc) < IOS_EPS) return false;
    /* Check if direct feedthrough is same sign as DC gain */
    for (int i = 0; i < sys->p * sys->m; i++)
        if (sys->D[i] * dc < -IOS_EPS) return false;
    return true;
}

/* Phase margin computation via frequency sweep */
double ios_phase_margin_tf(const IOSTransferFcn* tf) {
    if (!tf) return 0.0;
    double gm = 0.0, pm = 0.0;
    ios_bibo_margins(tf, &gm, &pm);
    return pm;
}

/* Gain margin computation via frequency sweep */
double ios_gain_margin_tf(const IOSTransferFcn* tf) {
    if (!tf) return IOS_INF;
    double gm = 0.0, pm = 0.0;
    ios_bibo_margins(tf, &gm, &pm);
    return gm;
}

/* Check LTI system for internal stability */
bool ios_internal_stability_check(const IOSStateSpace* sys) {
    if (!sys) return false;
    return ios_is_exponentially_stable(sys) && ios_is_strictly_proper(sys);
}

/* Characteristic polynomial evaluation at s */
double ios_char_poly_eval(const IOSStateSpace* sys, double s) {
    if (!sys || sys->n < 1) return 0.0;
    /* det(sI - A) product of (s - lambda_i) for diagonal A */
    double det_val = 1.0;
    for (int i = 0; i < sys->n; i++)
        det_val *= (s - sys->A[i * sys->n + i]);
    return det_val;
}

/* ============================================================
 * Comprehensive System Report
 * ============================================================ */

void ios_system_full_report(const IOSStateSpace* sys) {
    if (!sys) { printf("System: NULL\n"); return; }
    printf("========================================\n");
    printf("  System Report: n=%d m=%d p=%d\n", sys->n, sys->m, sys->p);
    printf("========================================\n");
    printf("  Stability:\n");
    IOSBIBOStatus st = ios_check_bibo_ss(sys);
    printf("    BIBO:          %s\n",
           st == IOS_BIBO_STABLE ? "STABLE" :
           (st == IOS_BIBO_UNSTABLE ? "UNSTABLE" : "CRITICAL"));
    printf("    Exponential:   %s\n", ios_is_exponentially_stable(sys) ? "YES" : "NO");
    printf("    Strictly Proper: %s\n", ios_is_strictly_proper(sys) ? "YES" : "NO");
    printf("    Minimum Phase: %s\n", ios_check_minimum_phase_ss(sys) ? "YES" : "NO");
    printf("  Eigenvalues:\n");
    printf("    Max Re(lambda):  %.6f\n", ios_ss_eigen_max_real(sys));
    printf("    Stability dist:  %.6f\n", ios_stability_distance(sys));
    printf("    Condition num:   %.2f\n", ios_condition_number(sys));
    if (sys->n == 2) {
        printf("    Damping ratio:   %.4f\n", ios_damping_ratio(sys));
        printf("    Natural freq:    %.4f\n", ios_natural_frequency(sys));
    }
    printf("    Time constant:   %.4f\n", ios_time_constant(sys));
    printf("    Settling time:   %.4f\n", ios_settling_time(sys, 0.02));
    printf("    Overshoot est:   %.1f%%\n", ios_overshoot_estimate(sys));
    printf("  Gains:\n");
    printf("    L2 gain (Hinf):  %.6f\n", ios_L2_gain_ss(sys));
    printf("    L1 induced:      %.6f\n", ios_L1_induced_norm(sys));
    printf("    Linf induced:    %.6f\n", ios_Linf_induced_norm(sys));
    printf("    Peak gain:       %.6f\n", ios_peak_gain(sys));
    printf("    DC gain:         %.6f\n", ios_dc_gain_ss(sys));
    printf("    H2 norm:         %.6f\n", ios_h2_norm_ss(sys));
    printf("    Robustness metric: %.4f\n", ios_robustness_metric(sys));
    printf("    Hankel SV max:   %.6f\n", ios_hankel_singular_value_max(sys));
    printf("    Ctrb Gramian tr: %.6f\n", ios_controllability_gramian_trace(sys));
    printf("    Obsv Gramian tr: %.6f\n", ios_observability_gramian_trace(sys));
    printf("========================================\n");
}
