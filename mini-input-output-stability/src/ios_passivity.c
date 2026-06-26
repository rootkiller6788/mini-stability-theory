#include "ios_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================
 * Passivity Verification
 *
 * A system is passive if (u, y) >= 0 for all u (net energy absorbed).
 * Strict passivity requires (u, y) >= epsilon * ||u||^2 + delta * ||y||^2.
 *
 * For LTI system (A,B,C,D):
 * - Positive Real Lemma: if there exists P = P^T > 0 such that
 *     [A^T*P + P*A   P*B - C^T]
 *     [B^T*P - C     -(D + D^T)]  <= 0
 *   then the system is passive.
 *
 * Simplified tests based on transfer function:
 * - G(s) is PR if Re[G(jw)] >= 0 for all w
 * ============================================================ */

IOSPassivityResult ios_check_passivity_ss(const IOSStateSpace* sys) {
    IOSPassivityResult r = {false, false, false, 0.0};
    if (!sys || sys->n < 1 || sys->m != sys->p) return r;

    /* Check stability first -- unstable systems cannot be passive */
    IOSBIBOStatus st = ios_check_bibo_ss(sys);
    if (st != IOS_BIBO_STABLE) return r;

    /* For SISO: check Re[G(jw)] >= 0 */
    bool all_pr = true;
    double min_re = IOS_INF;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double Re = 0.0;
        if (sys->n == 1) {
            double a = sys->A[0], b = sys->B[0], c = sys->C[0], d = sys->D[0];
            double denom = a * a + w * w;
            Re = d + c * b * a / fmax(denom, IOS_EPS);
        } else {
            /* Diagonal approximation */
            Re = sys->D[0];
            for (int i = 0; i < sys->n; i++) {
                double a_ii = sys->A[i * sys->n + i];
                Re += sys->C[i] * sys->B[i] * a_ii
                      / fmax(a_ii * a_ii + w * w, IOS_EPS);
            }
        }
        if (Re < min_re) min_re = Re;
        if (Re < -IOS_EPS) all_pr = false;
    }

    r.passive = all_pr;
    r.dissipation = (min_re < 0.0) ? 0.0 : min_re;
    r.strictly_input = (min_re > 0.1);
    r.strictly_output = (min_re > 0.01);
    return r;
}

IOSPassivityResult ios_check_passivity_tf(const IOSTransferFcn* tf) {
    IOSPassivityResult r = {false, false, false, 0.0};
    if (!tf) return r;

    /* Check if Re[G(jw)] >= 0 for all w */
    bool all_pr = true;
    double min_re = IOS_INF;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        double G_re = mag * cos(phase);
        if (G_re < min_re) min_re = G_re;
        if (G_re < -IOS_EPS) all_pr = false;
    }

    r.passive = all_pr;
    r.dissipation = (min_re < 0.0) ? 0.0 : min_re;
    r.strictly_input = (min_re > 0.1);
    r.strictly_output = (min_re > 0.01);
    return r;
}

bool ios_is_passive(const IOSSystem* sys) {
    if (!sys || !sys->sys) return false;
    return ios_check_passivity_ss(sys->sys).passive;
}

bool ios_is_strictly_passive(const IOSSystem* sys, double epsilon) {
    if (!sys || !sys->sys) return false;
    IOSPassivityResult r = ios_check_passivity_ss(sys->sys);
    return r.passive && r.dissipation > epsilon;
}

/* ============================================================
 * Popov Criterion
 *
 * For Lur'e system with nonlinearity phi in [0, k]:
 * Stable if there exists eta such that Re[(1 + j*w*eta) * G(jw)] + 1/k > 0.
 * ============================================================ */

bool ios_popov_criterion(const IOSStateSpace* sys, double k) {
    if (!sys || sys->n < 1 || k <= IOS_EPS) return false;

    for (int i = 0; i < IOS_FREQ_POINTS; i++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)i / (double)IOS_FREQ_POINTS);
        double Re = 0.0;

        if (sys->n == 1) {
            double a = sys->A[0], b = sys->B[0], c = sys->C[0], d = sys->D[0];
            double denom = a * a + w * w;
            if (denom < IOS_EPS) continue;
            Re = d + c * b * a / denom;
        } else if (sys->n == 2) {
            double a11 = sys->A[0], a12 = sys->A[1];
            double a21 = sys->A[2], a22 = sys->A[3];
            double b1 = sys->B[0], b2 = sys->B[1];
            double c1 = sys->C[0], c2 = sys->C[1];
            double mat_det = a11 * a22 - a12 * a21 - w * w;
            double mat_tr_w = (a11 + a22) * w;
            double denom = mat_det * mat_det + mat_tr_w * mat_tr_w;
            if (denom < IOS_EPS) continue;
            double num_re = c1 * (b1 * a22 - b2 * a12) + c2 * (-b1 * a21 + b2 * a11);
            double num_im = -w * (c1 * b1 + c2 * b2);
            Re = (num_re * mat_det + num_im * mat_tr_w) / denom + sys->D[0];
        } else {
            Re = sys->D[0];
            for (int j = 0; j < sys->n; j++) {
                double a_jj = sys->A[j * sys->n + j];
                double denom_j = a_jj * a_jj + w * w;
                if (denom_j > IOS_EPS) {
                    Re += sys->C[j] * sys->B[j] * a_jj / denom_j;
                }
            }
        }

        /* Popov: Re[(1 + j*w*eta)*G(jw)] + 1/k > 0 */
        if (Re + 1.0 / k < -IOS_EPS) return false;
    }
    return true;
}

bool ios_popov_criterion_tf(const IOSTransferFcn* tf, double k) {
    if (!tf) return false;
    IOSStateSpace* ss = ios_ss_create(tf->n_den, 1, 1);
    if (!ss) return false;
    /* Companion form */
    for (int i = 0; i < ss->n - 1; i++)
        ss->A[i * ss->n + i + 1] = 1.0;
    for (int i = 0; i < ss->n; i++)
        ss->A[(ss->n - 1) * ss->n + i] = -tf->den[tf->n_den - 1 - i];
    ss->B[ss->n - 1] = 1.0;
    for (int i = 0; i < tf->n_num && i < ss->n; i++)
        ss->C[i] = tf->num[tf->n_num - 1 - i];
    bool result = ios_popov_criterion(ss, k);
    ios_ss_free(ss);
    return result;
}

double ios_popov_sector_bound(const IOSStateSpace* sys) {
    if (!sys) return IOS_INF;
    /* Binary search for maximum k satisfying Popov */
    double lo = 0.01, hi = 100.0;
    for (int i = 0; i < 50; i++) {
        double mid = (lo + hi) / 2.0;
        if (ios_popov_criterion(sys, mid)) lo = mid;
        else hi = mid;
    }
    return lo;
}

/* ============================================================
 * Dissipation Inequality
 *
 * V(x(t1)) - V(x(t0)) <= integral_0^t w(u(t), y(t)) dt
 * ============================================================ */

double ios_dissipation_inequality(const IOSStateSpace* sys, const IOSSignal* u, const IOSSignal* y) {
    if (!sys || !u || !y) return IOS_INF;
    double s = 0.0;
    int n = u->length < y->length ? u->length : y->length;
    for (int i = 0; i < n; i++)
        s += u->data[i] * y->data[i];
    return s;
}

double ios_storage_function(const IOSStateSpace* sys, const double* x) {
    if (!sys || !x) return 0.0;
    /* Quadratic storage: V(x) = 0.5 * x^T * P * x, default P = I */
    double V = 0.0;
    for (int i = 0; i < sys->n; i++)
        V += x[i] * x[i];
    return 0.5 * V;
}

bool ios_is_dissipative(const IOSStateSpace* sys, double supply_rate) {
    if (!sys) return false;
    return ios_check_passivity_ss(sys).dissipation >= supply_rate;
}

double ios_lyapunov_derivative_ss(const IOSStateSpace* sys, const double* x) {
    if (!sys || !x) return 0.0;
    double Vdot = 0.0;
    for (int i = 0; i < sys->n; i++)
        for (int j = 0; j < sys->n; j++)
            Vdot += x[i] * sys->A[i * sys->n + j] * x[j];
    return Vdot;
}

/* ============================================================
 * Conic Sector Computation
 *
 * A system G is inside sector [a, b] if a * ||u||^2 <= (u, Gu) <= b * ||u||^2.
 * For LTI, this corresponds to Nyquist plot of G(jw) lying within cone [a, b].
 * ============================================================ */

IOSSector ios_compute_sector_ss(const IOSStateSpace* sys) {
    IOSSector sec = {0.0, 0.0, IOS_SECTOR_NONE};
    if (!sys) return sec;

    double gmin = IOS_INF, gmax = -IOS_INF;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double G = 0.0;
        /* DC-like gain computed from diagonal of A */
        if (sys->n == 1) {
            G = sys->C[0] * sys->B[0] / fmax(fabs(sys->A[0]), IOS_EPS);
        } else {
            G = 0.0;
            for (int i = 0; i < sys->p; i++)
                for (int j = 0; j < sys->m; j++)
                    G += fabs(sys->D[i * sys->m + j]);
        }
        if (G < gmin) gmin = G;
        if (G > gmax) gmax = G;
    }
    sec.a = gmin;
    sec.b = gmax;
    sec.type = (gmin > -1e90 && gmax < 1e90) ? IOS_SECTOR_FINITE : IOS_SECTOR_INFINITE;
    return sec;
}

IOSSector ios_compute_sector_tf(const IOSTransferFcn* tf) {
    IOSSector sec = {0.0, 0.0, IOS_SECTOR_NONE};
    if (!tf) return sec;
    /* For TF: compute frequency-domain bounds */
    double gmin = IOS_INF, gmax = -IOS_INF;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        double G = mag * cos(phase);
        if (G < gmin) gmin = G;
        if (G > gmax) gmax = G;
    }
    sec.a = gmin;
    sec.b = gmax;
    sec.type = IOS_SECTOR_FINITE;
    return sec;
}

bool ios_in_conic_sector(const IOSTransferFcn* tf, double a, double b) {
    if (!tf) return false;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        double G = mag * cos(phase);
        if (G < a - IOS_EPS || G > b + IOS_EPS) return false;
    }
    return true;
}

/* ============================================================
 * Positive Real Condition
 * ============================================================ */

double ios_positive_real_condition(const IOSTransferFcn* tf) {
    if (!tf) return IOS_INF;
    double mn = IOS_INF;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        double G_re = mag * cos(phase);
        if (G_re < mn) mn = G_re;
    }
    return mn;
}

bool ios_strict_positive_real(const IOSTransferFcn* tf, double eps) {
    return ios_positive_real_condition(tf) >= eps;
}

double ios_passivity_index(const IOSTransferFcn* tf) {
    return ios_positive_real_condition(tf);
}

/* ============================================================
 * Passivity Metrics: Excess and Shortage
 * ============================================================ */

double ios_excess_passivity(const IOSStateSpace* sys) {
    if (!sys) return 0.0;
    return ios_check_passivity_ss(sys).dissipation;
}

double ios_shortage_passivity(const IOSStateSpace* sys) {
    if (!sys) return IOS_INF;
    double d = ios_excess_passivity(sys);
    return d < 0.0 ? -d : 0.0;
}

bool ios_is_very_strictly_passive(const IOSStateSpace* sys) {
    return ios_excess_passivity(sys) > 0.1;
}

bool ios_output_strictly_passive(const IOSStateSpace* sys) {
    if (!sys) return false;
    return ios_check_passivity_ss(sys).strictly_output;
}

/* ============================================================
 * Feedback Passivity Stability
 *
 * Negative feedback interconnection of two passive systems is stable.
 * If either is strictly passive, asymptotic stability follows.
 * ============================================================ */

bool ios_feedback_passivity_stable(const IOSStateSpace* G, const IOSStateSpace* H) {
    if (!G || !H) return false;
    IOSPassivityResult rG = ios_check_passivity_ss(G);
    IOSPassivityResult rH = ios_check_passivity_ss(H);
    return rG.passive && rH.passive;
}

double ios_passivity_gap(const IOSStateSpace* sys) {
    IOSPassivityResult r = ios_check_passivity_ss(sys);
    return r.passive ? 0.0 : 1.0 - r.dissipation;
}

/* ============================================================
 * Mixed Passivity Index
 * ============================================================ */

double ios_mixed_passivity_index(const IOSStateSpace* sys, double alpha) {
    if (!sys) return 0.0;
    double gain = ios_L2_gain_ss(sys);
    double pr_cond = fmax(ios_excess_passivity(sys), 0.0);
    return alpha * pr_cond + (1.0 - alpha) * (gain > IOS_EPS ? 1.0 / gain : 0.0);
}

double ios_passivity_input_index(const IOSStateSpace* sys) {
    IOSPassivityResult r = ios_check_passivity_ss(sys);
    return r.strictly_input ? 1.0 : 0.0;
}

double ios_passivity_output_index(const IOSStateSpace* sys) {
    IOSPassivityResult r = ios_check_passivity_ss(sys);
    return r.strictly_output ? 1.0 : 0.0;
}

/* ============================================================
 * Feedforward Passivity Index
 * ============================================================ */

double ios_feedforward_passivity_index(const IOSTransferFcn* tf) {
    if (!tf) return -IOS_INF;
    return ios_positive_real_condition(tf);
}

/* ============================================================
 * Supply Rate
 * ============================================================ */

double ios_supply_rate_quadratic_form(const IOSStateSpace* sys, const double* x, const double* u) {
    if (!sys || !x || !u) return 0.0;
    double w = 0.0;
    for (int i = 0; i < sys->p; i++) {
        double y = 0.0;
        for (int j = 0; j < sys->n; j++)
            y += sys->C[i * sys->n + j] * x[j];
        w += u[i] * y;
    }
    return w;
}

double ios_dissipation_rate(const IOSStateSpace* sys, const double* x, const double* u) {
    return ios_supply_rate_quadratic_form(sys, x, u);
}

/* ============================================================
 * Positive Real Lemma Check
 * ============================================================ */

bool ios_positive_real_lemma_check(const IOSStateSpace* sys) {
    if (!sys) return false;
    return ios_check_passivity_ss(sys).passive;
}

int ios_passivity_compare(const IOSStateSpace* a, const IOSStateSpace* b) {
    double pa = ios_excess_passivity(a);
    double pb = ios_excess_passivity(b);
    return pa > pb ? 1 : (pa < pb ? -1 : 0);
}

/* ============================================================
 * Reporting Functions
 * ============================================================ */

void ios_passivity_report(const IOSStateSpace* sys) {
    IOSPassivityResult r = ios_check_passivity_ss(sys);
    printf("Passivity: %s  diss=%.4f  strict_in=%s  strict_out=%s\n",
           r.passive ? "PASSIVE" : "NOT PASSIVE",
           r.dissipation,
           r.strictly_input ? "YES" : "NO",
           r.strictly_output ? "YES" : "NO");
}

void ios_passivity_full_report(const IOSStateSpace* sys) {
    ios_passivity_report(sys);
}

void ios_passivity_margin_report_full(const IOSStateSpace* sys) {
    ios_passivity_report(sys);
}

void ios_passivity_connect_report(const IOSStateSpace* G, const IOSStateSpace* H) {
    bool ok = ios_feedback_passivity_stable(G, H);
    printf("Passivity Feedback: ");
    printf("G=%s H=%s -> %s\n",
           ios_check_passivity_ss(G).passive ? "PASSIVE" : "NOT",
           ios_check_passivity_ss(H).passive ? "PASSIVE" : "NOT",
           ok ? "STABLE" : "INCONCLUSIVE");
}

void ios_sector_report(const IOSTransferFcn* tf) {
    IOSSector s = ios_compute_sector_tf(tf);
    printf("Sector: [%.4f, %.4f] type=%s\n",
           s.a, s.b,
           s.type == IOS_SECTOR_FINITE ? "FINITE" :
           (s.type == IOS_SECTOR_INFINITE ? "INFINITE" : "NONE"));
}

void ios_sector_analysis_report(const IOSTransferFcn* tf) {
    ios_sector_report(tf);
}

void ios_popov_sweep(const IOSTransferFcn* tf) {
    if (!tf) return;
    printf("=== Popov Criterion Sweep ===\n");
    for (double k = 0.1; k <= 10.0; k *= 2.0) {
        bool stable = ios_popov_criterion_tf(tf, k);
        printf("  k=%.1f: %s\n", k, stable ? "STABLE" : "INCONCLUSIVE");
    }
}

void ios_passivity_all_metrics(const IOSStateSpace* sys) {
    if (!sys) return;
    IOSPassivityResult r = ios_check_passivity_ss(sys);
    printf("Passivity Metrics: passive=%s strict_in=%s strict_out=%s excess=%.4f shortage=%.4f\n",
           r.passive ? "YES" : "NO",
           r.strictly_input ? "YES" : "NO",
           r.strictly_output ? "YES" : "NO",
           ios_excess_passivity(sys), ios_shortage_passivity(sys));
}
