#include "ios_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================
 * Small-Gain Theorem (Zames, 1966)
 *
 * For two interconnected systems H1, H2 with L2 gains gamma1, gamma2:
 *   gamma1 * gamma2 < 1  =>  feedback interconnection is L2 stable
 * ============================================================ */

IOSSmallGainResult ios_small_gain_test(double gamma1, double gamma2) {
    IOSSmallGainResult r;
    r.gamma = gamma1 * gamma2;
    r.satisfied = (r.gamma < 1.0);
    r.margin = 1.0 - r.gamma;
    return r;
}

bool ios_small_gain_stable(const IOSSystem* sys1, const IOSSystem* sys2) {
    if (!sys1 || !sys2) return false;
    double g1 = ios_L2_gain_ss(sys1->sys);
    double g2 = ios_L2_gain_ss(sys2->sys);
    return (g1 * g2 < 1.0);
}

bool ios_scaled_small_gain(const IOSSystem* sys1, const IOSSystem* sys2, double scale) {
    if (!sys1 || !sys2 || scale <= 0.0) return false;
    double g1 = ios_L2_gain_ss(sys1->sys) * scale;
    double g2 = ios_L2_gain_ss(sys2->sys) / scale;
    return (g1 * g2 < 1.0);
}

/* ============================================================
 * Circle Criterion (Lur'e problem)
 *
 * Feedback: y = G*u, u = -phi(y) with phi in sector [k1, k2].
 * Stable if Nyquist plot of G(jw) does not enter the critical disk
 * centered at -(1/k1 + 1/k2)/2 with radius (1/k1 - 1/k2)/2.
 * ============================================================ */

bool ios_circle_criterion(const IOSStateSpace* sys, double k1, double k2) {
    if (!sys || sys->n < 1 || fabs(k2 - k1) < IOS_EPS) return false;

    double center = -(1.0 / k1 + 1.0 / k2) / 2.0;
    double radius = fabs(1.0 / k1 - 1.0 / k2) / 2.0;

    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double Re = 0.0, Im = 0.0;

        /* Compute G(jw) from state-space */
        if (sys->n == 1) {
            double a = sys->A[0], b = sys->B[0], c = sys->C[0];
            double denom = a * a + w * w;
            if (denom < IOS_EPS) continue;
            Re = c * b * a / denom + sys->D[0];
            Im = -c * b * w / denom;
        } else if (sys->n == 2) {
            double a11 = sys->A[0], a12 = sys->A[1];
            double a21 = sys->A[2], a22 = sys->A[3];
            double b1 = sys->B[0], b2 = sys->B[1];
            double c1 = sys->C[0], c2 = sys->C[1];
            double d = sys->D[0];
            /* det(jw*I - A) */
            double mat_det = (a11 * a22 - a12 * a21 - w * w);
            double mat_tr_w = (a11 + a22) * w;
            double denom = mat_det * mat_det + mat_tr_w * mat_tr_w;
            if (denom < IOS_EPS) continue;
            /* G(jw) = C*(jwI-A)^{-1}*B + D */
            double num_re = c1 * (b1 * a22 - b2 * a12) + c2 * (-b1 * a21 + b2 * a11);
            double num_im = -w * (c1 * b1 + c2 * b2);
            Re = (num_re * mat_det + num_im * mat_tr_w) / denom + d;
            Im = (num_im * mat_det - num_re * mat_tr_w) / denom;
        } else {
            /* For higher order, approximate using pole/w decomposition */
            Re = sys->D[0];
            for (int i = 0; i < sys->n; i++) {
                double a_ii = sys->A[i * sys->n + i];
                double denom_i = a_ii * a_ii + w * w;
                if (denom_i > IOS_EPS) {
                    Re += sys->C[i] * sys->B[i] * a_ii / denom_i;
                    Im -= sys->C[i] * sys->B[i] * w / denom_i;
                }
            }
        }

        double dist = sqrt((Re - center) * (Re - center) + Im * Im);
        if (dist < radius - IOS_EPS) return false;
    }
    return true;
}

bool ios_circle_criterion_tf(const IOSTransferFcn* tf, double k1, double k2) {
    if (!tf) return false;
    /* Convert TF to SS and use circle criterion */
    IOSStateSpace* ss = ios_ss_create(tf->n_den, 1, 1);
    if (!ss) return false;
    /* Companion form: A = [[0 I]; [-a_n ... -a_1]], B = [0 ... 1]^T, C = [b_n ... b_1] */
    for (int i = 0; i < ss->n - 1; i++)
        ss->A[i * ss->n + i + 1] = 1.0;
    for (int i = 0; i < ss->n; i++)
        ss->A[(ss->n - 1) * ss->n + i] = -tf->den[tf->n_den - 1 - i];
    ss->B[ss->n - 1] = 1.0;
    for (int i = 0; i < tf->n_num && i < ss->n; i++)
        ss->C[i] = tf->num[tf->n_num - 1 - i];
    bool result = ios_circle_criterion(ss, k1, k2);
    ios_ss_free(ss);
    return result;
}

IOSCircleResult ios_circle_analysis(const IOSTransferFcn* tf, double k1, double k2) {
    IOSCircleResult r;
    r.center = -(1.0 / k1 + 1.0 / k2) / 2.0;
    r.radius = fabs(1.0 / k1 - 1.0 / k2) / 2.0;
    r.stable = ios_circle_criterion_tf(tf, k1, k2);
    return r;
}

/* ============================================================
 * Nyquist Distance / Disk Margin
 * ============================================================ */

double ios_nyquist_distance(const IOSTransferFcn* tf, double w) {
    if (!tf) return IOS_INF;
    double mag = 0.0, phase = 0.0;
    ios_frequency_response(tf, w, &mag, &phase);
    double G_re = mag * cos(phase);
    double G_im = mag * sin(phase);
    return sqrt((1.0 + G_re) * (1.0 + G_re) + G_im * G_im);
}

double ios_disk_margin(const IOSTransferFcn* tf) {
    if (!tf) return 0.0;
    double min_dist = IOS_INF;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double d = ios_nyquist_distance(tf, w);
        if (d < min_dist) min_dist = d;
    }
    return min_dist;
}

double ios_modulus_margin(const IOSTransferFcn* tf) {
    return ios_disk_margin(tf);
}

/* ============================================================
 * Sensitivity Functions
 *
 * S(s) = 1 / (1 + L(s))  -- sensitivity
 * T(s) = L(s) / (1 + L(s)) -- complementary sensitivity
 * ============================================================ */

double ios_sensitivity_peak(const IOSTransferFcn* L) {
    if (!L) return IOS_INF;
    double maxS = 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(L, w, &mag, &phase);
        double G_re = mag * cos(phase);
        double G_im = mag * sin(phase);
        double S = 1.0 / sqrt((1.0 + G_re) * (1.0 + G_re) + G_im * G_im);
        if (S > maxS) maxS = S;
    }
    return maxS;
}

double ios_complementary_sensitivity(const IOSTransferFcn* L) {
    if (!L) return IOS_INF;
    double maxT = 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(L, w, &mag, &phase);
        double G_re = mag * cos(phase);
        double G_im = mag * sin(phase);
        double denom = (1.0 + G_re) * (1.0 + G_re) + G_im * G_im;
        double T = denom > IOS_EPS ? sqrt((G_re * G_re + G_im * G_im) / denom) : IOS_INF;
        if (T > maxT) maxT = T;
    }
    return maxT;
}

double ios_sensitivity_peak_tf(const IOSTransferFcn* tf) {
    return ios_sensitivity_peak(tf);
}

double ios_complementary_sensitivity_tf(const IOSTransferFcn* tf) {
    return ios_complementary_sensitivity(tf);
}

/* ============================================================
 * Gap Metric (Georgiou & Smith, 1990)
 *
 * gap(G1, G2) = sup_w kappa(G1(jw), G2(jw))
 * where kappa is the chordal distance between points on Riemann sphere
 * ============================================================ */

double ios_gap_metric(const IOSTransferFcn* G1, const IOSTransferFcn* G2) {
    if (!G1 || !G2) return IOS_INF;
    double max_gap = 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag1 = 0.0, phase1 = 0.0;
        double mag2 = 0.0, phase2 = 0.0;
        ios_frequency_response(G1, w, &mag1, &phase1);
        ios_frequency_response(G2, w, &mag2, &phase2);
        double G1_re = mag1 * cos(phase1), G1_im = mag1 * sin(phase1);
        double G2_re = mag2 * cos(phase2), G2_im = mag2 * sin(phase2);
        /* Chordal distance: |G1 - G2| / sqrt((1+|G1|^2)(1+|G2|^2)) */
        double diff_sq = (G1_re - G2_re) * (G1_re - G2_re) + (G1_im - G2_im) * (G1_im - G2_im);
        double norm1_sq = 1.0 + G1_re * G1_re + G1_im * G1_im;
        double norm2_sq = 1.0 + G2_re * G2_re + G2_im * G2_im;
        double gap = norm1_sq * norm2_sq > IOS_EPS ?
                     sqrt(diff_sq / (norm1_sq * norm2_sq)) : 0.0;
        if (gap > max_gap) max_gap = gap;
    }
    return max_gap;
}

/* ============================================================
 * Structured Singular Value (mu) approximation
 * ============================================================ */

double ios_structured_singular_value(const IOSStateSpace* sys) {
    if (!sys) return IOS_INF;
    /* Simple bound: mu <= max singular value */
    return ios_L2_gain_ss(sys);
}

/* ============================================================
 * Incremental Gain (for nonlinear systems)
 * ============================================================ */

double ios_incremental_gain(const IOSMapFunction f, void* params, int n) {
    if (!f || n < 1) return IOS_INF;
    double* u1 = (double*)calloc((size_t)n, sizeof(double));
    double* u2 = (double*)calloc((size_t)n, sizeof(double));
    double* y1 = (double*)calloc((size_t)n, sizeof(double));
    double* y2 = (double*)calloc((size_t)n, sizeof(double));
    if (!u1 || !u2 || !y1 || !y2) {
        free(u1); free(u2); free(y1); free(y2);
        return IOS_INF;
    }
    u2[0] = 1.0;  /* small perturbation */
    f(u1, y1, n, params);
    f(u2, y2, n, params);
    double du = 0.0, dy = 0.0;
    for (int i = 0; i < n; i++) {
        du += (u1[i] - u2[i]) * (u1[i] - u2[i]);
        dy += (y1[i] - y2[i]) * (y1[i] - y2[i]);
    }
    free(u1); free(u2); free(y1); free(y2);
    return (du > IOS_EPS) ? sqrt(dy / du) : IOS_INF;
}

/* ============================================================
 * Zames-Falb Criterion
 * (Multiplier-based stability for monotone nonlinearities)
 * ============================================================ */

bool ios_zames_falb_criterion(const IOSTransferFcn* tf, double k) {
    if (!tf || k <= 0.0) return false;
    /* Zames-Falb: If there exists a multiplier Z(s) such that
       Re[Z(jw)*G(jw)] + 1/k > 0 for all w, system is stable.
       Here we check the simplest case: Z(s) = 1 (circle criterion). */
    double min_re = IOS_INF;
    for (int j = 0; j < IOS_FREQ_POINTS; j++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)j / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        double G_re = mag * cos(phase);
        if (G_re + 1.0 / k < min_re)
            min_re = G_re + 1.0 / k;
    }
    return min_re > IOS_EPS;
}

/* ============================================================
 * Nichols Peak
 * ============================================================ */

double ios_nichols_peak(const IOSTransferFcn* tf) {
    if (!tf) return 0.0;
    double max_mag = 0.0;
    for (int k = 0; k < IOS_FREQ_POINTS; k++) {
        double w = 0.01 * pow(10.0, 4.0 * (double)k / (double)IOS_FREQ_POINTS);
        double mag = 0.0, phase = 0.0;
        ios_frequency_response(tf, w, &mag, &phase);
        double G_re = mag * cos(phase);
        double G_im = mag * sin(phase);
        double denom = (1.0 + G_re) * (1.0 + G_re) + G_im * G_im;
        double T = denom > IOS_EPS ? mag / sqrt(denom) : IOS_INF;
        if (T > max_mag) max_mag = T;
    }
    return max_mag;
}

/* ============================================================
 * Robust Stability Margin
 * ============================================================ */

double ios_robust_stability_margin(const IOSTransferFcn* G) {
    return G ? ios_disk_margin(G) : 0.0;
}

/* ============================================================
 * Multi-loop Small Gain
 * ============================================================ */

bool ios_multiloop_small_gain_test(const IOSTransferFcn** G, int n, double* gammas) {
    if (!G || !gammas || n < 1) return false;
    double prod = 1.0;
    for (int i = 0; i < n; i++) {
        gammas[i] = G[i] ? ios_Hinf_norm_tf(G[i]) : IOS_INF;
        prod *= gammas[i];
    }
    return prod < 1.0;
}

/* ============================================================
 * Reporting Functions
 * ============================================================ */

void ios_sg_report(double g1, double g2) {
    IOSSmallGainResult r = ios_small_gain_test(g1, g2);
    printf("SG: %.4f * %.4f = %.4f  %s  margin=%.4f\n",
           g1, g2, r.gamma, r.satisfied ? "OK" : "FAIL", r.margin);
}

void ios_smallgain_report_full(double g1, double g2) {
    printf("=== Small-Gain Analysis ===\n");
    ios_sg_report(g1, g2);
}

void ios_circle_report(const IOSTransferFcn* tf, double k1, double k2) {
    IOSCircleResult r = ios_circle_analysis(tf, k1, k2);
    printf("Circle Criterion [%.2f, %.2f]: center=%.4f radius=%.4f stable=%s\n",
           k1, k2, r.center, r.radius, r.stable ? "YES" : "NO");
}

void ios_circle_criterion_sweep(const IOSTransferFcn* tf) {
    if (!tf) return;
    printf("=== Circle Criterion Sweep ===\n");
    double k_pairs[][2] = {{0.1, 1.0}, {0.0, 5.0}, {0.0, 10.0}, {0.5, 2.0}, {1.0, 5.0}};
    for (int i = 0; i < 5; i++) {
        bool stable = ios_circle_criterion_tf(tf, k_pairs[i][0], k_pairs[i][1]);
        printf("  [%.1f, %.1f]: %s\n", k_pairs[i][0], k_pairs[i][1],
               stable ? "STABLE" : "INCONCLUSIVE");
    }
}

void ios_smallgain_parameter_sweep(double g1_min, double g1_max, int n1,
                                     double g2_min, double g2_max, int n2) {
    printf("=== Small-Gain Parameter Sweep ===\n");
    double dg1 = (n1 > 1) ? (g1_max - g1_min) / (double)(n1 - 1) : 0.0;
    double dg2 = (n2 > 1) ? (g2_max - g2_min) / (double)(n2 - 1) : 0.0;
    for (int i = 0; i < n1; i++) {
        double g1 = g1_min + (double)i * dg1;
        for (int j = 0; j < n2; j++) {
            double g2 = g2_min + (double)j * dg2;
            IOSSmallGainResult r = ios_small_gain_test(g1, g2);
            if (r.satisfied)
                printf("  Stable: g1=%.2f g2=%.2f product=%.4f margin=%.4f\n",
                       g1, g2, r.gamma, r.margin);
        }
    }
}

void ios_smallgain_sweep(double* g1_vals, int n1, double* g2_vals, int n2, double** results) {
    if (!g1_vals || !g2_vals || !results) return;
    for (int i = 0; i < n1; i++)
        for (int j = 0; j < n2; j++)
            results[i][j] = (g1_vals[i] * g2_vals[j] < 1.0) ? 0.0 : 1.0;
}

void ios_smallgain_print_summary(double g1, double g2) {
    ios_sg_report(g1, g2);
}

void ios_robustness_report(const IOSTransferFcn* tf) {
    if (!tf) return;
    printf("=== Robustness Report ===\n");
    printf("  Disk margin:    %.4f\n", ios_disk_margin(tf));
    printf("  Modulus margin: %.4f\n", ios_modulus_margin(tf));
    printf("  Nichols peak:   %.4f\n", ios_nichols_peak(tf));
}

double ios_gap_metric_L2(const IOSTransferFcn* G1, const IOSTransferFcn* G2) {
    return ios_gap_metric(G1, G2);
}

/* ============================================================
 * Mu-Analysis (Structured Singular Value) Approximations
 *
 * mu(M) = 1 / min{sigma_bar(Delta) : det(I - M*Delta) = 0}
 * Upper bound: mu(M) <= sigma_max(M)
 * Lower bound (for 1-block): mu(M) >= rho(M) (spectral radius)
 * ============================================================ */

/* Upper bound: maximum singular value */
double ios_mu_upper_bound(const IOSStateSpace* sys) {
    return ios_L2_gain_ss(sys);
}

/* Lower bound: spectral radius of (A + B*K*C) for certain K */
double ios_mu_lower_bound(const IOSStateSpace* sys) {
    if (!sys || sys->n < 1) return 0.0;
    double rho = 0.0;
    for (int i = 0; i < sys->n; i++) {
        double a_ii = sys->A[i * sys->n + i];
        double bc_sum = 0.0;
        for (int j = 0; j < sys->m; j++)
            bc_sum += fabs(sys->B[i * sys->m + j]);
        for (int j = 0; j < sys->p; j++)
            bc_sum *= fabs(sys->C[j * sys->n + i]);
        double val = fabs(a_ii) + bc_sum;
        if (val > rho) rho = val;
    }
    return rho;
}

/* Mixed mu bound (D-scale approximation) */
double ios_mu_d_scale(const IOSStateSpace* sys, double d) {
    if (!sys || d <= IOS_EPS) return IOS_INF;
    /* D-scaled mu: inf_D sigma_max(D*M*D^{-1}) */
    double scaled = ios_L2_gain_ss(sys);
    /* Simple heuristic: scale output channels */
    return d > 1.0 ? scaled / d : scaled * d;
}

/* ============================================================
 * IQC (Integral Quadratic Constraint) Analysis
 *
 * Simplest IQC: Re[F(jw)*Delta(jw)] >= 0 for multiplier F(s)
 * ============================================================ */

/* Check if Delta satisfies IQC with multiplier F */
bool ios_iqc_check(const IOSTransferFcn* G, const IOSTransferFcn* Delta, double tau) {
    if (!G || !Delta || tau <= 0.0) return false;
    /* Simple check: ||G|| * ||Delta|| < 1 */
    double gG = ios_Hinf_norm_tf(G);
    double gD = ios_Hinf_norm_tf(Delta);
    IOSSmallGainResult sg = ios_small_gain_test(gG, gD);
    return sg.satisfied;
}

/* ============================================================
 * Robust Performance Analysis
 *
 * Combines robust stability with nominal performance.
 * SSV (structured singular value) framework: mu <= 1 for RP.
 * ============================================================ */

/* Robust performance criterion (simplified) */
bool ios_robust_performance_check(const IOSTransferFcn* G, const IOSTransferFcn* Wp,
                                    const IOSTransferFcn* Wu, double gamma) {
    if (!G || !Wp || !Wu) return false;
    /* ||Wp*S|| + ||Wu*T|| < 1 for robust performance */
    double Ms = ios_sensitivity_peak(G);
    double Mt = ios_complementary_sensitivity(G);
    double wp = ios_Hinf_norm_tf(Wp);
    double wu = ios_Hinf_norm_tf(Wu);
    return (wp * Ms + wu * Mt < gamma);
}

/* Nominal performance check */
bool ios_nominal_performance(const IOSTransferFcn* G, double tol) {
    if (!G) return false;
    double Ms = ios_sensitivity_peak(G);
    return Ms < tol;
}

/* ============================================================
 * Vinnicombe Nu-Gap Metric
 *
 * nu(G1, G2) = sup_w kappa(G1(jw), G2(jw)) if winding number condition holds
 * ============================================================ */

double ios_nu_gap_metric(const IOSTransferFcn* G1, const IOSTransferFcn* G2) {
    /* Nu-gap reduces to gap metric when winding number condition satisfied */
    double gap = ios_gap_metric(G1, G2);
    /* Check winding number condition (simplified): both stable */
    if (ios_tf_is_stable(G1) && ios_tf_is_stable(G2))
        return gap;
    /* If either is unstable, nu-gap may be 1.0 (max) */
    return 1.0;
}

/* ============================================================
 * Frequency-Domain Robustness Analysis
 * ============================================================ */

/* Maximum modulus of sensitivity (Ms) */
double ios_ms_peak(const IOSTransferFcn* G) {
    return ios_sensitivity_peak(G);
}

/* Maximum modulus of complementary sensitivity (Mt) */
double ios_mt_peak(const IOSTransferFcn* G) {
    return ios_complementary_sensitivity(G);
}

/* Robustness index: min(Ms, Mt) */
double ios_robustness_index(const IOSTransferFcn* G) {
    if (!G) return 0.0;
    double Ms = ios_sensitivity_peak(G);
    double Mt = ios_complementary_sensitivity(G);
    return Ms < Mt ? Ms : Mt;
}

/* ============================================================
 * Describing Function Analysis (simplified)
 * ============================================================ */

/* Describing function for saturation nonlinearity */
double ios_saturation_df(double A, double sat_level) {
    if (A <= sat_level + IOS_EPS) return 1.0;
    double ratio = sat_level / A;
    return (2.0 / IOS_PI) * (asin(ratio) + ratio * sqrt(1.0 - ratio * ratio));
}

/* Describing function for deadzone nonlinearity */
double ios_deadzone_df(double A, double deadzone) {
    if (A <= deadzone + IOS_EPS) return 0.0;
    double ratio = deadzone / A;
    return 1.0 - (2.0 / IOS_PI) * (asin(ratio) + ratio * sqrt(1.0 - ratio * ratio));
}

/* Describing function for relay nonlinearity */
double ios_relay_df(double A, double relay_amp) {
    if (A < IOS_EPS) return IOS_INF;
    return 4.0 * relay_amp / (IOS_PI * A);
}

/* ============================================================
 * Comprehensive Small-Gain Report
 * ============================================================ */

void ios_smallgain_comprehensive_report(const IOSTransferFcn* G1, const IOSTransferFcn* G2) {
    if (!G1 || !G2) return;
    printf("========================================\n");
    printf("  Small-Gain Comprehensive Report\n");
    printf("========================================\n");
    double g1 = ios_Hinf_norm_tf(G1);
    double g2 = ios_Hinf_norm_tf(G2);
    IOSSmallGainResult sg = ios_small_gain_test(g1, g2);
    printf("  System 1 Hinf:   %.6f\n", g1);
    printf("  System 2 Hinf:   %.6f\n", g2);
    printf("  Product:         %.6f\n", sg.gamma);
    printf("  Small-Gain:      %s\n", sg.satisfied ? "SATISFIED (stable)" : "NOT SATISFIED");
    printf("  Margin:          %.6f\n", sg.margin);
    printf("  Gap metric:      %.6f\n", ios_gap_metric(G1, G2));
    printf("  Nu-gap:          %.6f\n", ios_nu_gap_metric(G1, G2));
    printf("========================================\n");
}

/* ============================================================
 * Circle Criterion Comprehensive Analysis
 * ============================================================ */

void ios_circle_comprehensive_report(const IOSTransferFcn* G, double k1, double k2) {
    if (!G) return;
    printf("========================================\n");
    printf("  Circle Criterion Comprehensive Report\n");
    printf("========================================\n");
    printf("  Sector bounds:   [%.4f, %.4f]\n", k1, k2);
    IOSCircleResult cr = ios_circle_analysis(G, k1, k2);
    printf("  Critical disk:\n");
    printf("    Center:        %.6f\n", cr.center);
    printf("    Radius:        %.6f\n", cr.radius);
    printf("  Circle criterion: %s\n", cr.stable ? "SATISFIED" : "NOT SATISFIED");
    double dm = ios_disk_margin(G);
    printf("  Disk margin:      %.6f\n", dm);
    double Ms = ios_sensitivity_peak(G);
    printf("  Sensitivity peak: %.6f\n", Ms);
    printf("========================================\n");
}
