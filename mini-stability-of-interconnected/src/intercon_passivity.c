#include "intercon_passivity.h"
#include "intercon_small_gain.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Passivity Test via KYP Lemma
 *
 * LTI system (A,B,C,D) is passive iff there exists P > 0:
 *   [A^T*P + P*A    P*B - C^T  ] <= 0
 *   [B^T*P - C      -(D + D^T)  ]
 *
 * For SISO 1st-order system: passive if D > 0 and C*B > 0
 * ============================================================== */

PassivityResult* ic_passivity_test(const ICSubsystem* s) {
    PassivityResult* pr = calloc(1, sizeof(PassivityResult));
    if (!pr) return NULL;

    /* Check basic necessary conditions */
    if (s->m != s->p) {
        pr->type = PASSIVE_ACTIVE;
        pr->is_passive = false;
        return pr;
    }

    /* For SISO 1D: passive iff D >= 0 and C*B >= 0 (positive real) */
    if (s->n == 1 && s->m == 1) {
        double a = s->A[0][0];
        double b = s->B[0][0];
        double c = s->C[0][0];
        double d = s->D[0][0];

        /* Positive real condition: Re[G(jw)] >= 0 for all w */
        bool pos_real = true;
        for (int i = 0; i < 100; i++) {
            double w = 0.01 * pow(10.0, 4.0 * i / 99.0);
            double mag, phase;
            ic_subsys_freq_response(s, w, &mag, &phase);
            double re_G = mag * cos(phase);
            if (re_G < -IC_EPSILON) { pos_real = false; break; }
        }

        pr->is_passive = pos_real;
        if (pos_real && d > IC_EPSILON) {
            pr->type = PASSIVE_STRICT;
            pr->excess = d;
        } else if (pos_real) {
            pr->type = PASSIVE_LOSSY;
            pr->excess = 0.0;
        } else {
            pr->type = PASSIVE_ACTIVE;
        }
        pr->dissipation = d;
        return pr;
    }

    /* For 2D system: solve KYP via LMI (check with grid over P) */
    if (s->n == 2 && s->m == 1) {
        bool found_passive = false;
        for (int pi = 0; pi < 20; pi++) {
            for (int pj = 0; pj < 20; pj++) {
                double p11 = 0.1 + 9.9 * pi / 19.0;
                double p22 = 0.1 + 9.9 * pj / 19.0;
                double p12 = -2.0 + 4.0 * (pi + pj) / 38.0;

                /* Compute KYP matrix */
                double M[3][3] = {{0}};
                /* M11 = A^T*P + P*A */
                M[0][0] = 2.0*s->A[0][0]*p11 + 2.0*s->A[1][0]*p12;
                M[0][1] = s->A[0][1]*p11 + s->A[0][0]*p12
                        + s->A[1][1]*p12 + s->A[1][0]*p22;
                M[1][0] = M[0][1];
                M[1][1] = 2.0*s->A[0][1]*p12 + 2.0*s->A[1][1]*p22;

                /* M12 = P*B - C^T */
                M[0][2] = p11*s->B[0][0] + p12*s->B[1][0] - s->C[0][0];
                M[1][2] = p12*s->B[0][0] + p22*s->B[1][0] - s->C[0][1];
                M[2][0] = M[0][2];
                M[2][1] = M[1][2];

                /* M22 = -(D + D^T) */
                M[2][2] = -2.0 * s->D[0][0];

                /* Check negative semidefinite: all principal minors <= 0? */
                /* Simplified: check trace <= 0 and det >= 0 (for 3x3) */
                double det = M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1])
                           - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0])
                           + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]);
                if (M[0][0] <= IC_EPSILON && M[0][0]+M[1][1] <= IC_EPSILON
                    && det >= -IC_EPSILON) {
                    found_passive = true;
                    break;
                }
            }
            if (found_passive) break;
        }
        pr->is_passive = found_passive;
        pr->type = found_passive ? PASSIVE_LOSSY : PASSIVE_ACTIVE;
        return pr;
    }

    /* Default: check via frequency response */
    bool all_pos_real = true;
    double min_re = 1e6;
    for (int i = 0; i < 100; i++) {
        double w = 0.01 * pow(10.0, 4.0 * i / 99.0);
        double mag, phase;
        ic_subsys_freq_response(s, w, &mag, &phase);
        double re = mag * cos(phase);
        if (re < min_re) min_re = re;
        if (re < -IC_EPSILON) all_pos_real = false;
    }
    pr->is_passive = all_pos_real;
    pr->type = all_pos_real ? PASSIVE_LOSSY : PASSIVE_ACTIVE;
    pr->excess = min_re;
    pr->shortage = (min_re < 0) ? -min_re : 0.0;
    return pr;
}

void ic_passivity_result_free(PassivityResult* pr) { free(pr); }

void ic_passivity_result_print(const PassivityResult* pr) {
    printf("PassivityResult: passive=%s  type=",
           pr->is_passive ? "YES" : "NO");
    switch (pr->type) {
        case PASSIVE_STRICT:  printf("STRICT"); break;
        case PASSIVE_LOSSY:   printf("LOSSY"); break;
        case PASSIVE_ACTIVE:  printf("ACTIVE"); break;
        default: printf("NONE"); break;
    }
    printf("  excess=%.4f  dissipation=%.4f\n", pr->excess, pr->dissipation);
}

/* ==============================================================
 * Passivity-Based Feedback Stability
 *
 * Theorem: Feedback interconnection of two passive systems
 * is L2-stable if at least one is strictly passive.
 * ============================================================== */

bool ic_passivity_feedback_stable(const ICSubsystem* s1,
                                    const ICSubsystem* s2) {
    PassivityResult* pr1 = ic_passivity_test(s1);
    PassivityResult* pr2 = ic_passivity_test(s2);
    if (!pr1 || !pr2) {
        ic_passivity_result_free(pr1);
        ic_passivity_result_free(pr2);
        return false;
    }
    bool stable = (pr1->is_passive && pr2->is_passive)
               && (pr1->type == PASSIVE_STRICT || pr2->type == PASSIVE_STRICT);
    ic_passivity_result_free(pr1);
    ic_passivity_result_free(pr2);
    return stable;
}

double ic_passivity_index(const ICSubsystem* s) {
    /* Find largest nu such that G(s) + nu is passive */
    double lo = -10.0, hi = 10.0;
    double best = -1e6;
    for (int iter = 0; iter < 30; iter++) {
        double mid = (lo + hi) / 2.0;
        /* Check if G(s) + mid is passive: test Re[G(jw)] + mid >= 0 */
        bool passive = true;
        for (int i = 0; i < 50; i++) {
            double w = 0.01 * pow(10.0, 2.0 * i / 49.0);
            double mag, phase;
            ic_subsys_freq_response(s, w, &mag, &phase);
            double re = mag * cos(phase);
            if (re + mid < -IC_EPSILON) { passive = false; break; }
        }
        if (passive) { best = mid; lo = mid; }
        else hi = mid;
    }
    return best;
}

double ic_output_feedback_passivity(const ICSubsystem* s) {
    /* Find u = -K*y + v making system passive */
    double best_K = 0.0;
    double best_index = -1e6;
    for (int ki = 0; ki < 50; ki++) {
        double K = 0.01 + 9.99 * ki / 49.0;
        /* With output feedback: G_cl(s) = G(s)/(1 + K*G(s)) */
        bool passive = true;
        for (int i = 0; i < 50; i++) {
            double w = 0.01 * pow(10.0, 2.0 * i / 49.0);
            double mag, phase;
            ic_subsys_freq_response(s, w, &mag, &phase);
            /* Simple check: with large K, feedforward term dominates */
            double re = mag * cos(phase);
            if (re + K < -IC_EPSILON) { passive = false; break; }
        }
        if (passive && K > best_K) best_K = K;
    }
    return best_K;
}

/* ==============================================================
 * Storage Function
 * ============================================================== */

void ic_storage_function(const ICSubsystem* s, double* P, int max_dim) {
    int n = s->n;
    if (n > max_dim) return;
    /* Solve KYP for P: A^T*P + P*A = -Q where Q = C^T*C (for passivity) */
    double Q[IC_MAX_DIM * IC_MAX_DIM] = {0};
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < s->p; k++)
                Q[i*n+j] += s->C[k][i] * s->C[k][j];
    ic_solve_lyapunov((double*)s->A, Q, P, n);
}

double ic_available_storage(const ICSubsystem* s) {
    /* Available storage = sup integral of -u^T*y dt
     * For stable system, approx by 0.5*x^T*P*x */
    double P[IC_MAX_DIM * IC_MAX_DIM];
    ic_storage_function(s, P, IC_MAX_DIM);
    double V = 0.0;
    for (int i = 0; i < s->n; i++)
        for (int j = 0; j < s->n; j++)
            V += s->x[i] * P[i*s->n+j] * s->x[j];
    return 0.5 * V;
}

double ic_required_supply(const ICSubsystem* s, const double* x) {
    /* Required supply to reach state x from origin */
    double P[IC_MAX_DIM * IC_MAX_DIM];
    ic_storage_function(s, P, IC_MAX_DIM);
    double V = 0.0;
    for (int i = 0; i < s->n; i++)
        for (int j = 0; j < s->n; j++)
            V += x[i] * P[i*s->n+j] * x[j];
    return 0.5 * V;
}

/* ==============================================================
 * Dissipativity Framework
 * ============================================================== */

SupplyRate* ic_supply_rate_create(int dim, const double* Q, const double* S,
                                    const double* R) {
    SupplyRate* sr = calloc(1, sizeof(SupplyRate));
    if (!sr) return NULL;
    sr->dim = dim;
    if (Q) memcpy(sr->Q, Q, (size_t)(dim*dim) * sizeof(double));
    if (S) memcpy(sr->S, S, (size_t)(dim*dim) * sizeof(double));
    if (R) memcpy(sr->R, R, (size_t)(dim*dim) * sizeof(double));
    return sr;
}

void ic_supply_rate_free(SupplyRate* sr) { free(sr); }

DissipativityResult* ic_dissipativity_test(const ICSubsystem* s,
                                             const SupplyRate* sr) {
    DissipativityResult* dr = calloc(1, sizeof(DissipativityResult));
    if (!dr) return NULL;
    dr->n = s->n;

    /* For SISO 1D: dissipativity means the system satisfies
     * the frequency-domain inequality. Check via Nyquist. */
    bool dissipative = true;
    for (int i = 0; i < 50; i++) {
        double w = 0.01 * pow(10.0, 3.0 * i / 49.0);
        double mag, phase;
        ic_subsys_freq_response(s, w, &mag, &phase);
        double re = mag * cos(phase), im = mag * sin(phase);
        /* Supply rate: w(u,y) = Q*y^2 + 2*S*y*u + R*u^2
         * For u = 1, y = G(jw): w = Q*|G|^2 + 2*S*Re(G) + R >= 0 */
        double w_val = sr->Q[0][0] * mag * mag
                     + 2.0 * sr->S[0][0] * re
                     + sr->R[0][0];
        if (w_val < -IC_EPSILON) { dissipative = false; break; }
    }
    dr->is_dissipative = dissipative;
    if (dissipative) dr->P[0][0] = 1.0;
    return dr;
}

void ic_dissipativity_result_free(DissipativityResult* dr) { free(dr); }

void ic_dissipativity_result_print(const DissipativityResult* dr) {
    printf("DissipativityResult: dissipative=%s  n=%d\n",
           dr->is_dissipative ? "YES" : "NO", dr->n);
}

bool ic_dissipativity_intercon_stable(const SupplyRate* sr1,
                                        const SupplyRate* sr2) {
    /* Interconnection of (Q1,S1,R1) and (Q2,S2,R2) with u1=-y2+v1, u2=y1+v2
     * is dissipative wrt combined supply if:
     * [Q1+R2   S1   ] > 0 (positive definite)
     * [S1^T    R1+Q2]
     */
    double M[2][2] = {
        {sr1->Q[0][0] + sr2->R[0][0], sr1->S[0][0]},
        {sr1->S[0][0], sr1->R[0][0] + sr2->Q[0][0]}
    };
    /* Positive definiteness: M[0][0] > 0 and det(M) > 0 */
    double det = M[0][0] * M[1][1] - M[0][1] * M[1][0];
    return M[0][0] > IC_EPSILON && det > IC_EPSILON;
}
double ic_passivity_l2_gain_bound(const ICSubsystem* s) {
    PassivityResult* pr = ic_passivity_test(s);
    if (!pr || !pr->is_passive) { ic_passivity_result_free(pr); return 1e6; }
    double gain = (pr->type == PASSIVE_STRICT) ? 1.0/fmax(pr->excess, 1e-6) : 1e6;
    ic_passivity_result_free(pr);
    return gain;
}
bool ic_is_finite_gain_stable(const ICSubsystem* s) {
    return ic_passivity_l2_gain_bound(s) < 1e5;
}
bool ic_mixed_passivity_small_gain(const ICSubsystem* s1, const ICSubsystem* s2) {
    double nu = ic_passivity_index(s1);
    double gamma = ic_l2_gain_hinf(s2);
    return (nu + 1.0/fmax(gamma, 1e-6)) > 0.0;
}
bool ic_incremental_passivity_test(const ICSubsystem* s) {
    PassivityResult* pr = ic_passivity_test(s);
    if (!pr) return false;
    bool inc_p = pr->is_passive && (pr->type >= PASSIVE_LOSSY);
    ic_passivity_result_free(pr);
    return inc_p;
}

/* Check if system is output strictly passive with margin rho */
bool ic_output_strict_passivity_margin(const ICSubsystem* s, double* rho) {
    PassivityResult* pr = ic_passivity_test(s);
    if (!pr) return false;
    bool osp = (pr->type == PASSIVE_OUTPUT_STRICT || pr->type == PASSIVE_STRICT);
    if (rho) *rho = osp ? pr->excess : 0.0;
    ic_passivity_result_free(pr);
    return osp;
}

/* Compute the largest sector for which the system is passive */
double ic_passivity_sector_max(const ICSubsystem* s) {
    double max_k = 0.0;
    for (double k = 0.1; k < 100.0; k *= 1.5) {
        double nu = ic_passivity_index(s);
        if (nu + k > 0.0) max_k = k;
        else break;
    }
    return max_k;
}

/* Check if negative feedback preserves passivity */
bool ic_feedback_preserves_passivity(const ICSubsystem* s1, const ICSubsystem* s2) {
    PassivityResult* pr1 = ic_passivity_test(s1);
    PassivityResult* pr2 = ic_passivity_test(s2);
    if (!pr1 || !pr2) { ic_passivity_result_free(pr1); ic_passivity_result_free(pr2); return false; }
    bool preserve = pr1->is_passive && pr2->is_passive;
    ic_passivity_result_free(pr1);
    ic_passivity_result_free(pr2);
    return preserve;
}
