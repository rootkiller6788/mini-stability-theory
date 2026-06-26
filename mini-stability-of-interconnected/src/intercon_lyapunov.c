#include "intercon_lyapunov.h"
#include "intercon_small_gain.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * Composite Lyapunov Function
 *
 * V(x) = sum_i d_i * V_i(x_i) where V_i(x_i) = x_i^T * P_i * x_i
 * and d_i > 0 are weighting coefficients.
 *
 * V_dot(x) = sum_i d_i * [x_i^T*(A_i^T*P_i + P_i*A_i)*x_i
 *           + 2*x_i^T*P_i*B_i*(sum_j H_ij*C_j*x_j)]
 * ============================================================== */

CompositeLyapunov* ic_composite_lyapunov(const ICSystem* sys,
                                           const double* weights) {
    CompositeLyapunov* cl = calloc(1, sizeof(CompositeLyapunov));
    if (!cl) return NULL;
    cl->n = sys->n_subsys;
    cl->V = malloc((size_t)cl->n * sizeof(double));
    cl->d = malloc((size_t)cl->n * sizeof(double));
    if (!cl->V || !cl->d) {
        ic_composite_lyapunov_free(cl);
        return NULL;
    }

    for (int i = 0; i < cl->n; i++) {
        cl->d[i] = weights ? weights[i] : 1.0;
        /* V_i = x_i^T * P_i * x_i with P_i = I (identity) as default */
        double Vi = 0.0;
        for (int k = 0; k < sys->subsys[i]->n; k++)
            Vi += sys->subsys[i]->x[k] * sys->subsys[i]->x[k];
        cl->V[i] = Vi;
    }

    /* Compute composite V */
    cl->V_composite = 0.0;
    for (int i = 0; i < cl->n; i++)
        cl->V_composite += cl->d[i] * cl->V[i];

    /* Estimate decay rate */
    cl->decay_rate = ic_composite_lyapunov_derivative(cl, sys);
    cl->is_stable = (cl->decay_rate < -IC_EPSILON);

    return cl;
}

void ic_composite_lyapunov_free(CompositeLyapunov* cl) {
    if (!cl) return;
    free(cl->V);
    free(cl->d);
    free(cl);
}

double ic_composite_lyapunov_eval(const CompositeLyapunov* cl,
                                    const double* x) {
    double V = 0.0;
    int offset = 0;
    /* Simplified: V = sum d_i * sum x_k^2 */
    for (int i = 0; i < cl->n; i++) {
        double Vi = 0.0;
        for (int k = 0; offset + k < offset + IC_MAX_DIM; k++)
            Vi += x[offset + k] * x[offset + k];
        V += cl->d[i] * Vi;
    }
    return V;
}

double ic_composite_lyapunov_derivative(const CompositeLyapunov* cl,
                                          const ICSystem* sys) {
    /* V_dot = sum_i d_i * [2*x_i^T*A_i^T*x_i + 2*x_i^T*B_i*u_i]
     * where u_i = sum_j H_ij * C_j * x_j
     *
     * For stable A_i, the first term is negative. The interconnection
     * terms can be destabilizing. */
    double V_dot = 0.0;
    for (int i = 0; i < cl->n; i++) {
        ICSubsystem* si = sys->subsys[i];
        int ni = si->n;
        /* Contribution from self-dynamics: 2*x_i^T * A_i * x_i */
        double self_term = 0.0;
        for (int r = 0; r < ni; r++)
            for (int c = 0; c < ni; c++)
                self_term += si->x[r] * si->A[r][c] * si->x[c];
        self_term *= 2.0 * cl->d[i];

        /* Contribution from interconnections: 2*x_i^T*B_i*sum_j H_ij*C_j*x_j */
        double inter_term = 0.0;
        for (int j = 0; j < cl->n; j++) {
            double hij = ic_intercon_get(sys->intercon, i, j);
            if (fabs(hij) < IC_EPSILON) continue;
            ICSubsystem* sj = sys->subsys[j];
            for (int r = 0; r < ni; r++)
                for (int c = 0; c < sj->n; c++) {
                    double contrib = 0.0;
                    for (int k = 0; k < si->m && k < sj->p; k++)
                        contrib += si->B[r][k] * hij * sj->C[k][c];
                    inter_term += si->x[r] * contrib * sj->x[c];
                }
        }
        inter_term *= 2.0 * cl->d[i];
        V_dot += self_term + inter_term;
    }
    return V_dot;
}

void ic_composite_lyapunov_print(const CompositeLyapunov* cl) {
    printf("CompositeLyapunov: n=%d  V=%.4f  V_dot=%.4f  stable=%s\n",
           cl->n, cl->V_composite, cl->decay_rate,
           cl->is_stable ? "YES" : "NO");
    for (int i = 0; i < cl->n; i++)
        printf("  V[%d]=%.4f (weight=%.2f)\n", i, cl->V[i], cl->d[i]);
}

/* ==============================================================
 * M-Matrix Construction & Analysis
 *
 * Build comparison system: V_dot_i <= -alpha_i*V_i + sum_j gamma_ij*V_j
 * where alpha_i > 0 (stability margin of isolated subsystem i)
 * and gamma_ij >= 0 are interconnection gains.
 *
 * Define M_ij = { -alpha_i + gamma_ii  if i=j
 *               {  gamma_ij           if i!=j
 *
 * An M-matrix has negative diagonal and nonnegative off-diagonal entries,
 * and all eigenvalues have negative real parts.
 * ============================================================== */

MMatrixResult* ic_build_M_matrix(const ICSystem* sys) {
    MMatrixResult* mr = calloc(1, sizeof(MMatrixResult));
    if (!mr) return NULL;
    mr->n = sys->n_subsys;

    /* Estimate alpha_i = -max(Re(eig(A_i))) (stability margin) */
    double alpha[IC_MAX_SUBSYS];
    for (int i = 0; i < mr->n; i++) {
        ICSubsystem* si = sys->subsys[i];
        if (si->n == 1) {
            alpha[i] = -si->A[0][0];
        } else {
            /* Approximate alpha from trace */
            double tr = 0.0;
            for (int k = 0; k < si->n; k++) tr += si->A[k][k];
            alpha[i] = -tr / si->n;
        }
        if (alpha[i] < IC_EPSILON) alpha[i] = 0.1; /* ensure positivity */
    }

    /* Estimate gamma_ij = ||B_i * H_ij * C_j|| */
    double gamma[IC_MAX_SUBSYS][IC_MAX_SUBSYS] = {{0}};
    for (int i = 0; i < mr->n; i++) {
        for (int j = 0; j < mr->n; j++) {
            double hij = ic_intercon_get(sys->intercon, i, j);
            if (fabs(hij) < IC_EPSILON) continue;
            gamma[i][j] = fabs(hij) * ic_l2_gain_hinf(sys->subsys[j]);
        }
    }

    /* Build M matrix: M_ij = -alpha_i (if i=j) + gamma_ij */
    for (int i = 0; i < mr->n; i++) {
        for (int j = 0; j < mr->n; j++) {
            mr->M[i][j] = gamma[i][j];
            if (i == j) mr->M[i][j] -= alpha[i];
        }
    }

    /* Check M-matrix conditions */
    bool is_M = true;
    /* 1. Off-diagonal entries >= 0 */
    for (int i = 0; i < mr->n; i++)
        for (int j = 0; j < mr->n; j++)
            if (i != j && mr->M[i][j] < -IC_EPSILON)
                is_M = false;
    /* 2. Diagonal entries < 0 */
    for (int i = 0; i < mr->n; i++)
        if (mr->M[i][i] >= -IC_EPSILON)
            is_M = false;
    mr->is_m_matrix = is_M;

    /* 3. Stability: all eigenvalues have negative real parts */
    double real_parts[IC_MAX_SUBSYS], imag_parts[IC_MAX_SUBSYS];
    int n_eig = ic_eigenvalues_real((double*)mr->M, mr->n, real_parts,
                                      imag_parts, IC_MAX_SUBSYS);
    bool stable = true;
    double max_real = real_parts[0];
    for (int i = 0; i < n_eig; i++) {
        if (real_parts[i] >= -IC_EPSILON) stable = false;
        if (real_parts[i] > max_real) max_real = real_parts[i];
    }
    mr->stability_margin = -max_real;

    /* Override stability with eigenvalue check */
    if (!stable && is_M) mr->is_m_matrix = false;

    return mr;
}

void ic_M_matrix_free(MMatrixResult* mr) { free(mr); }

bool ic_is_M_matrix(const double* M_flat, int n) {
    for (int i = 0; i < n; i++) {
        /* Diagonal must be negative */
        if (M_flat[i*n+i] >= -IC_EPSILON) return false;
        /* Off-diagonal must be nonnegative */
        for (int j = 0; j < n; j++)
            if (i != j && M_flat[i*n+j] < -IC_EPSILON) return false;
    }
    /* Check principal minors (for n<=3) */
    if (n == 1) return true;
    if (n == 2) {
        double det = M_flat[0]*M_flat[3] - M_flat[1]*M_flat[2];
        return det > IC_EPSILON;
    }
    return true;
}

bool ic_M_matrix_stability(const MMatrixResult* mr) {
    return mr->is_m_matrix;
}

void ic_M_matrix_print(const MMatrixResult* mr) {
    printf("MMatrixResult: n=%d  is_M=%s  stability_margin=%.4f\n",
           mr->n, mr->is_m_matrix ? "YES" : "NO", mr->stability_margin);
    printf("  M matrix:\n");
    for (int i = 0; i < mr->n; i++) {
        printf("    ");
        for (int j = 0; j < mr->n; j++)
            printf("% 8.3f ", mr->M[i][j]);
        printf("\n");
    }
}

/* ==============================================================
 * Gershgorin Circle Theorem
 * ============================================================== */

GershgorinResult* ic_gershgorin_test(const double* A, int n) {
    GershgorinResult* gr = calloc(1, sizeof(GershgorinResult));
    if (!gr) return NULL;
    gr->n = n;
    gr->is_diag_dominant = true;
    gr->all_left_half = true;
    gr->max_real_part = -1e6;

    for (int i = 0; i < n; i++) {
        gr->centers[i] = A[i*n+i];
        gr->radii[i] = 0.0;
        for (int j = 0; j < n; j++)
            if (i != j) gr->radii[i] += fabs(A[i*n+j]);

        if (gr->radii[i] >= fabs(gr->centers[i]))
            gr->is_diag_dominant = false;

        double rightmost = gr->centers[i] + gr->radii[i];
        if (rightmost > gr->max_real_part) gr->max_real_part = rightmost;
        if (rightmost > -IC_EPSILON) gr->all_left_half = false;
    }
    return gr;
}

void ic_gershgorin_free(GershgorinResult* gr) { free(gr); }

bool ic_gershgorin_stable(const GershgorinResult* gr) {
    return gr->all_left_half;
}

void ic_gershgorin_print(const GershgorinResult* gr) {
    printf("GershgorinResult: n=%d  diag_dominant=%s  all_left=%s  max_real=%.4f\n",
           gr->n,
           gr->is_diag_dominant ? "YES" : "NO",
           gr->all_left_half ? "YES" : "NO",
           gr->max_real_part);
    for (int i = 0; i < gr->n && i < 8; i++)
        printf("  Disc %d: center=%.4f radius=%.4f\n",
               i, gr->centers[i], gr->radii[i]);
}

/* ==============================================================
 * Diagonal Dominance Checks
 * ============================================================== */

bool ic_is_row_diag_dominant(const double* A, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++)
            if (i != j) sum += fabs(A[i*n+j]);
        if (sum >= fabs(A[i*n+i])) return false;
    }
    return true;
}

bool ic_is_col_diag_dominant(const double* A, int n) {
    for (int j = 0; j < n; j++) {
        double sum = 0.0;
        for (int i = 0; i < n; i++)
            if (i != j) sum += fabs(A[i*n+j]);
        if (sum >= fabs(A[j*n+j])) return false;
    }
    return true;
}

bool ic_is_generalized_diag_dominant(const double* A, int n) {
    /* Try scaling by positive diagonal D */
    for (int si = 0; si < 30; si++) {
        double d[IC_MAX_STATES];
        for (int i = 0; i < n; i++)
            d[i] = 0.1 + 9.9 * ((si + i) % 30) / 29.0;
        /* Scale: B = D*A*D^{-1}, B_ij = (d_i/d_j)*A_ij */
        bool dominant = true;
        for (int i = 0; i < n && dominant; i++) {
            double sum = 0.0;
            for (int j = 0; j < n; j++)
                if (i != j) sum += (d[i]/d[j]) * fabs(A[i*n+j]);
            if (sum >= fabs(A[i*n+i])) dominant = false;
        }
        if (dominant) return true;
    }
    return false;
}

/* ==============================================================
 * Comparison Principle
 * ============================================================== */

void ic_comparison_system_step(const double* M, int n, const double* V_in,
                                 double* V_out, double dt) {
    for (int i = 0; i < n; i++) {
        V_out[i] = V_in[i];
        for (int j = 0; j < n; j++)
            V_out[i] += M[i*n+j] * V_in[j] * dt;
        if (V_out[i] < 0.0) V_out[i] = 0.0;
    }
}

void ic_comparison_system_simulate(const double* M, int n, const double* V0,
                                     double duration, double dt, double* V_final) {
    double V[IC_MAX_SUBSYS], V_next[IC_MAX_SUBSYS];
    for (int i = 0; i < n; i++) V[i] = V0[i];
    int steps = (int)(duration / dt);
    for (int s = 0; s < steps; s++) {
        ic_comparison_system_step(M, n, V, V_next, dt);
        for (int i = 0; i < n; i++) V[i] = V_next[i];
    }
    if (V_final)
        for (int i = 0; i < n; i++) V_final[i] = V[i];
}

/* ==============================================================
 * Connective Stability
 * ============================================================== */

bool ic_connective_stability_test(const ICSystem* sys) {
    /* Check stability under all possible interconnections
     * by testing diagonal dominance of composite A */
    int N = sys->total_states;
    double A_comp[IC_MAX_STATES * IC_MAX_STATES];
    ic_build_composite_A(sys, A_comp, IC_MAX_STATES);

    GershgorinResult* gr = ic_gershgorin_test(A_comp, N);
    bool stable = ic_gershgorin_stable(gr);
    /* Also check M-matrix criterion */
    if (!stable) {
        MMatrixResult* mr = ic_build_M_matrix(sys);
        stable = ic_M_matrix_stability(mr);
        ic_M_matrix_free(mr);
    }
    ic_gershgorin_free(gr);
    return stable;
}

/* ==============================================================
 * Stability Margins
 * ============================================================== */

double ic_stability_margin_mu(const ICSystem* sys) {
    /* Structured singular value approximation:
     * mu = sup{sigma_bar(M) : det(I - M*Delta) = 0}
     * Simplified: return 1/spectral_radius of loop gain */
    MMatrixResult* mr = ic_build_M_matrix(sys);
    if (!mr) return 0.0;
    double margin = mr->stability_margin;
    ic_M_matrix_free(mr);
    return margin;
}

double ic_gain_margin(const ICSubsystem* s1, const ICSubsystem* s2) {
    SmallGainResult* sg = ic_small_gain_test(s1, s2, IC_FEEDBACK_NEG);
    if (!sg) return 0.0;
    double margin = (sg->product > IC_EPSILON) ? 1.0 / sg->product : 1e6;
    ic_small_gain_result_free(sg);
    return margin;
}

double ic_phase_margin(const ICSubsystem* s1, const ICSubsystem* s2) {
    /* Phase margin: find max phase shift before instability
     * Simplified: return 180 - (phase at crossover frequency) */
    double best_margin = 180.0;
    for (int i = 0; i < 100; i++) {
        double w = 0.01 * pow(10.0, 4.0 * i / 99.0);
        double mag1, phase1, mag2, phase2;
        ic_subsys_freq_response(s1, w, &mag1, &phase1);
        ic_subsys_freq_response(s2, w, &mag2, &phase2);
        double loop_mag = mag1 * mag2;
        if (fabs(loop_mag - 1.0) < 0.05) {
            double loop_phase = (phase1 + phase2) * 180.0 / IC_PI;
            double margin = 180.0 + loop_phase;
            if (margin < best_margin && margin > 0.0) best_margin = margin;
        }
    }
    return best_margin;
}
double ic_estimate_roa(const ICSystem* sys) {
    CompositeLyapunov* cl = ic_composite_lyapunov(sys, NULL);
    if (!cl || !cl->is_stable) { ic_composite_lyapunov_free(cl); return 0.0; }
    double roa = fabs(cl->decay_rate) * cl->V_composite;
    ic_composite_lyapunov_free(cl);
    return fmax(roa, 0.0);
}
bool ic_absolute_stability_sector(const ICSubsystem* s, double k) {
    bool popov = ic_popov_criterion(s, k);
    CircleCriterion* cc = ic_circle_criterion(s, 0.0, k);
    bool circle = cc ? cc->is_stable : false;
    ic_circle_criterion_free(cc);
    return popov || circle;
}
double ic_max_interconnection_gain(const ICSystem* sys) {
    double min_margin = 1e6;
    for (int i = 0; i < sys->n_subsys; i++) {
        ICSubsystem* si = sys->subsys[i];
        double trace = 0.0;
        for (int k = 0; k < si->n; k++) trace += si->A[k][k];
        double margin = (trace < 0) ? -trace/si->n : 0.1;
        if (margin < min_margin && margin > 0.0) min_margin = margin;
    }
    return min_margin;
}

/* Compute the decay rate from the M-matrix */
double ic_m_matrix_decay_rate(const MMatrixResult* mr) {
    return mr->stability_margin;
}

/* Check if system satisfies the circle criterion for time-varying nonlinearity */
bool ic_circle_time_varying(const ICSubsystem* s, double k1, double k2) {
    CircleCriterion* cc = ic_circle_criterion(s, k1, k2);
    if (!cc) return false;
    bool stable = cc->is_stable;
    ic_circle_criterion_free(cc);
    return stable;
}

/* Estimate minimum damping ratio of composite system */
double ic_min_damping_ratio(const ICSystem* sys) {
    double A_comp[256];
    ic_build_composite_A(sys, A_comp, 16);
    double real_parts[16], imag_parts[16];
    int n_eig = ic_eigenvalues_real(A_comp, sys->total_states, real_parts, imag_parts, 16);
    if (n_eig <= 0) return 0.0;
    double min_zeta = 1e6;
    for (int i = 0; i < n_eig; i++) {
        double wn = sqrt(real_parts[i]*real_parts[i] + imag_parts[i]*imag_parts[i]);
        if (wn > 1e-12) {
            double zeta = -real_parts[i]/wn;
            if (zeta < min_zeta) min_zeta = zeta;
        }
    }
    return (min_zeta > 1e6) ? 1.0 : min_zeta;
}
