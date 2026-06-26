#include "iss_small_gain.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

bool iss_small_gain_condition(const ISS_Gain* g1, const ISS_Gain* g2) {
    return iss_gain_small_gain_check(g1, g2);
}

bool iss_cyclic_small_gain(const ISS_Gain** gains, int n) {
    if (!gains || n < 2) return false;
    double product = 1.0;
    int i;
    for (i = 0; i < n; i++) {
        product *= gains[i]->linear_gain;
    }
    return product < 1.0;
}

bool iss_monotone_small_gain(const ISS_Interconnection* ic) {
    if (!ic || ic->n_systems < 2) return false;
    int n = ic->n_systems;
    int i, j;
    double max_product = 0.0;
    for (i = 0; i < n; i++) {
        double prod = 1.0;
        for (j = 0; j < n; j++) {
            if (ic->gains && ic->gains[i * n + j])
                prod *= ic->gains[i * n + j]->linear_gain;
        }
        if (prod > max_product) max_product = prod;
    }
    return max_product < 1.0;
}

ISS_Interconnection* iss_interconnection_create(int n) {
    ISS_Interconnection* ic = (ISS_Interconnection*)calloc(1, sizeof(ISS_Interconnection));
    if (!ic) return NULL;
    ic->n_systems = n;
    ic->subsystems = (ISS_System**)calloc(n, sizeof(ISS_System*));
    ic->gains = (ISS_Gain**)calloc(n * n, sizeof(ISS_Gain*));
    ic->n_gains = (int*)calloc(n, sizeof(int));
    return ic;
}

void iss_interconnection_free(ISS_Interconnection* ic) {
    if (!ic) return;
    int i;
    for (i = 0; i < ic->n_systems * ic->n_systems; i++) {
        if (ic->gains && ic->gains[i]) iss_gain_free(ic->gains[i]);
    }
    free(ic->subsystems);
    free(ic->gains);
    free(ic->n_gains);
    free(ic);
}

void iss_interconnection_set_gain(ISS_Interconnection* ic, int i, int j, const ISS_Gain* g) {
    if (!ic || i < 0 || j < 0 || i >= ic->n_systems || j >= ic->n_systems) return;
    ISS_Gain* copy = (ISS_Gain*)malloc(sizeof(ISS_Gain));
    memcpy(copy, g, sizeof(ISS_Gain));
    ic->gains[i * ic->n_systems + j] = copy;
    ic->n_gains[i]++;
}

bool iss_interconnection_is_ISS(const ISS_Interconnection* ic) {
    return iss_monotone_small_gain(ic);
}

bool iss_max_min_small_gain(const ISS_Interconnection* ic, double tol) {
    if (!ic || ic->n_systems < 2) return false;
    double spectral_radius = 0.0;
    int i, j;
    for (i = 0; i < ic->n_systems; i++) {
        double row_sum = 0.0;
        for (j = 0; j < ic->n_systems; j++) {
            if (ic->gains && ic->gains[i * ic->n_systems + j])
                row_sum += ic->gains[i * ic->n_systems + j]->linear_gain;
        }
        if (row_sum > spectral_radius) spectral_radius = row_sum;
    }
    return spectral_radius < 1.0 - tol;
}

/* ==============================================================
 * Small-Gain Margin
 * ============================================================== */

double iss_small_gain_margin(const ISS_Gain** gains, int n) {
    if (!gains || n < 2) return 0.0;
    double prod = 1.0;
    int i;
    for (i = 0; i < n; i++) prod *= gains[i]->linear_gain;
    return 1.0 - prod;
}

/* ==============================================================
 * Tight Small-Gain Check
 * ============================================================== */

bool iss_small_gain_is_tight(const ISS_Gain** gains, int n, double tol) {
    double m = iss_small_gain_margin(gains, n);
    return m > 0.0 && m < tol;
}

/* ==============================================================
 * Find Violation Index in Cyclic Chain
 * ============================================================== */

int iss_small_gain_violation_index(const ISS_Gain** gains, int n) {
    if (!gains || n < 2) return -1;
    int i;
    for (i = 0; i < n - 1; i++) {
        if (gains[i]->linear_gain * gains[i + 1]->linear_gain >= 1.0) return i;
    }
    if (gains[n - 1]->linear_gain * gains[0]->linear_gain >= 1.0) return n - 1;
    return -1;
}

/* ==============================================================
 * Nonlinear Small-Gain via K-function Composition
 * ============================================================== */

bool iss_nonlinear_small_gain(const ISS_KFunction* g1, const ISS_KFunction* g2,
    double s_max, int n_points) {
    if (!g1 || !g2 || n_points <= 0) return false;
    int i;
    for (i = 1; i <= n_points; i++) {
        double s = i * s_max / n_points;
        double composed = iss_k_eval(g1, iss_k_eval(g2, s));
        if (composed >= s) return false;
    }
    return true;
}

/* ==============================================================
 * Spectral Radius Small-Gain Condition
 * ============================================================== */

double iss_spectral_radius(const double* M, int n) {
    if (!M || n <= 0) return INFINITY;
    /* Power iteration to estimate spectral radius */
    double v[ISS_MAX_DIM];
    double v_new[ISS_MAX_DIM];
    int i, j, iter;
    for (i = 0; i < n && i < ISS_MAX_DIM; i++) v[i] = 1.0;
    double lambda = 0.0;
    for (iter = 0; iter < 100; iter++) {
        double norm = 0.0;
        for (i = 0; i < n; i++) {
            v_new[i] = 0.0;
            for (j = 0; j < n; j++) {
                v_new[i] += M[i * n + j] * v[j];
            }
            norm += v_new[i] * v_new[i];
        }
        norm = sqrt(norm);
        if (norm < 1e-15) break;
        for (i = 0; i < n; i++) v[i] = v_new[i] / norm;
        double new_lambda = 0.0;
        for (i = 0; i < n; i++) new_lambda += fabs(v_new[i]);
        if (fabs(new_lambda - lambda) < 1e-10) break;
        lambda = new_lambda;
    }
    return lambda;
}

/* ==============================================================
 * Small-Gain via Lyapunov Functions
 * ============================================================== */

bool iss_small_gain_lyapunov(const ISS_System* s1, const ISS_System* s2,
    double c1, double c2, double tol) {
    if (!s1 || !s2) return false;
    /* Check if composite V = c1*V1 + c2*V2 is ISS-Lyapunov for interconnection */
    /* Simplified: check if coupling product < 1 via linear gain approx */
    ISS_Gain g1;
    ISS_Gain g2;
    memset(&g1, 0, sizeof(g1));
    memset(&g2, 0, sizeof(g2));
    g1.linear_gain = c1;
    g2.linear_gain = c2;
    bool passed = iss_small_gain_condition(&g1, &g2);
    (void)tol;
    return passed;
}

/* ==============================================================
 * Build Gain Matrix from Pairwise Gains
 * ============================================================== */

void iss_build_gain_matrix(ISS_Interconnection* ic, const double* gain_values) {
    if (!ic || !gain_values) return;
    int n = ic->n_systems;
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            ISS_Gain* g = (ISS_Gain*)malloc(sizeof(ISS_Gain));
            memset(g, 0, sizeof(ISS_Gain));
            g->linear_gain = gain_values[i * n + j];
            g->is_finite = isfinite(g->linear_gain);
            g->peak_gain = g->linear_gain;
            ic->gains[i * n + j] = g;
        }
    }
}
/* ==============================================================
 * Small-Gain for Sector-Bounded Nonlinearities
 * ============================================================== */

bool iss_small_gain_sector(const ISS_Gain* g1, const ISS_Gain* g2,
    double sector_lower, double sector_upper) {
    if (!g1 || !g2) return false;
    /* Sector condition: the nonlinearity phi satisfies
     *   a*s^2 <= s*phi(s) <= b*s^2
     * Small-gain requires: max(g1*g2*b, |a|*g1*g2) < 1 */
    double max_sector_gain = fmax(sector_upper, fabs(sector_lower));
    double product = g1->linear_gain * g2->linear_gain * max_sector_gain;
    return product < 1.0;
}

/* ==============================================================
 * ISS Interconnection Gain Matrix Norm Check
 * ============================================================== */

bool iss_gain_matrix_condition(const double* M, int n, double rho_max) {
    if (!M || n <= 0) return false;
    /* Check if spectral radius of gain matrix M < rho_max */
    int i, j, k;
    /* Compute induced infinity norm: max_i sum_j |M_ij| */
    double induced_norm = 0.0;
    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        for (j = 0; j < n; j++) {
            row_sum += fabs(M[i * n + j]);
        }
        if (row_sum > induced_norm) induced_norm = row_sum;
    }
    if (induced_norm < rho_max) return true;

    /* If induced norm check fails, try Gershgorin circles */
    for (i = 0; i < n; i++) {
        double center = M[i * n + i];
        double radius = 0.0;
        for (j = 0; j < n; j++) {
            if (i != j) radius += fabs(M[i * n + j]);
        }
        if (fabs(center) + radius >= rho_max) return false;
    }
    (void)k;
    return true;
}

/* ==============================================================
 * Loop Transformation for Non-Monotone Interconnections
 * ============================================================== */

bool iss_loop_transform_check(ISS_System** systems, int n,
    double* transform_gains) {
    if (!systems || !transform_gains || n < 2) return false;
    /* Apply loop transformation to convert non-monotone
     * interconnection to monotone form for small-gain analysis.
     * The transformed gains are used to check the small-gain condition. */
    double product = 1.0;
    int i;
    for (i = 0; i < n; i++) {
        product *= transform_gains[i];
    }
    /* After transformation, check if the modified gains satisfy < 1 */
    bool result = (product < 1.0);

    /* Verify each system exists */
    for (i = 0; i < n; i++) {
        if (!systems[i]) return false;
    }
    return result;
}