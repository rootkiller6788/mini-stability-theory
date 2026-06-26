#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fts_core.h"
#include "fts_homogeneous.h"

/* ============================================================
 * fts_homogeneous.c -- Homogeneity-Based Finite-Time Stability
 *
 * A vector field f(x) is homogeneous of degree d with weights
 * (r1,...,rn) if for all lambda > 0:
 *   f_i(lambda^{r1} x1, ..., lambda^{rn} xn) = lambda^{d+ri} f_i(x)
 *
 * FTS Theorem: If origin is asymptotically stable and f has
 * negative homogeneity degree (d < 0), then origin is FTS.
 *
 * Ref: Bhat & Bernstein (2005), Mathematics of Control,
 * Signals, and Systems, 17:101-127.
 * ============================================================ */

/* --- Lifecycle --- */

FTSHomogeneousWeights* fts_hom_weights_create(int dim, double* w) {
    FTSHomogeneousWeights* hw = calloc(1, sizeof(FTSHomogeneousWeights));
    if (!hw) return NULL;
    hw->n = dim;
    hw->weights = calloc(dim, sizeof(double));
    if (!hw->weights) { free(hw); return NULL; }
    if (w) {
        memcpy(hw->weights, w, dim * sizeof(double));
    } else {
        for (int i = 0; i < dim; i++) hw->weights[i] = 1.0;
    }
    hw->degree = 0.0;
    hw->is_negative = false;
    return hw;
}

void fts_hom_weights_free(FTSHomogeneousWeights* hw) {
    if (hw) { free(hw->weights); free(hw); }
}

/* --- Dilation Operator --- */

FTSState fts_hom_dilate(FTSState x, FTSHomogeneousWeights* hw,
                         double lambda) {
    /* Delta_lambda(x) = (lambda^{r1}*x1, ..., lambda^{rn}*xn)
     *
     * The dilation forms a one-parameter transformation group.
     * A function f is homogeneous of degree d w.r.t. weights r_i if:
     *   f(Delta_lambda(x)) = lambda^d * Delta_lambda(f(x))
     * which component-wise is:
     *   f_i(Delta_lambda(x)) = lambda^{d+ri} * f_i(x)
     *
     * Key property: Delta_{a} o Delta_{b} = Delta_{a*b}
     */
    FTSState result = x;
    if (!hw) return result;
    for (int i = 0; i < x.dim && i < hw->n; i++) {
        double scaled = pow(fabs(lambda), hw->weights[i]);
        result.x[i] = scaled * x.x[i];
    }
    return result;
}

void fts_hom_dilate_inplace(FTSState* x, FTSHomogeneousWeights* hw,
                             double lambda) {
    if (!x || !hw) return;
    for (int i = 0; i < x->dim && i < hw->n; i++)
        x->x[i] = pow(fabs(lambda), hw->weights[i]) * x->x[i];
}

/* --- Homogeneity Verification --- */

bool fts_hom_is_homogeneous(FTSSystem* sys, FTSHomogeneousWeights* hw,
                             double deg) {
    if (!sys || !hw) return false;
    FTSState x0; x0.dim = sys->dim;
    for (int i = 0; i < sys->dim; i++) x0.x[i] = 1.0;
    FTSState f_x0; f_x0.dim = sys->dim;
    sys->rhs(sys, x0, &f_x0);
    double lambdas[] = {0.5, 0.7, 1.0, 1.5, 2.0};
    double max_residual = 0.0;
    for (int k = 0; k < 5; k++) {
        double lam = lambdas[k];
        FTSState x_dil = fts_hom_dilate(x0, hw, lam);
        FTSState f_dil; f_dil.dim = sys->dim;
        sys->rhs(sys, x_dil, &f_dil);
        for (int i = 0; i < sys->dim; i++) {
            double expected = pow(lam, deg + hw->weights[i]) * f_x0.x[i];
            double residual = fabs(f_dil.x[i] - expected);
            double scale = fabs(expected) + 1e-6;
            if (residual / scale > max_residual)
                max_residual = residual / scale;
        }
    }
    return max_residual < 1e-3;
}

bool fts_hom_verify_field(FTSSystem* sys, FTSHomogeneousWeights* hw,
                           FTSState xt, double deg, double* residual) {
    if (!sys || !hw || !residual) return false;
    *residual = 0.0;
    FTSState fx; fx.dim = sys->dim;
    sys->rhs(sys, xt, &fx);
    double lambdas[] = {0.5, 1.0, 2.0};
    for (int k = 0; k < 3; k++) {
        double lam = lambdas[k];
        FTSState xd = fts_hom_dilate(xt, hw, lam);
        FTSState fd; fd.dim = sys->dim;
        sys->rhs(sys, xd, &fd);
        for (int i = 0; i < sys->dim; i++) {
            double exp = pow(lam, deg + hw->weights[i]) * fx.x[i];
            double err = fabs(fd.x[i] - exp);
            double rel = err / (fabs(exp) + 1e-9);
            if (rel > *residual) *residual = rel;
        }
    }
    return (*residual < 1e-3);
}

double fts_hom_degree_estimate(FTSSystem* sys, FTSState x0, int n_tests) {
    if (!sys || n_tests < 1) return 0.0;
    double best_deg = 0.0, best_err = 1e99;
    double weights[10];
    for (int i = 0; i < sys->dim && i < 10; i++) weights[i] = 1.0;
    for (int d_int = -20; d_int <= 20; d_int++) {
        double d = d_int * 0.1;
        FTSHomogeneousWeights* hw = fts_hom_weights_create(sys->dim, weights);
        if (!hw) continue;
        double total_err = 0.0;
        FTSState xt = x0;
        int actual = n_tests < 5 ? n_tests : 5;
        for (int t = 0; t < actual; t++) {
            double residual;
            fts_hom_verify_field(sys, hw, xt, d, &residual);
            total_err += residual;
            for (int i = 0; i < sys->dim; i++)
                xt.x[i] = x0.x[i] * (1.0 + 0.1*(t+1));
        }
        if (total_err < best_err) { best_err = total_err; best_deg = d; }
        fts_hom_weights_free(hw);
    }
    return best_deg;
}

/* --- FTS Detection --- */

bool fts_hom_is_finite_time(FTSHomogeneousWeights* hw, double deg) {
    return hw && deg < 0.0;
}

bool fts_hom_has_negative_degree(FTSHomogeneousWeights* hw) {
    return hw && hw->is_negative;
}

/* --- Homogeneous Norms --- */

double fts_hom_norm(FTSState x, FTSHomogeneousWeights* hw) {
    if (!hw) return fts_norm(x);
    double sum = 0.0;
    for (int i = 0; i < x.dim && i < hw->n; i++) {
        if (hw->weights[i] < 1e-12) continue;
        double exponent = 2.0 / hw->weights[i];
        sum += pow(fabs(x.x[i]), exponent);
    }
    return pow(sum, 0.5);
}

double fts_hom_weighted_norm(FTSState x, FTSHomogeneousWeights* hw,
                              double p) {
    if (!hw || p < 1e-12) return fts_norm(x);
    double sum = 0.0;
    for (int i = 0; i < x.dim && i < hw->n; i++) {
        if (hw->weights[i] < 1e-12) continue;
        sum += pow(fabs(x.x[i]), p / hw->weights[i]);
    }
    return pow(sum, 1.0 / p);
}

double fts_hom_euclidean_norm(FTSState x, FTSHomogeneousWeights* hw) {
    if (!hw) return fts_norm(x);
    double sum = 0.0;
    for (int i = 0; i < x.dim && i < hw->n; i++)
        sum += hw->weights[i] * x.x[i] * x.x[i];
    return sqrt(sum);
}

/* --- Convergence Testing --- */

FTSHomogeneousTest* fts_hom_test_convergence(FTSSystem* sys,
                                              FTSState x0,
                                              FTSHomogeneousWeights* hw,
                                              double T, double dt) {
    if (!sys) return NULL;
    FTSHomogeneousTest* ht = calloc(1, sizeof(FTSHomogeneousTest));
    if (!ht) return NULL;
    ht->x0 = x0;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) { free(ht); return NULL; }
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    ht->converged = false;
    ht->settling_time = T;
    int max_steps = (int)(T / dt);
    if (max_steps <= 0) max_steps = 1000;
    const double tol = 1e-3;
    for (int i = 0; i < max_steps; i++) {
        fts_step_rk4(cp);
        ht->steps++;
        double nrm = hw ? fts_hom_norm(cp->state, hw)
                        : fts_norm(cp->state);
        if (nrm < tol) {
            ht->converged = true;
            ht->settling_time = cp->t;
            ht->xf = cp->state;
            ht->accuracy = nrm;
            break;
        }
    }
    if (!ht->converged) {
        ht->xf = cp->state;
        ht->accuracy = hw ? fts_hom_norm(cp->state, hw)
                          : fts_norm(cp->state);
    }
    fts_free(cp);
    return ht;
}

void fts_hom_test_free(FTSHomogeneousTest* ht) { free(ht); }

/* --- Bounds --- */

double fts_hom_settling_time_bound(FTSHomogeneousWeights* hw,
                                    double deg, FTSState x0) {
    if (!hw || deg >= 0.0) return 1e99;
    double nrm = fts_hom_norm(x0, hw);
    if (nrm < 1e-12) return 0.0;
    return pow(nrm, -deg);
}

double fts_hom_settling_time_upper(FTSHomogeneousWeights* hw,
                                    double deg, FTSState x0, double k) {
    if (!hw || deg >= 0.0) return 1e99;
    double nrm = fts_hom_norm(x0, hw);
    if (nrm < 1e-12) return 0.0;
    return k * pow(nrm, -deg);
}

/* --- Output --- */

void fts_hom_print(FTSHomogeneousWeights* hw) {
    if (!hw) return;
    printf("Homogeneous: n=%d deg=%.4f neg=%s\n",
           hw->n, hw->degree, hw->is_negative ? "yes" : "no");
}

void fts_hom_test_print(FTSHomogeneousTest* ht) {
    if (!ht) return;
    printf("Hom Test: conv=%s Ts=%.4f acc=%.2e steps=%d\n",
           ht->converged ? "yes" : "no",
           ht->settling_time, ht->accuracy, ht->steps);
}

int fts_hom_compare_weights(FTSHomogeneousWeights* a,
                             FTSHomogeneousWeights* b) {
    if (!a && !b) return 0;
    if (!a) return -1; if (!b) return 1;
    int n = (a->n < b->n) ? a->n : b->n;
    for (int i = 0; i < n; i++) {
        if (a->weights[i] < b->weights[i] - 1e-12) return -1;
        if (a->weights[i] > b->weights[i] + 1e-12) return 1;
    }
    if (a->n < b->n) return -1;
    if (a->n > b->n) return 1;
    return 0;
}

/* --- Presets --- */

FTSHomogeneousWeights* fts_hom_standard_weights(int dim) {
    double* w = calloc(dim, sizeof(double));
    if (!w) return NULL;
    for (int i = 0; i < dim; i++)
        w[i] = (double)(dim - i);
    FTSHomogeneousWeights* hw = fts_hom_weights_create(dim, w);
    free(w);
    return hw;
}

FTSHomogeneousWeights* fts_hom_uniform_weights(int dim, double wv) {
    double* w = calloc(dim, sizeof(double));
    if (!w) return NULL;
    for (int i = 0; i < dim; i++) w[i] = wv;
    FTSHomogeneousWeights* hw = fts_hom_weights_create(dim, w);
    free(w);
    return hw;
}