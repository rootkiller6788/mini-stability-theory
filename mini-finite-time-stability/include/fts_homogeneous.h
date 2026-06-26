#ifndef FTS_HOMOGENEOUS_H
#define FTS_HOMOGENEOUS_H

#include "fts_core.h"

/* ============================================================
 * Homogeneity-Based Finite-Time Stability
 *
 * A vector field f(x) is homogeneous of degree d with weights
 * (r1,...,rn) if:
 *   f_i(lambda^r1*x1, ..., lambda^rn*xn) = lambda^(d+ri) * f_i(x)
 *
 * FTS condition: origin asymptotically stable + f has negative
 * homogeneity degree (d < 0).
 *
 * Ref: Bhat & Bernstein (2005), "Geometric homogeneity with
 * applications to finite-time stability"
 * ============================================================ */

typedef struct {
    double* weights;    /* Homogeneity weights r_i > 0 */
    int n;              /* Number of weights (dimension) */
    double degree;      /* Homogeneity degree d */
    bool is_negative;   /* True if degree < 0 */
} FTSHomogeneousWeights;

typedef struct {
    FTSState x0;            /* Initial condition */
    FTSState xf;            /* Final state */
    double settling_time;   /* Settling time */
    bool converged;         /* Convergence detected */
    double accuracy;        /* Final accuracy */
    int steps;              /* Integration steps taken */
} FTSHomogeneousTest;

/* --- Lifecycle --- */
FTSHomogeneousWeights* fts_hom_weights_create(int dim, double* w);
void fts_hom_weights_free(FTSHomogeneousWeights* hw);

/* --- Dilation --- */
FTSState fts_hom_dilate(FTSState x, FTSHomogeneousWeights* hw, double lambda);
void fts_hom_dilate_inplace(FTSState* x, FTSHomogeneousWeights* hw,
                             double lambda);

/* --- Homogeneity Verification --- */
bool fts_hom_is_homogeneous(FTSSystem* sys, FTSHomogeneousWeights* hw,
                             double deg);
double fts_hom_degree_estimate(FTSSystem* sys, FTSState x0, int n_tests);
bool fts_hom_verify_field(FTSSystem* sys, FTSHomogeneousWeights* hw,
                           FTSState x_test, double deg, double* residual);

/* --- FTS Detection --- */
bool fts_hom_is_finite_time(FTSHomogeneousWeights* hw, double deg);
bool fts_hom_has_negative_degree(FTSHomogeneousWeights* hw);

/* --- Norms --- */
double fts_hom_norm(FTSState x, FTSHomogeneousWeights* hw);
double fts_hom_weighted_norm(FTSState x, FTSHomogeneousWeights* hw, double p);
double fts_hom_euclidean_norm(FTSState x, FTSHomogeneousWeights* hw);

/* --- Convergence Testing --- */
FTSHomogeneousTest* fts_hom_test_convergence(FTSSystem* sys, FTSState x0,
                                              FTSHomogeneousWeights* hw,
                                              double T, double dt);
void fts_hom_test_free(FTSHomogeneousTest* ht);

/* --- Bounds --- */
double fts_hom_settling_time_bound(FTSHomogeneousWeights* hw,
                                    double deg, FTSState x0);
double fts_hom_settling_time_upper(FTSHomogeneousWeights* hw,
                                    double deg, FTSState x0, double k);

/* --- Utility --- */
void fts_hom_print(FTSHomogeneousWeights* hw);
void fts_hom_test_print(FTSHomogeneousTest* ht);
int  fts_hom_compare_weights(FTSHomogeneousWeights* a,
                              FTSHomogeneousWeights* b);

/* --- Standard Weight Presets --- */
FTSHomogeneousWeights* fts_hom_standard_weights(int dim);
FTSHomogeneousWeights* fts_hom_uniform_weights(int dim, double w);

#endif /* FTS_HOMOGENEOUS_H */
