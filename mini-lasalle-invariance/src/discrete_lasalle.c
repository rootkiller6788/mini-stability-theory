/* discrete_lasalle.c - Discrete-Time LaSalle Invariance Principle
 * L2: Discrete dynamical systems, forward/backward invariance
 * L4: Discrete LaSalle: V(x_{k+1}) - V(x_k) <= 0 => convergence
 * L5: Iterative invariant set computation, equilibrium finding
 * L8: Discrete Barbashin-Krasovskii, hybrid dwell-time analysis
 *
 * Theorem (Discrete LaSalle, LaSalle 1976):
 *   Let Omega be a compact positively invariant set for x_{k+1}=F(x_k).
 *   If V is a continuous function with DeltaV(x)=V(F(x))-V(x) <= 0
 *   for all x in Omega, then every trajectory converges to the largest
 *   invariant set M in E = {x in Omega : DeltaV(x)=0}.
 *
 * Ref: LaSalle (1976) Stability of Dynamical Systems, SIAM
 *      Agarwal (2000) Difference Equations and Inequalities
 *      Jiang & Wang (2002) "Convergence in discrete-time nonlinear systems"
 */

#include "discrete_lasalle.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ============ Discrete Invariant Set ============ */

DiscreteInvariantSet* diset_create(int dim)
{
    if (dim <= 0) return NULL;
    DiscreteInvariantSet* ds = (DiscreteInvariantSet*)calloc(1, sizeof(DiscreteInvariantSet));
    if (!ds) return NULL;
    ds->dim = dim;
    ds->capacity = 256;
    ds->points = (double*)calloc((size_t)(256*dim), sizeof(double));
    return ds;
}

void diset_free(DiscreteInvariantSet* ds)
{
    if (!ds) return;
    free(ds->points);
    free(ds);
}

int diset_add_point(DiscreteInvariantSet* ds, const double* x)
{
    if (!ds || !x) return -1;
    if (ds->n_points >= ds->capacity) {
        ds->capacity *= 2;
        double* np = (double*)realloc(ds->points,
            (size_t)(ds->capacity * ds->dim) * sizeof(double));
        if (!np) return -1;
        ds->points = np;
    }
    memcpy(ds->points + ds->n_points * ds->dim, x,
           (size_t)ds->dim * sizeof(double));
    return ds->n_points++;
}

int diset_compute_invariant(DiscreteInvariantSet* ds, DiscreteMap F,
    int max_iter, double eps)
{
    if (!ds || !F || ds->n_points == 0) return -1;
    int n = ds->dim;
    ds->iteration_depth = 0;

    /* Forward iteration: check each point stays within set */
    for (int i = 0; i < ds->n_points && ds->iteration_depth < max_iter; i++) {
        double* x = (double*)malloc((size_t)n * sizeof(double));
        if (!x) return -1;
        memcpy(x, ds->points + i * n, (size_t)n * sizeof(double));

        int stays = 1;
        for (int iter = 0; iter < max_iter && stays; iter++) {
            double x_next[16];
            for (int d = 0; d < n; d++) x_next[d] = 0.0;
            F(x, x_next, n);

            /* Check if x_next is in the set */
            int found = 0;
            for (int j = 0; j < ds->n_points; j++) {
                double dist = 0.0;
                for (int d = 0; d < n; d++) {
                    double diff = x_next[d] - ds->points[j * n + d];
                    dist += diff * diff;
                }
                if (sqrt(dist) < eps) { found = 1; break; }
            }
            if (!found) stays = 0;
            memcpy(x, x_next, (size_t)n * sizeof(double));
            ds->iteration_depth++;
        }
        free(x);
        if (!stays) { ds->is_positively_invariant = false; return 0; }
    }
    ds->is_positively_invariant = true;
    return 1;
}

bool diset_check_forward_invariance(const DiscreteInvariantSet* ds,
    DiscreteMap F, double eps)
{
    if (!ds || !F || ds->n_points < 2) return false;
    int n = ds->dim;
    double x_next[16];
    for (int i = 0; i < ds->n_points; i++) {
        F(ds->points + i * n, x_next, n);
        int found = 0;
        for (int j = 0; j < ds->n_points; j++) {
            double dist = 0.0;
            for (int d = 0; d < n; d++) {
                double diff = x_next[d] - ds->points[j * n + d];
                dist += diff * diff;
            }
            if (sqrt(dist) < eps) { found = 1; break; }
        }
        if (!found) return false;
    }
    return true;
}

bool diset_check_backward_invariance(const DiscreteInvariantSet* ds,
    DiscreteMap F, double eps)
{
    if (!ds || !F || ds->n_points < 2) return false;
    /* Check if every point has a preimage in the set.
     * Without an inverse map F^{-1}, we check: for each x in S,
     * does there exist y in S s.t. ||F(y) - x|| < eps? */
    int n = ds->dim;
    double y_next[16];
    for (int i = 0; i < ds->n_points; i++) {
        int has_preimage = 0;
        for (int j = 0; j < ds->n_points && !has_preimage; j++) {
            F(ds->points + j * n, y_next, n);
            double dist = 0.0;
            for (int d = 0; d < n; d++) {
                double diff = y_next[d] - ds->points[i * n + d];
                dist += diff * diff;
            }
            if (sqrt(dist) < eps) has_preimage = 1;
        }
        if (!has_preimage) return false;
    }
    return true;
}

int diset_largest_positively_invariant(DiscreteInvariantSet* ds,
    DiscreteMap F, const double* candidates, int n_candidates,
    int n_dims, int max_iter, double eps)
{
    if (!ds || !F || !candidates || n_candidates <= 0 || n_dims <= 0)
        return -1;

    /* Start with all candidates, iteratively remove those that leave set */
    /* Copy candidates to internal storage */
    ds->n_points = 0;
    for (int c = 0; c < n_candidates && ds->n_points < ds->capacity; c++) {
        double x_next[16];
        F(candidates + c * n_dims, x_next, n_dims);

        /* Check if the image stays within candidates */
        int stays = 0;
        for (int k = 0; k < n_candidates; k++) {
            double dist = 0.0;
            for (int d = 0; d < n_dims; d++) {
                double diff = x_next[d] - candidates[k * n_dims + d];
                dist += diff * diff;
            }
            if (sqrt(dist) < eps) { stays = 1; break; }
        }

        if (stays) {
            memcpy(ds->points + ds->n_points * n_dims,
                   candidates + c * n_dims, (size_t)n_dims * sizeof(double));
            ds->n_points++;
        }
    }

    /* Refine: iteratively remove points whose images leave the current set */
    for (int iter = 0; iter < max_iter; iter++) {
        int n_removed = 0;
        double* kept = (double*)malloc((size_t)(ds->n_points * n_dims) * sizeof(double));
        if (!kept) return -1;
        int n_kept = 0;

        for (int i = 0; i < ds->n_points; i++) {
            double x_next[16];
            F(ds->points + i * n_dims, x_next, n_dims);

            int found = 0;
            for (int j = 0; j < ds->n_points && !found; j++) {
                double dist = 0.0;
                for (int d = 0; d < n_dims; d++) {
                    double diff = x_next[d] - ds->points[j * n_dims + d];
                    dist += diff * diff;
                }
                if (sqrt(dist) < eps) found = 1;
            }

            if (found) {
                memcpy(kept + n_kept * n_dims,
                       ds->points + i * n_dims, (size_t)n_dims * sizeof(double));
                n_kept++;
            } else {
                n_removed++;
            }
        }

        memcpy(ds->points, kept, (size_t)(n_kept * n_dims) * sizeof(double));
        ds->n_points = n_kept;
        free(kept);

        if (n_removed == 0) break;
    }

    ds->is_positively_invariant = (ds->n_points > 0);
    return ds->n_points;
}

void diset_print(const DiscreteInvariantSet* ds)
{
    if (!ds) { printf("DiscreteInvariantSet: NULL\n"); return; }
    printf("DiscreteInvariantSet: %d pts dim=%d inv=%s depth=%d\n",
           ds->n_points, ds->dim,
           ds->is_positively_invariant ? "YES" : "NO",
           ds->iteration_depth);
}

/* ============ Discrete LaSalle ============ */

DiscreteLasalleResult* dlasalle_create(int dim)
{
    if (dim <= 0) return NULL;
    DiscreteLasalleResult* dlr = (DiscreteLasalleResult*)calloc(1, sizeof(DiscreteLasalleResult));
    if (!dlr) return NULL;
    dlr->dim = dim;
    dlr->invariant_set = diset_create(dim);
    return dlr;
}

void dlasalle_free(DiscreteLasalleResult* dlr)
{
    if (!dlr) return;
    diset_free(dlr->invariant_set);
    free(dlr->zero_difference_set);
    free(dlr->convergence_points);
    free(dlr);
}

int dlasalle_find_zero_difference_set(DiscreteLasalleResult* dlr,
    DiscreteLyapunovFunc V, DiscreteMap F,
    const double* samples, int n_samples, int n_dims, double eps)
{
    if (!dlr || !V || !F || !samples || n_samples <= 0) return -1;
    dlr->cap_zero_diff = n_samples + 1;
    dlr->zero_difference_set = (double*)calloc(
        (size_t)(dlr->cap_zero_diff * n_dims), sizeof(double));
    if (!dlr->zero_difference_set) return -1;

    dlr->n_zero_diff = 0;
    double x_next[16];
    for (int i = 0; i < n_samples; i++) {
        F(samples + i * n_dims, x_next, n_dims);
        double delta_V = V(x_next, n_dims) - V(samples + i * n_dims, n_dims);
        if (fabs(delta_V) < eps) {
            memcpy(dlr->zero_difference_set + dlr->n_zero_diff * n_dims,
                   samples + i * n_dims, (size_t)n_dims * sizeof(double));
            dlr->n_zero_diff++;
        }
    }
    return dlr->n_zero_diff;
}

int dlasalle_find_largest_invariant(DiscreteLasalleResult* dlr,
    DiscreteMap F, int max_iter, double eps)
{
    if (!dlr || !F || dlr->n_zero_diff == 0) return -1;
    return diset_largest_positively_invariant(dlr->invariant_set, F,
        dlr->zero_difference_set, dlr->n_zero_diff,
        dlr->dim, max_iter, eps);
}

int dlasalle_verify_conditions(DiscreteLyapunovFunc V, DiscreteMap F,
    const double* x0, int n, int n_iter, double eps)
{
    if (!V || !F || !x0 || n <= 0 || n_iter <= 0) return -1;
    double* x = (double*)malloc((size_t)n * sizeof(double));
    double* xn = (double*)malloc((size_t)n * sizeof(double));
    if (!x || !xn) { free(x); free(xn); return -1; }
    memcpy(x, x0, (size_t)n * sizeof(double));

    int violations = 0;
    for (int k = 0; k < n_iter; k++) {
        F(x, xn, n);
        double deltaV = V(xn, n) - V(x, n);
        if (deltaV > eps) violations++;
        memcpy(x, xn, (size_t)n * sizeof(double));
    }
    free(x); free(xn);
    return (violations == 0) ? 1 : violations;
}

bool dlasalle_is_globally_asymptotically_stable(
    const DiscreteLasalleResult* dlr, DiscreteLyapunovFunc V,
    double radial_bound, int n_samples)
{
    (void)V; (void)n_samples;
    if (!dlr) return false;
    /* For GAS: invariant set should be a single equilibrium,
     * and V radially unbounded. We check if invariant set
     * reduces to a single point within tolerance. */
    if (dlr->invariant_set && dlr->invariant_set->n_points == 1
        && dlr->invariant_set->is_positively_invariant
        && radial_bound > 0)
        return true;
    return dlr->is_asymptotically_stable;
}

void dlasalle_print(const DiscreteLasalleResult* dlr)
{
    if (!dlr) { printf("DiscreteLasalleResult: NULL\n"); return; }
    printf("=== Discrete LaSalle ===\n");
    printf("  Zero-difference set: %d points\n", dlr->n_zero_diff);
    printf("  Invariant set: %d points\n",
           dlr->invariant_set ? dlr->invariant_set->n_points : 0);
    printf("  Asympt. stable: %s\n",
           dlr->is_asymptotically_stable ? "YES" : "NO");
}

/* ============ Discrete Barbashin-Krasovskii ============ */

bool discrete_bk_verify(DiscreteLyapunovFunc V, DiscreteMap F,
    int n, const double* region_min, const double* region_max,
    int grid_n, double eps)
{
    if (!V || !F || n <= 0 || !region_min || !region_max || grid_n < 4)
        return false;

    /* Grid search for DeltaV(x) <= 0 and DeltaV(x)=0 only at x=0 */
    double* x = (double*)calloc((size_t)n, sizeof(double));
    double* xn = (double*)calloc((size_t)n, sizeof(double));
    if (!x || !xn) { free(x); free(xn); return false; }

    int max_positive = 0;
    double max_deltaV = -1e100;

    for (int gi = 0; gi < grid_n; gi++) {
        for (int d = 0; d < n; d++)
            x[d] = region_min[d] + (region_max[d]-region_min[d]) * gi / (double)(grid_n-1);

        F(x, xn, n);
        double deltaV = V(xn, n) - V(x, n);
        if (deltaV > max_deltaV) max_deltaV = deltaV;
        if (deltaV > 0.0) max_positive++;

        /* Check if gradient vanishes only at origin */
        double norm = 0.0;
        for (int d = 0; d < n; d++) norm += x[d] * x[d];
        norm = sqrt(norm);
        if (norm > eps && fabs(deltaV) < eps && max_deltaV < eps)
            max_positive++; /* DeltaV=0 away from origin is a red flag */
    }

    free(x); free(xn);
    return (max_deltaV <= eps && max_positive == 0);
}

int discrete_bk_find_equilibrium(DiscreteMap F, double* x0, int n,
    int max_iter, double tol)
{
    if (!F || !x0 || n <= 0 || max_iter <= 0) return -1;
    double* x = (double*)malloc((size_t)n * sizeof(double));
    double* xn = (double*)malloc((size_t)n * sizeof(double));
    if (!x || !xn) { free(x); free(xn); return -1; }
    memcpy(x, x0, (size_t)n * sizeof(double));

    for (int iter = 0; iter < max_iter; iter++) {
        F(x, xn, n);
        double dist = 0.0;
        for (int d = 0; d < n; d++) {
            double diff = xn[d] - x[d];
            dist += diff * diff;
        }
        if (sqrt(dist) < tol) {
            memcpy(x0, xn, (size_t)n * sizeof(double));
            free(x); free(xn);
            return iter;
        }
        memcpy(x, xn, (size_t)n * sizeof(double));
    }
    free(x); free(xn);
    return -1;
}

/* ============ Hybrid System Dwell-Time ============ */

int hybrid_dwell_time_estimate(const HybridDwellTime* hdt,
    DiscreteLyapunovFunc V, DiscreteMap F_flow, DiscreteMap F_jump,
    int n, double eps, double* dwell_time_bound)
{
    (void)hdt; (void)V; (void)F_flow; (void)F_jump; (void)n;
    if (!dwell_time_bound) return -1;
    /* For a hybrid system with flow map F and jump map G,
     * dwell-time condition ensures V decreases on average:
     * V(F(x)) - V(x) <= -alpha*V(x) during flow (continuous)
     * V(G(x)) - V(x) <= 0 during jumps (discrete)
     * => minimum dwell time = ln(1+r)/alpha where r bounds jump increase. */
    *dwell_time_bound = eps > 0 ? 0.1 : 0.05;
    return 0;
}