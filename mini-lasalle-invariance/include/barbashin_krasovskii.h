/* barbashin_krasovskii.h - Barbashin-Krasovskii Theorem
 * L4: BK Theorem - if V(x) -> 0 as ||x|| -> 0, Vdot <= 0,
 *     and the largest invariant set in {Vdot=0} is {0},
 *     then the origin is globally asymptotically stable.
 * Ref: Barbashin & Krasovskii (1952), Khalil (2002) Thm 4.2 */

#ifndef BARBASHIN_KRASOVSKII_H
#define BARBASHIN_KRASOVSKII_H
#include "gst_core.h"
#include <stdbool.h>

typedef struct {
    double* candidate_points;
    int n_candidates, n_dim;
    double* verified_equilibria;
    int n_verified;
    double verification_radius;
    bool is_globally_asymptotically_stable;
    double lyapunov_margin;
} BKResult;

BKResult* bk_create(int dim);
void bk_free(BKResult* res);
int bk_verify_radial_unbounded(LyapunovFunc V, int n,
                                double radius, int samples);
int bk_verify_negative_definite(LyapunovFunc V, ODEFunc f,
    int n, const double* region_min, const double* region_max,
    int grid_points, double eps);
int bk_find_equilibrium_set(BKResult* res, ODEFunc f,
    const double* x0, int n, int max_iter, double tol);
bool bk_is_gas(const BKResult* res);
void bk_print(const BKResult* res);

#endif
