#ifndef LDM_REGION_H
#define LDM_REGION_H
#include "ldm_core.h"

/* ==============================================================
 * ldm_region.h - Region of Attraction Estimation
 *
 * For asymptotically stable equilibria, the region of attraction
 * (RoA) is the set of initial conditions that converge to it.
 *
 * Lyapunov-based estimation:
 *   RoA contains {x : V(x) <= c} where c = min V(x) on boundary
 *   of region where dV/dt < 0.
 *
 * Methods:
 *   1. Sublevel set: V(x) <= c_max = max{c : dV/dt<0 in V<=c}
 *   2. Sampling-based: Monte Carlo simulation from grid
 *   3. Zubov method: solve V' = -phi(x)*(1-V) for RoA boundary
 *   4. Chetaev: identify instability regions
 *
 * References:
 *   Chiang et al. (1988) IEEE TAC - RoA estimation
 *   Zubov (1964) Methods of A.M. Lyapunov
 * ============================================================== */

typedef struct {
    double* P;
    int n;
    double c_max;
    double volume_estimate;
    double* boundary_points;
    int n_boundary;
    bool is_conservative;
} LDM_RegionEstimate;

/* Sublevel set methods */
LDM_RegionEstimate ldm_roa_sublevel(const LDM_System* sys, const double* P, int n, double c_try);
double ldm_roa_max_c(const LDM_System* sys, const double* P, int n, double c_min, double c_max, int n_grid);

/* Monte Carlo methods */
LDM_RegionEstimate ldm_roa_sampling(const LDM_System* sys, int n, int n_samples, double range, double t_final);
LDM_RegionEstimate ldm_roa_sampling_filtered(const LDM_System* sys, int n, int n_samples, double range, double t_final, const double* P, double c_max);

/* Zubov's method */
LDM_RegionEstimate ldm_roa_zubov(const LDM_System* sys, int n, double c_level, int n_grid);

/* Lifecycle */
void ldm_region_free(LDM_RegionEstimate* re);
void ldm_region_print(const LDM_RegionEstimate* re);

/* Utilities */
double ldm_roa_volume_ratio(const LDM_RegionEstimate* re, double total_volume);
bool ldm_point_in_roa(const LDM_System* sys, const double* x0, int n, double t_final, double tol);
double ldm_roa_ellipsoid_volume(const double* P, int n, double c);
double ldm_roa_consistency(const LDM_RegionEstimate* re1, const LDM_RegionEstimate* re2);

/* Instability */
double ldm_roa_chetaev_boundary(const LDM_System* sys, int n);

#endif
