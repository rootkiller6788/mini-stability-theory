/* region_attraction.h - Region of Attraction Estimation
 * L4: Lyapunov ROA theorem - sublevel sets of V estimate ROA
 * L5: Bisection ROA, Monte Carlo verification, SOS relaxation
 * L6: Canonical: Van der Pol oscillator ROA, pendulum ROA
 * Ref: Chesi (2011) Domain of Attraction, Khalil (2002) Ch.8 */

#ifndef REGION_ATTRACTION_H
#define REGION_ATTRACTION_H
#include "gst_core.h"
#include <stdbool.h>

typedef struct {
    double* sublevel_points;
    int n_points, cap_points, dim;
    double lyapunov_value;
    double max_safe_level;
    double* boundary_points;
    int n_boundary;
    double volume_estimate;
} SublevelSet;

typedef struct {
    int dim;
    SublevelSet* inner_approximation;
    double* roa_boundary;
    int n_boundary_points;
    double estimated_volume;
    double confidence;
    double* samples;
    double* conv_flags;
    int n_samples, n_conv;
    char method_name[64];
} RegionAttraction;

SublevelSet* sublevel_create(int dim);
void sublevel_free(SublevelSet* ss);
int sublevel_compute(SublevelSet* ss, LyapunovFunc V, int dim,
    const double* center, const double* bbox_min,
    const double* bbox_max, int grid_n, double level);
int sublevel_find_max_safe(SublevelSet* ss, LyapunovFunc V,
    ODEFunc f, int dim, const double* center,
    double max_level, int bisection_iters, double eps);
void sublevel_print(const SublevelSet* ss);

RegionAttraction* roa_create(int dim);
void roa_free(RegionAttraction* roa);
int roa_estimate_lyapunov(RegionAttraction* roa, LyapunovFunc V,
    ODEFunc f, int dim, const double* equilibrium,
    double max_radius, int grid_n, double eps);
int roa_monte_carlo_verify(RegionAttraction* roa, ODEFunc f,
    int dim, const double* center, double radius,
    int num_samples, double dt, double dur, double tol);
int roa_sampling_grow(RegionAttraction* roa, ODEFunc f, int dim,
    const double* center, int max_iters, double dt, double tol);
bool roa_contains(const RegionAttraction* roa, const double* x);
void roa_print(const RegionAttraction* roa);

#endif
