/* invariant_sets.h - Invariant Set Computation
 * L2: positively invariant, negatively invariant, omega-limit sets
 * L4: LaSalle invariance theorem - largest invariant set
 * L5: Forward/backward reachability, omega-limit approximation
 * Ref: LaSalle (1960), Bhatia & Szego (1970) Stability Theory */

#ifndef INVARIANT_SETS_H
#define INVARIANT_SETS_H
#include "gst_core.h"
#include <stdbool.h>

typedef struct {
    double* points;
    int n_points, capacity, dim;
    bool is_positively_invariant;
    bool is_negatively_invariant;
    bool is_invariant;
    double diameter;
    double volume_estimate;
} InvariantSet;

typedef struct {
    double* omega_points;
    int n_omega, cap_omega, dim;
    double* distances_to_limit;
    int n_distances;
    double approach_rate;
    bool is_singleton;
    bool is_limit_cycle;
} OmegaLimitSet;

InvariantSet* iset_create(int dim);
void iset_free(InvariantSet* is);
int iset_add_point(InvariantSet* is, const double* x);
bool iset_check_positive_invariance(InvariantSet* is, ODEFunc f,
    double dt, double duration, double eps);
bool iset_check_negative_invariance(InvariantSet* is, ODEFunc f,
    double dt, double duration, double eps);
bool iset_find_fixed_points(const InvariantSet* is, ODEFunc f,
    double* fp, int max_fp, int* n_found, double eps);
int iset_grid_compute(InvariantSet* is, ODEFunc f,
    const double* region_min, const double* region_max,
    int grid_n, double dt, double dur, double eps);
double iset_diameter(const InvariantSet* is);
bool iset_is_subset(const InvariantSet* a, const InvariantSet* b,
                     double tol);
void iset_print(const InvariantSet* is);

OmegaLimitSet* omega_create(int dim);
void omega_free(OmegaLimitSet* ol);
int omega_compute(OmegaLimitSet* ol, ODEFunc f, const double* x0,
    int dim, double dt, double dur, int num_samples, double eps);
double omega_convergence_rate(const OmegaLimitSet* ol);
bool omega_is_singleton(const OmegaLimitSet* ol, double tol);
bool omega_contains(const OmegaLimitSet* ol, const double* x, double tol);
void omega_print(const OmegaLimitSet* ol);

#endif
