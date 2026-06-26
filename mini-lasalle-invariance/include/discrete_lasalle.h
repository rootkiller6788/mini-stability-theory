/* discrete_lasalle.h - Discrete-Time LaSalle Invariance Principle
 * L2: Discrete dynamical systems, forward/backward invariance
 * L4: Discrete LaSalle theorem: V(f(x))-V(x)<=0 => convergence to
 *     largest invariant set in {x: V(f(x))=V(x)}
 * L5: Discrete Lyapunov function construction, invariant set iteration
 * L8: Discrete Barbashin-Krasovskii, hybrid system dwell-time conditions
 * Ref: LaSalle (1976) The Stability of Dynamical Systems
 *      Agarwal (2000) Difference Equations and Inequalities
 *      Goebel, Sanfelice & Teel (2012) Hybrid Dynamical Systems
 */

#ifndef DISCRETE_LASALLE_H
#define DISCRETE_LASALLE_H
#include "gst_core.h"
#include <stdbool.h>

/* Discrete map type: x_{k+1} = F(x_k) */
typedef void (*DiscreteMap)(const double* x, double* x_next, int n);

/* ---------- Discrete Invariant Set ---------- */
typedef struct {
    double* points;
    int n_points, capacity, dim;
    bool is_positively_invariant;
    bool is_negatively_invariant;
    bool is_invariant;
    double diameter;
    int iteration_depth;
} DiscreteInvariantSet;

/* ---------- Discrete LaSalle Result ---------- */
typedef struct {
    int dim;
    DiscreteInvariantSet* invariant_set;
    double* zero_difference_set;
    int n_zero_diff, cap_zero_diff;
    double* convergence_points;
    int n_converged;
    double convergence_rate;
    bool is_asymptotically_stable;
    int iterations_to_converge;
} DiscreteLasalleResult;

/* ---------- Discrete-Time Lyapunov ---------- */
typedef double (*DiscreteLyapunovFunc)(const double* x, int n);

/* ---------- Core Functions ---------- */
DiscreteInvariantSet* diset_create(int dim);
void diset_free(DiscreteInvariantSet* ds);
int diset_add_point(DiscreteInvariantSet* ds, const double* x);
int diset_compute_invariant(DiscreteInvariantSet* ds, DiscreteMap F,
    int max_iter, double eps);
bool diset_check_forward_invariance(const DiscreteInvariantSet* ds,
    DiscreteMap F, double eps);
bool diset_check_backward_invariance(const DiscreteInvariantSet* ds,
    DiscreteMap F, double eps);
int diset_largest_positively_invariant(DiscreteInvariantSet* ds,
    DiscreteMap F, const double* candidate_set, int n_candidates,
    int n_dims, int max_iter, double eps);
void diset_print(const DiscreteInvariantSet* ds);

/* ---------- Discrete LaSalle ---------- */
DiscreteLasalleResult* dlasalle_create(int dim);
void dlasalle_free(DiscreteLasalleResult* dlr);
int dlasalle_find_zero_difference_set(DiscreteLasalleResult* dlr,
    DiscreteLyapunovFunc V, DiscreteMap F,
    const double* samples, int n_samples, int n_dims, double eps);
int dlasalle_find_largest_invariant(DiscreteLasalleResult* dlr,
    DiscreteMap F, int max_iter, double eps);
int dlasalle_verify_conditions(DiscreteLyapunovFunc V, DiscreteMap F,
    const double* x0, int n, int n_iter, double eps);
bool dlasalle_is_globally_asymptotically_stable(
    const DiscreteLasalleResult* dlr, DiscreteLyapunovFunc V,
    double radial_bound, int n_samples);
void dlasalle_print(const DiscreteLasalleResult* dlr);

/* ---------- Discrete Barbashin-Krasovskii ---------- */
bool discrete_bk_verify(DiscreteLyapunovFunc V, DiscreteMap F,
    int n, const double* region_min, const double* region_max,
    int grid_n, double eps);
int discrete_bk_find_equilibrium(DiscreteMap F, double* x0, int n,
    int max_iter, double tol);

/* ---------- Hybrid System Dwell-Time ---------- */
typedef struct {
    double* flow_map_params;
    double* jump_map_params;
    double* guard_set;
    double* flow_set;
    int n_params;
    double min_dwell_time;
    double max_dwell_time;
} HybridDwellTime;

int hybrid_dwell_time_estimate(const HybridDwellTime* hdt,
    DiscreteLyapunovFunc V, DiscreteMap F_flow, DiscreteMap F_jump,
    int n, double eps, double* dwell_time_bound);

#endif