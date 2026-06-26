/* limit_cycle_analysis.h - Limit Cycle Detection and Orbital Stability
 * L4: Poincare-Bendixson theorem implications for planar systems
 * L5: Poincare map computation, Floquet multiplier analysis
 * L6: Van der Pol, FitzHugh-Nagumo, BZ reaction limit cycles
 * L8: Monodromy matrix, orbital stability via Floquet theory
 * Ref: Guckenheimer & Holmes (1983) Nonlinear Oscillations
 *      Nayfeh & Balachandran (1995) Applied Nonlinear Dynamics
 */

#ifndef LIMIT_CYCLE_ANALYSIS_H
#define LIMIT_CYCLE_ANALYSIS_H
#include "gst_core.h"
#include <stdbool.h>

/* ---------- Poincare Section ---------- */
typedef struct {
    int dim;
    double* section_plane_normal;
    double section_plane_offset;
    double* intersection_points;
    int n_intersections, cap_intersections;
    double* intersection_times;
    double period_estimate;
    double period_variance;
} PoincareSection;

/* ---------- Limit Cycle Result ---------- */
typedef struct {
    int dim;
    double* cycle_points;
    int n_cycle_points, cap_cycle_points;
    double period;
    double* floquet_multipliers;
    int n_multipliers;
    double orbital_stability_margin;
    bool is_stable;
    bool is_orbitally_stable;
    double lyapunov_exponent;
    double amplitude;
    char system_name[64];
} LimitCycleResult;

/* ---------- Poincare Map ---------- */
PoincareSection* poincare_create(int dim, const double* plane_normal,
                                  double offset);
void poincare_free(PoincareSection* ps);
int poincare_compute_intersections(PoincareSection* ps, ODEFunc f,
    const double* x0, int n, double dt, double t_max, double eps);
double poincare_period_estimate(const PoincareSection* ps);
int poincare_map_linearize(const PoincareSection* ps, ODEFunc f,
    int n, double dt, double eps, double* monodromy_matrix);

/* ---------- Limit Cycle Detection ---------- */
LimitCycleResult* lc_result_create(int dim);
void lc_result_free(LimitCycleResult* lc);
int lc_detect_from_trajectory(LimitCycleResult* lc,
    const double* traj, int n_steps, int n_dims,
    double dt, double eps);
int lc_detect_poincare(LimitCycleResult* lc, ODEFunc f,
    const double* x0, int n, double dt, double t_max, double eps);
int lc_compute_floquet(LimitCycleResult* lc, ODEFunc f,
    int n, double dt, double eps);
bool lc_is_orbitally_stable(const LimitCycleResult* lc);
double lc_phase_response(const LimitCycleResult* lc, double phase);
int lc_period_estimate_from_fft(const LimitCycleResult* lc,
    double dt, double* dominant_freq);
void lc_print(const LimitCycleResult* lc);

/* ---------- Bendixson-Dulac Criterion ---------- */
bool bendixson_dulac_criterion_2d(ODEFunc f, const double* region_min,
    const double* region_max, int grid_n, double eps);
int bendixson_find_limit_cycle_region_2d(ODEFunc f,
    const double* bbox_min, const double* bbox_max,
    int grid_n, double* inner_radius, double* outer_radius);

/* ---------- FitzHugh-Nagumo Neural Excitation ---------- */
void fitzhugh_nagumo(const double* x, double* dxdt, int n, double t);
void fitzhugh_nagumo_param(const double* x, double* dxdt, int n, double t,
    double I_ext, double a, double b, double tau);

#endif