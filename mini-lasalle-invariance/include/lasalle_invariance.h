#ifndef LASALLE_INVARIANCE_H
#define LASALLE_INVARIANCE_H
#include "gst_core.h"
#include <stdbool.h>

/* LaSalle's Invariance Principle (1960): For x'=f(x) with V continuously
 * differentiable, if V'(x) <= 0 on a compact positively invariant set Omega,
 * every trajectory starting in Omega approaches the largest invariant set
 * M contained in E = {x in Omega : V'(x) = 0}. */

typedef struct {
    double* omega_points;
    int n_omega, cap_omega;
    double* zero_derivative_set;
    int n_zero, cap_zero;
    double* largest_invariant;
    int n_invariant, cap_invariant;
    double convergence_time;
    double convergence_tolerance;
    bool is_converged;
} LasalleResult;

typedef struct {
    double* x_trajectory;
    int n_steps, n_dims;
    double* V_values;
    double* Vdot_values;
    double final_V, final_Vdot;
    bool approaches_invariant;
    double distance_to_invariant;
} InvarianceAnalysis;

LasalleResult* lasalle_create(void);
void lasalle_free(LasalleResult* lr);
int lasalle_find_zero_derivative_set(LasalleResult* lr, LyapunovFunc V,
    ODEFunc f, const double* x0, int n, double dt,
    double duration, double eps);
int lasalle_find_largest_invariant(LasalleResult* lr, ODEFunc f,
    const double* x0, int n, int max_iter, double dt, double eps);
bool lasalle_verify_conditions(LyapunovFunc V, ODEFunc f,
    const double* x_sample, int n, int n_samples,
    double radius, double eps);
double lasalle_convergence_time_estimate(const LasalleResult* lr);
bool lasalle_is_asymptotically_stable(const LasalleResult* lr);
void lasalle_print(const LasalleResult* lr);

InvarianceAnalysis* invar_analyze(ODEFunc f, LyapunovFunc V,
    const double* x0, int n, double dt, double duration);
void invar_free(InvarianceAnalysis* ia);
double invar_compute_distance(const InvarianceAnalysis* ia,
    const double* invariant_set, int n_invariant, int n_dims);
bool invar_within_tolerance(const InvarianceAnalysis* ia, double tol);
void invar_print(const InvarianceAnalysis* ia);

/* L5: Trajectory energy profile -- V(t) and Vdot(t) evolution */
int lasalle_trajectory_energy_profile(ODEFunc f, LyapunovFunc V,
    const double* x0, int n, double dt, double duration,
    double* V_history, double* Vdot_history, int max_steps);

/* L8: Lasalle-Yoshizawa: non-autonomous system verification */
int lasalle_yoshizawa_check(ODEFunc f, LyapunovFunc V,
    LyapunovFunc W3, const double* samples, int n_samples,
    int n_dims, double eps);

/* L8: Check if a system is a gradient system x_dot = -grad V(x) */
bool lasalle_check_gradient_system(ODEFunc f, int n,
    const double* x_test, double eps);
#endif

