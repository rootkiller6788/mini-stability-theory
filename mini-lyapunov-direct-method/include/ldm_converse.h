#ifndef LDM_CONVERSE_H
#define LDM_CONVERSE_H
#include "ldm_core.h"

/* ==============================================================
 * ldm_converse.h - Converse Lyapunov Theorems
 *
 * Massera (1949): For asymptotically stable systems, a smooth
 * Lyapunov function always exists (converse theorem).
 *
 * Construction methods:
 *   1. Integral: V(x) = int_0^inf ||phi(t,x)||^2 dt
 *   2. Yoshizawa: V(x) = sup_{t>=0} ||phi(t,x)|| * g(t)
 *   3. Kurzweil: V(x) = int_0^inf w(||phi(t,x)||) dt
 *
 * Practical construction via trajectory sampling.
 *
 * References:
 *   Massera (1949) Annals of Math
 *   Yoshizawa (1966) Stability Theory by Liapunov's Second Method
 * ============================================================== */

typedef struct {
    double* V_grid;
    double* dV_grid;
    int nx, ny;
    double x_min, x_max, y_min, y_max;
    bool is_valid;
} LDM_ConverseResult;

/* Construction methods */
LDM_ConverseResult ldm_converse_integral(const LDM_System* sys, int n, double t_max, double dt, double* grid_params);
LDM_ConverseResult ldm_converse_sampling(const LDM_System* sys, int n, int n_samples, double range);
LDM_ConverseResult ldm_converse_kurzweil(const LDM_System* sys, int n, double t_max, double dt, double p);

/* Lifecycle */
void ldm_converse_free(LDM_ConverseResult* cr);
void ldm_converse_print(const LDM_ConverseResult* cr);

/* Validation and analysis */
double ldm_converse_accuracy(const LDM_ConverseResult* cr);
double ldm_converse_decay_rate(const LDM_ConverseResult* cr, int ix, int iy);
bool ldm_converse_is_positive_definite(const LDM_ConverseResult* cr);

/* Analysis utilities */
double ldm_converse_max_value(const LDM_ConverseResult* cr);
double ldm_converse_min_value(const LDM_ConverseResult* cr);
double ldm_converse_average_dV(const LDM_ConverseResult* cr);
double ldm_converse_monotonicity(const LDM_ConverseResult* cr);
double ldm_converse_lipschitz_estimate(const LDM_ConverseResult* cr);
double ldm_converse_unit_circle_norm(const LDM_ConverseResult* cr);

/* Theoretical conditions */
bool ldm_massera_condition(const LDM_System* sys, int n);
double ldm_yoshizawa_bound(const LDM_System* sys, int n, double t_max);

#endif
