#ifndef ISS_SMALL_GAIN_H
#define ISS_SMALL_GAIN_H
#include "iss_core.h"
#include "iss_gain.h"

/* ==============================================================
 * iss_small_gain.h ? Small-Gain Theorem for ISS
 *
 * For interconnected systems:
 *   dx1/dt = f1(x1,x2,u1), ISS with gain gamma1
 *   dx2/dt = f2(x2,x1,u2), ISS with gain gamma2
 *
 * If gamma1(gamma2(s)) < s for all s > 0, the interconnection is ISS.
 *
 * Cyclic small-gain for N interconnected systems.
 * ============================================================== */

/* Interconnected system descriptor */
typedef struct {
    ISS_System** subsystems;
    int n_systems;
    ISS_Gain** gains;          /* n x n gain matrix: gamma_{ij} */
    int* n_gains;
} ISS_Interconnection;

/* Small-gain condition check */
bool iss_small_gain_condition(const ISS_Gain* g1, const ISS_Gain* g2);

/* Cyclic small-gain for N systems */
bool iss_cyclic_small_gain(const ISS_Gain** gains, int n);

/* Monotone small-gain (for MIMO interconnections) */
bool iss_monotone_small_gain(const ISS_Interconnection* ic);

/* Construct interconnection */
ISS_Interconnection* iss_interconnection_create(int n_systems);
void iss_interconnection_free(ISS_Interconnection* ic);
void iss_interconnection_set_gain(ISS_Interconnection* ic, int i, int j, const ISS_Gain* g);

/* Check if interconnection is ISS */
bool iss_interconnection_is_ISS(const ISS_Interconnection* ic);

/* Max-min small gain theorem */
bool iss_max_min_small_gain(const ISS_Interconnection* ic, double tol);

/* Extended small-gain */
double iss_small_gain_margin(const ISS_Gain** gains, int n);
bool iss_small_gain_is_tight(const ISS_Gain** gains, int n, double tol);
int iss_small_gain_violation_index(const ISS_Gain** gains, int n);
bool iss_nonlinear_small_gain(const ISS_KFunction* g1, const ISS_KFunction* g2,
    double s_max, int n_points);
double iss_spectral_radius(const double* M, int n);
bool iss_small_gain_lyapunov(const ISS_System* s1, const ISS_System* s2,
    double c1, double c2, double tol);
void iss_build_gain_matrix(ISS_Interconnection* ic, const double* gain_values);

#endif




























































