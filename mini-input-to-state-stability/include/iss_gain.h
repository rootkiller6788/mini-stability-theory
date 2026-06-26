#ifndef ISS_GAIN_H
#define ISS_GAIN_H
#include "iss_core.h"
#include "iss_lyapunov.h"

/* ==============================================================
 * iss_gain.h ? ISS Gain Computation
 *
 * The ISS gain gamma determines how input affects the state:
 *   |x(t)| <= max{ beta(|x(0)|,t), gamma(sup_{tau<=t}|u(tau)|) }
 *
 * Gain properties:
 *   - Zero input -> asymptotic stability
 *   - Bounded input -> bounded state
 *   - Converging input -> converging state
 * ============================================================== */

/* ISS gain function: maps input norm to state bound */
typedef struct {
    ISS_KFunction k_function;
    double linear_gain;       /* gamma(s) = K * s (linear approximation) */
    double peak_gain;
    double l2_gain;           /* L2-to-L2 gain estimate */
    bool is_finite;
    double margin;            /* Stability margin */
} ISS_Gain;

/* Compute ISS gain from Lyapunov function */
ISS_Gain iss_gain_from_lyapunov(const ISS_Certificate* cert);

/* Compute ISS gain via simulation */
ISS_Gain iss_gain_from_simulation(ISS_System* sys, double t_final, double dt,
    double* input_amplitudes, int n_inputs, double state_threshold);

/* Linear gain estimate from Jacobian */
double iss_gain_linear_estimate(ISS_System* sys, const ISS_State* equilibrium);

/* Check small-gain condition: gamma1 * gamma2 < 1 */
bool iss_gain_small_gain_check(const ISS_Gain* g1, const ISS_Gain* g2);

/* Compute H_infinity norm estimate as ISS gain */
double iss_gain_hinf_norm(ISS_System* sys, double omega_min, double omega_max, int n_freq);

/* Print gain analysis */
void iss_gain_print(const ISS_Gain* gain);
void iss_gain_free(ISS_Gain* gain);

/* Extended gain computation */
double iss_gain_asymptotic_estimate(ISS_System* s, double t_final, double dt);
ISS_Gain iss_gain_composite(const ISS_Gain* g1, const ISS_Gain* g2);
double iss_gain_sensitivity(const ISS_Gain* g, double delta);
ISS_Gain iss_gain_monte_carlo(ISS_System* sys, int n_trials, double t_final, double dt, double u_max);
double iss_gain_frequency_response(ISS_System* sys, double omega, double dt);
double iss_gain_margin(const ISS_Gain* g);
double iss_gain_lti(const double* A, const double* B, const double* C, int n, int m, int p);

#endif




























































