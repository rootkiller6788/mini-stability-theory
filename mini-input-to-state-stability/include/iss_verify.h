#ifndef ISS_VERIFY_H
#define ISS_VERIFY_H
#include "iss_core.h"
#include "iss_lyapunov.h"
#include "iss_gain.h"

/* ==============================================================
 * iss_verify.h ? Numerical ISS Verification Suite
 *
 * Full ISS verification pipeline:
 *   1. Lyapunov function candidate check
 *   2. ISS condition verification
 *   3. Gain estimation
 *   4. Small-gain analysis (for interconnected)
 * ============================================================== */

/* Verification configuration */
typedef struct {
    double state_radius;
    double input_radius;
    int grid_resolution;
    int monte_carlo_samples;
    double tolerance;
    bool use_adaptive_sampling;
    bool check_radial_unboundedness;
} ISS_VerifyConfig;

/* Full verification result */
typedef struct {
    bool is_ISS;
    ISS_LyapunovResult lyap_result;
    ISS_Gain gain;
    ISS_Certificate certificate;
    double verification_time_ms;
    int total_checks;
    int passed_checks;
} ISS_Verification;

/* Default configuration */
ISS_VerifyConfig iss_verify_config_default(void);

/* Run full ISS verification */
ISS_Verification iss_verify(ISS_System* sys, const ISS_VerifyConfig* config);

/* Quick ISS check */
bool iss_verify_quick(ISS_System* sys);

/* Verify K_inf property of a function */
bool iss_verify_kinf(double (*func)(double), double tol,
    int n_samples, double max_s);

/* Verify KL property of a function */
bool iss_verify_kl(double (*func)(double, double), double tol,
    int n_s, int n_t, double s_max, double t_max);

/* Print verification result */
void iss_verification_print(const ISS_Verification* ver);

/* Convergence test: check if trajectory approaches equilibrium */
bool iss_verify_convergence(ISS_System* sys, double t_final, double dt,
    int n_initial_conditions);

/* Extended verification */
double iss_verify_stability_margin(ISS_System* s);
int iss_verify_pass_count(const ISS_Verification* v);
bool iss_verify_is_definitive(const ISS_Verification* v);
ISS_Verification iss_verify_adaptive(ISS_System* sys, const ISS_VerifyConfig* config);
ISS_Certificate iss_verify_full_pipeline(ISS_System* sys, const ISS_VerifyConfig* config);
int iss_verify_batch(ISS_System** systems, int n_systems, ISS_Verification* results);
bool iss_verify_robustness(ISS_System* sys, double param_delta, int n_perturbations);
double iss_verify_gain_empirical(ISS_System* sys, double dt, double t_final, double u_magnitude);

#endif




























































