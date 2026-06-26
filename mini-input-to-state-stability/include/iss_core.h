#ifndef ISS_CORE_H
#define ISS_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * iss_core.h ? ISS Core Types
 *
 * Input-to-State Stability (Sontag, 1989):
 *   |x(t)| <= beta(|x(0)|, t) + gamma(sup|u|)
 *
 * where beta in KL class, gamma in K class.
 *
 * ISS Lyapunov function: V'(x) <= -alpha(|x|) + sigma(|u|)
 *   with alpha, sigma in K_inf class.
 *
 * References:
 *   Sontag (1989) Smooth Stabilization Implies Coprime Factorization
 *   Sontag & Wang (1995) On Characterizations of ISS
 *   Khalil (2002) Nonlinear Systems, Ch.4
 * ============================================================== */

#define ISS_MAX_DIM 10
#define ISS_MAX_POINTS 10000

/* Comparison function classes */
typedef enum { ISS_CLASS_K = 0, ISS_CLASS_KINF = 1, ISS_CLASS_KL = 2 } ISS_ClassType;

/* A K-class function: gamma(0)=0, continuous, strictly increasing */
typedef struct ISS_KFunction {
    double (*eval)(const struct ISS_KFunction* self, double s);
    double (*inverse)(const struct ISS_KFunction* self, double s);
    void* param;
    char* name;
    ISS_ClassType type;
} ISS_KFunction;

/* A K_inf class function: K-class with gamma(s)->inf as s->inf */
typedef ISS_KFunction ISS_KInfFunction;

/* A KL class function: beta(s,t), decreasing in t, beta(.,t) in K */
typedef struct {
    double (*eval)(double s, double t);
    void* param;
    char* name;
} ISS_KLFunction;

/* State vector */
typedef struct {
    double x[ISS_MAX_DIM];
    int dim;
} ISS_State;

/* Nonlinear system: dx/dt = f(x) + g(x)*u */
typedef void (*ISS_Dynamics)(const double* x, int n, const double* u, int m, double* dxdt);

/* ISS Lyapunov function candidate V(x) */
typedef double (*ISS_Lyapunov)(const double* x, int n);

/* Gradient of Lyapunov function dV/dx */
typedef void (*ISS_LyapunovGrad)(const double* x, int n, double* grad);

/* ISS certificate */
typedef struct {
    double alpha_coeff;      /* alpha(s) = alpha_coeff * V(x) */
    double sigma_coeff;       /* sigma(s) = sigma_coeff * |u|^2 */
    double gamma_gain;        /* ISS gain gamma */
    double K_inf_alpha_coeff; /* alpha in K_inf class coefficient */
    bool is_ISS;
    bool is_iISS;             /* integral ISS */
    double max_input_norm;
    double lyapunov_decay_rate;
    int verification_samples;
} ISS_Certificate;

/* ISS system description */
typedef struct {
    char* name;
    int n_states;
    int n_inputs;
    ISS_Dynamics dynamics;
    ISS_Lyapunov lyapunov;
    ISS_LyapunovGrad lyapunov_grad;
    ISS_State initial_state;
    ISS_Certificate certificate;
} ISS_System;

/* Standard K-class functions */
ISS_KFunction iss_k_linear(double coeff);
ISS_KFunction iss_k_power(double coeff, double exponent);
ISS_KFunction iss_k_saturation(double coeff, double limit);
ISS_KInfFunction iss_kinf_linear(double coeff);
ISS_KInfFunction iss_kinf_power(double coeff, double exponent);

/* Standard KL-class functions */
ISS_KLFunction iss_kl_exponential(double lambda);
ISS_KLFunction iss_kl_polynomial(double coeff, double exponent);

/* Core utilities */
double iss_k_eval(const ISS_KFunction* k, double s);
double iss_kl_eval(const ISS_KLFunction* kl, double s, double t);
bool iss_k_is_valid(const ISS_KFunction* k);
void iss_k_free(ISS_KFunction* k);
void iss_state_copy(const ISS_State* src, ISS_State* dst);
double iss_state_norm(const ISS_State* s);
ISS_System* iss_system_create(const char* name, int n, int m);
void iss_system_free(ISS_System* sys);

/* Additional K-class builders */
ISS_KFunction iss_k_cubic(double coeff);
ISS_KFunction iss_k_exponential(double coeff);
ISS_KFunction iss_k_deadzone(double coeff, double threshold);
ISS_KFunction iss_k_sigmoid(double coeff, double midpoint, double steepness);

/* Additional KL-class builders */
ISS_KLFunction iss_kl_double_exponential(double l1, double l2);
ISS_KLFunction iss_kl_product(const ISS_KFunction* k, double lambda);

/* State operations */
void iss_state_set(ISS_State* s, const double* values, int n);
double iss_state_distance(const ISS_State* a, const ISS_State* b);
ISS_State iss_state_add(const ISS_State* a, const ISS_State* b);
ISS_State iss_state_scale(const ISS_State* s, double c);
double iss_state_inner_product(const ISS_State* a, const ISS_State* b);
void iss_state_zero(ISS_State* s);

/* Numerical integration */
typedef void (*ISS_VectorField)(const double* x, int n, double* dxdt);
void iss_rk4_step(ISS_VectorField f, double* x, int n, double dt, void* params);
void iss_rk4_integrate(ISS_VectorField f, double* x, int n, double dt, int steps, void* params);

/* ISS trajectory simulation */
typedef struct {
    double* t;
    double** x;
    int n_states;
    int n_steps;
} ISS_Trajectory;
ISS_Trajectory* iss_trajectory_simulate(ISS_System* sys, const ISS_State* x0, const double* u, double dt, int steps);
void iss_trajectory_free(ISS_Trajectory* traj);
double iss_trajectory_max_norm(const ISS_Trajectory* traj);

/* ISS bound computation: |x(t)| <= beta(|x0|, t) + gamma(||u||) */
double iss_bound_beta(const ISS_KLFunction* beta, double x0_norm, double t);
double iss_bound_gamma(const ISS_KFunction* gamma, double u_norm);
double iss_bound_full(const ISS_KLFunction* beta, const ISS_KFunction* gamma, double x0_norm, double t, double u_norm);

/* System dynamics helpers */
void iss_linear_dynamics(const double* x, int n, const double* u, int m, double* dxdt, const double* A, const double* B);
double iss_system_simulate_and_measure(ISS_System* sys, const ISS_State* x0, double* u_seq, int u_steps, double dt);

#endif








































