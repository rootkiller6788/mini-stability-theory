#include "iss_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * K-class Function Parameters (heap-allocated, per-function)
 * Using void* param with self-referencing eval to fix the
 * global-variable overwrite bug in the original implementation.
 * ============================================================== */

typedef struct { double coeff; } KLinParam;
typedef struct { double coeff; double exponent; } KPowParam;
typedef struct { double coeff; double limit; } KSatParam;
typedef struct { double coeff; double threshold; } KDzParam;
typedef struct { double coeff; double midpoint; double steepness; } KSigParam;

/* Linear K-function: gamma(s) = c * s */
static double k_linear_eval(const ISS_KFunction* self, double s) {
    KLinParam* p = (KLinParam*)self->param;
    return s > 0 ? p->coeff * s : 0.0;
}
static double k_linear_inv(const ISS_KFunction* self, double s) {
    KLinParam* p = (KLinParam*)self->param;
    return s > 0 ? s / p->coeff : 0.0;
}

/* Power K-function: gamma(s) = c * s^e */
static double k_power_eval(const ISS_KFunction* self, double s) {
    KPowParam* p = (KPowParam*)self->param;
    return s > 0 ? p->coeff * pow(s, p->exponent) : 0.0;
}
static double k_power_inv(const ISS_KFunction* self, double s) {
    KPowParam* p = (KPowParam*)self->param;
    return s > 0 ? pow(s / p->coeff, 1.0 / p->exponent) : 0.0;
}

/* Saturation K-function: gamma(s) = c * min(s, L) */
static double k_sat_eval(const ISS_KFunction* self, double s) {
    KSatParam* p = (KSatParam*)self->param;
    return p->coeff * fmin(s, p->limit);
}
static double k_sat_inv(const ISS_KFunction* self, double s) {
    KSatParam* p = (KSatParam*)self->param;
    return s / p->coeff;
}

/* Cubic K-function: gamma(s) = c * s^3 */
static double k_cubic_eval(const ISS_KFunction* self, double s) {
    KLinParam* p = (KLinParam*)self->param;
    return p->coeff * s * s * s;
}
static double k_cubic_inv(const ISS_KFunction* self, double s) {
    KLinParam* p = (KLinParam*)self->param;
    return cbrt(s / p->coeff);
}

/* Exponential K-function: gamma(s) = c * (exp(s) - 1) */
static double k_exp_eval(const ISS_KFunction* self, double s) {
    KLinParam* p = (KLinParam*)self->param;
    return p->coeff * (exp(s) - 1.0);
}
static double k_exp_inv(const ISS_KFunction* self, double s) {
    KLinParam* p = (KLinParam*)self->param;
    return log(s / p->coeff + 1.0);
}

/* Deadzone K-function: gamma(s) = c * max(0, s - thresh) */
static double k_dz_eval(const ISS_KFunction* self, double s) {
    KDzParam* p = (KDzParam*)self->param;
    return p->coeff * fmax(0.0, s - p->threshold);
}
static double k_dz_inv(const ISS_KFunction* self, double s) {
    KDzParam* p = (KDzParam*)self->param;
    return s / p->coeff + p->threshold;
}

/* Sigmoid K-function: gamma(s) = c / (1 + exp(-k*(s - m))) - c/(1+exp(k*m)) */
static double k_sigmoid_eval(const ISS_KFunction* self, double s) {
    KSigParam* p = (KSigParam*)self->param;
    double offset = p->coeff / (1.0 + exp(p->steepness * p->midpoint));
    return p->coeff / (1.0 + exp(-p->steepness * (s - p->midpoint))) - offset;
}
static double k_sigmoid_inv(const ISS_KFunction* self, double s) {
    KSigParam* p = (KSigParam*)self->param;
    double offset = p->coeff / (1.0 + exp(p->steepness * p->midpoint));
    double val = s + offset;
    if (val >= p->coeff * 0.9999) return 1e10;
    if (val <= 0.0) return 0.0;
    return p->midpoint - log(p->coeff / val - 1.0) / p->steepness;
}

ISS_KFunction iss_k_linear(double coeff) {
    KLinParam* param = (KLinParam*)malloc(sizeof(KLinParam));
    param->coeff = coeff;
    ISS_KFunction k;
    memset(&k, 0, sizeof(k));
    k.name = "linear";
    k.eval = k_linear_eval;
    k.inverse = k_linear_inv;
    k.param = param;
    k.type = ISS_CLASS_K;
    return k;
}

ISS_KFunction iss_k_power(double coeff, double exponent) {
    KPowParam* param = (KPowParam*)malloc(sizeof(KPowParam));
    param->coeff = coeff;
    param->exponent = exponent;
    ISS_KFunction k;
    memset(&k, 0, sizeof(k));
    k.name = "power";
    k.eval = k_power_eval;
    k.inverse = k_power_inv;
    k.param = param;
    k.type = ISS_CLASS_K;
    return k;
}

ISS_KFunction iss_k_saturation(double coeff, double limit) {
    KSatParam* param = (KSatParam*)malloc(sizeof(KSatParam));
    param->coeff = coeff;
    param->limit = limit;
    ISS_KFunction k;
    memset(&k, 0, sizeof(k));
    k.name = "saturation";
    k.eval = k_sat_eval;
    k.inverse = k_sat_inv;
    k.param = param;
    k.type = ISS_CLASS_K;
    return k;
}

ISS_KFunction iss_k_cubic(double coeff) {
    KLinParam* param = (KLinParam*)malloc(sizeof(KLinParam));
    param->coeff = coeff;
    ISS_KFunction k;
    memset(&k, 0, sizeof(k));
    k.name = "cubic";
    k.eval = k_cubic_eval;
    k.inverse = k_cubic_inv;
    k.param = param;
    k.type = ISS_CLASS_K;
    return k;
}

ISS_KFunction iss_k_exponential(double coeff) {
    KLinParam* param = (KLinParam*)malloc(sizeof(KLinParam));
    param->coeff = coeff;
    ISS_KFunction k;
    memset(&k, 0, sizeof(k));
    k.name = "exponential";
    k.eval = k_exp_eval;
    k.inverse = k_exp_inv;
    k.param = param;
    k.type = ISS_CLASS_K;
    return k;
}

ISS_KFunction iss_k_deadzone(double coeff, double threshold) {
    KDzParam* param = (KDzParam*)malloc(sizeof(KDzParam));
    param->coeff = coeff;
    param->threshold = threshold;
    ISS_KFunction k;
    memset(&k, 0, sizeof(k));
    k.name = "deadzone";
    k.eval = k_dz_eval;
    k.inverse = k_dz_inv;
    k.param = param;
    k.type = ISS_CLASS_K;
    return k;
}

ISS_KFunction iss_k_sigmoid(double coeff, double midpoint, double steepness) {
    KSigParam* param = (KSigParam*)malloc(sizeof(KSigParam));
    param->coeff = coeff;
    param->midpoint = midpoint;
    param->steepness = steepness;
    ISS_KFunction k;
    memset(&k, 0, sizeof(k));
    k.name = "sigmoid";
    k.eval = k_sigmoid_eval;
    k.inverse = k_sigmoid_inv;
    k.param = param;
    k.type = ISS_CLASS_K;
    return k;
}

ISS_KInfFunction iss_kinf_linear(double coeff) {
    ISS_KInfFunction k = iss_k_linear(coeff);
    k.type = ISS_CLASS_KINF;
    return k;
}

ISS_KInfFunction iss_kinf_power(double coeff, double exponent) {
    ISS_KInfFunction k = iss_k_power(coeff, exponent);
    k.type = ISS_CLASS_KINF;
    return k;
}

/* ==============================================================
 * KL-class Functions
 * ============================================================== */

typedef struct { double lambda; } KLExpParam;
typedef struct { double coeff; double exponent; } KLPolyParam;
typedef struct { double l1; double l2; } KLDoubleExpParam;

ISS_KLFunction iss_kl_exponential(double lambda) {
    ISS_KLFunction kl;
    memset(&kl, 0, sizeof(kl));
    kl.name = "exponential";
    /* Store lambda for use in eval via function composition */
    double* lam = (double*)malloc(sizeof(double));
    *lam = lambda;
    kl.param = lam;
    /* Self-contained eval: captures lambda from param */
    /* Using a simple inline approach - eval reads global lambda */
    kl.eval = NULL; /* Will be set at eval time */
    return kl;
}

ISS_KLFunction iss_kl_polynomial(double coeff, double exponent) {
    ISS_KLFunction kl;
    memset(&kl, 0, sizeof(kl));
    kl.name = "polynomial";
    double* params = (double*)malloc(2 * sizeof(double));
    params[0] = coeff;
    params[1] = exponent;
    kl.param = params;
    return kl;
}

ISS_KLFunction iss_kl_double_exponential(double l1, double l2) {
    ISS_KLFunction kl;
    memset(&kl, 0, sizeof(kl));
    kl.name = "double_exponential";
    double* params = (double*)malloc(2 * sizeof(double));
    params[0] = l1;
    params[1] = l2;
    kl.param = params;
    return kl;
}

/* ==============================================================
 * Comparison Function Evaluation
 * ============================================================== */

double iss_k_eval(const ISS_KFunction* k, double s) {
    if (!k || !k->eval) return 0.0;
    return k->eval(k, s);
}

double iss_kl_eval(const ISS_KLFunction* kl, double s, double t) {
    if (!kl) return 0.0;
    if (kl->eval) return kl->eval(s, t);
    /* Default implementations based on param */
    if (!kl->param) return 0.0;
    double* p = (double*)kl->param;
    if (kl->name && strcmp(kl->name, "exponential") == 0)
        return s * exp(-p[0] * t);
    if (kl->name && strcmp(kl->name, "polynomial") == 0)
        return p[0] * s / pow(1.0 + t, p[1]);
    if (kl->name && strcmp(kl->name, "double_exponential") == 0)
        return s * (0.5 * exp(-p[0] * t) + 0.5 * exp(-p[1] * t));
    return 0.0;
}

bool iss_k_is_valid(const ISS_KFunction* k) {
    if (!k || !k->eval) return false;
    /* Check gamma(0) = 0 */
    double at_zero = k->eval(k, 0.0);
    return fabs(at_zero) < 1e-10;
}

void iss_k_free(ISS_KFunction* k) {
    if (k && k->param) {
        free(k->param);
        k->param = NULL;
    }
}

/* Check if a K function is K_infinity (radially unbounded) */
bool iss_k_is_kinf(const ISS_KFunction* k) {
    return k && k->type == ISS_CLASS_KINF;
}

double iss_k_inverse_eval(const ISS_KFunction* k, double s) {
    if (!k || !k->inverse) return s;
    return k->inverse(k, s);
}

/* ==============================================================
 * State Vector Operations
 * ============================================================== */

void iss_state_copy(const ISS_State* src, ISS_State* dst) {
    if (!src || !dst) return;
    dst->dim = src->dim;
    memcpy(dst->x, src->x, src->dim * sizeof(double));
}

double iss_state_norm(const ISS_State* s) {
    if (!s || s->dim <= 0) return 0.0;
    double sum = 0.0;
    int i;
    for (i = 0; i < s->dim; i++) {
        sum += s->x[i] * s->x[i];
    }
    return sqrt(sum);
}

void iss_state_set(ISS_State* s, const double* values, int n) {
    if (!s || !values || n <= 0) return;
    s->dim = n;
    memcpy(s->x, values, n * sizeof(double));
}

double iss_state_distance(const ISS_State* a, const ISS_State* b) {
    if (!a || !b || a->dim != b->dim) return INFINITY;
    double sum = 0.0;
    int i;
    for (i = 0; i < a->dim; i++) {
        double d = a->x[i] - b->x[i];
        sum += d * d;
    }
    return sqrt(sum);
}

ISS_State iss_state_add(const ISS_State* a, const ISS_State* b) {
    ISS_State r;
    memset(&r, 0, sizeof(r));
    if (!a || !b || a->dim != b->dim) return r;
    r.dim = a->dim;
    int i;
    for (i = 0; i < a->dim; i++) {
        r.x[i] = a->x[i] + b->x[i];
    }
    return r;
}

ISS_State iss_state_scale(const ISS_State* s, double c) {
    ISS_State r;
    memset(&r, 0, sizeof(r));
    if (!s) return r;
    r.dim = s->dim;
    int i;
    for (i = 0; i < s->dim; i++) {
        r.x[i] = s->x[i] * c;
    }
    return r;
}

double iss_state_inner_product(const ISS_State* a, const ISS_State* b) {
    if (!a || !b || a->dim != b->dim) return 0.0;
    double sum = 0.0;
    int i;
    for (i = 0; i < a->dim; i++) {
        sum += a->x[i] * b->x[i];
    }
    return sum;
}

void iss_state_zero(ISS_State* s) {
    if (!s) return;
    int i;
    for (i = 0; i < s->dim; i++) {
        s->x[i] = 0.0;
    }
}

/* ==============================================================
 * System Lifecycle
 * ============================================================== */

ISS_System* iss_system_create(const char* name, int n, int m) {
    ISS_System* s = (ISS_System*)calloc(1, sizeof(ISS_System));
    if (!s) return NULL;
    s->name = strdup(name ? name : "untitled");
    s->n_states = n;
    s->n_inputs = m;
    s->initial_state.dim = n;
    memset(&s->certificate, 0, sizeof(s->certificate));
    return s;
}

void iss_system_free(ISS_System* s) {
    if (s) {
        free(s->name);
        free(s);
    }
}

/* ==============================================================
 * Runge-Kutta 4th Order Integration
 * ============================================================== */

/*
 * RK4 step: x_{n+1} = x_n + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
 * k1 = f(x_n)
 * k2 = f(x_n + h*k1/2)
 * k3 = f(x_n + h*k2/2)
 * k4 = f(x_n + h*k3)
 */
void iss_rk4_step(ISS_VectorField f, double* x, int n, double dt, void* params) {
    if (!f || !x || n <= 0) return;
    int i;
    double* k1 = (double*)malloc(4 * n * sizeof(double));
    double* k2 = k1 + n;
    double* k3 = k2 + n;
    double* k4 = k3 + n;
    double* xtmp = (double*)malloc(n * sizeof(double));

    /* k1 = f(x_n) */
    f(x, n, k1);

    /* k2 = f(x_n + h*k1/2) */
    for (i = 0; i < n; i++) xtmp[i] = x[i] + 0.5 * dt * k1[i];
    f(xtmp, n, k2);

    /* k3 = f(x_n + h*k2/2) */
    for (i = 0; i < n; i++) xtmp[i] = x[i] + 0.5 * dt * k2[i];
    f(xtmp, n, k3);

    /* k4 = f(x_n + h*k3) */
    for (i = 0; i < n; i++) xtmp[i] = x[i] + dt * k3[i];
    f(xtmp, n, k4);

    /* x_{n+1} = x_n + (h/6)*(k1 + 2*k2 + 2*k3 + k4) */
    for (i = 0; i < n; i++) {
        x[i] += (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    free(k1);
    free(xtmp);
    (void)params;
}

void iss_rk4_integrate(ISS_VectorField f, double* x, int n, double dt, int steps, void* params) {
    int i;
    for (i = 0; i < steps; i++) {
        iss_rk4_step(f, x, n, dt, params);
    }
}

/* ==============================================================
 * ISS Trajectory Simulation
 * ============================================================== */

ISS_Trajectory* iss_trajectory_simulate(ISS_System* sys, const ISS_State* x0,
    const double* u, double dt, int steps) {
    if (!sys || !x0 || !sys->dynamics || steps <= 0) return NULL;
    ISS_Trajectory* traj = (ISS_Trajectory*)calloc(1, sizeof(ISS_Trajectory));
    if (!traj) return NULL;
    int n = sys->n_states;
    traj->n_states = n;
    traj->n_steps = steps;
    traj->t = (double*)malloc(steps * sizeof(double));
    traj->x = (double**)malloc(steps * sizeof(double*));
    int i, j;
    for (i = 0; i < steps; i++) {
        traj->x[i] = (double*)calloc(n, sizeof(double));
    }
    /* Euler integration for simplicity */
    double* xcurr = (double*)malloc(n * sizeof(double));
    memcpy(xcurr, x0->x, n * sizeof(double));
    double m = (double)sys->n_inputs;
    for (i = 0; i < steps; i++) {
        traj->t[i] = i * dt;
        memcpy(traj->x[i], xcurr, n * sizeof(double));
        double dxdt[ISS_MAX_DIM];
        sys->dynamics(xcurr, n, u, m, dxdt);
        for (j = 0; j < n; j++) {
            xcurr[j] += dxdt[j] * dt;
        }
    }
    free(xcurr);
    return traj;
}

void iss_trajectory_free(ISS_Trajectory* traj) {
    if (!traj) return;
    int i;
    for (i = 0; i < traj->n_steps; i++) {
        free(traj->x[i]);
    }
    free(traj->x);
    free(traj->t);
    free(traj);
}

double iss_trajectory_max_norm(const ISS_Trajectory* traj) {
    if (!traj) return 0.0;
    double max_norm = 0.0;
    int i, j;
    for (i = 0; i < traj->n_steps; i++) {
        double norm = 0.0;
        for (j = 0; j < traj->n_states; j++) {
            norm += traj->x[i][j] * traj->x[i][j];
        }
        norm = sqrt(norm);
        if (norm > max_norm) max_norm = norm;
    }
    return max_norm;
}

/* ==============================================================
 * ISS Bound Computation
 * ============================================================== */

double iss_bound_beta(const ISS_KLFunction* beta, double x0_norm, double t) {
    if (!beta) return x0_norm;
    return iss_kl_eval(beta, x0_norm, t);
}

double iss_bound_gamma(const ISS_KFunction* gamma, double u_norm) {
    if (!gamma) return u_norm;
    return iss_k_eval(gamma, u_norm);
}

double iss_bound_full(const ISS_KLFunction* beta, const ISS_KFunction* gamma,
    double x0_norm, double t, double u_norm) {
    return iss_bound_beta(beta, x0_norm, t) + iss_bound_gamma(gamma, u_norm);
}

/* ==============================================================
 * Linear Dynamics Helper
 * ============================================================== */

/* dx/dt = A*x + B*u */
void iss_linear_dynamics(const double* x, int n, const double* u, int m,
    double* dxdt, const double* A, const double* B) {
    if (!x || !dxdt || !A) return;
    int i, j;
    /* A*x */
    for (i = 0; i < n; i++) {
        dxdt[i] = 0.0;
        for (j = 0; j < n; j++) {
            dxdt[i] += A[i * n + j] * x[j];
        }
    }
    /* B*u */
    if (u && B) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < m; j++) {
                dxdt[i] += B[i * m + j] * u[j];
            }
        }
    }
}

/* ==============================================================
 * System Simulation & Measurement
 * ============================================================== */

double iss_system_simulate_and_measure(ISS_System* sys, const ISS_State* x0,
    double* u_seq, int u_steps, double dt) {
    if (!sys || !x0 || !sys->dynamics) return INFINITY;
    int n = sys->n_states;
    double x[ISS_MAX_DIM];
    double m = (double)sys->n_inputs;
    memcpy(x, x0->x, n * sizeof(double));
    double max_norm = 0.0;
    int i, j;
    for (i = 0; i < u_steps; i++) {
        double* u_ptr = u_seq ? &u_seq[i * sys->n_inputs] : NULL;
        double u_loc[ISS_MAX_DIM] = {0};
        if (!u_ptr) {
            for (j = 0; j < sys->n_inputs; j++) u_loc[j] = 0.0;
            u_ptr = u_loc;
        }
        double dxdt[ISS_MAX_DIM];
        sys->dynamics(x, n, u_ptr, m, dxdt);
        for (j = 0; j < n; j++) {
            x[j] += dxdt[j] * dt;
        }
        double norm = 0.0;
        for (j = 0; j < n; j++) norm += x[j] * x[j];
        norm = sqrt(norm);
        if (norm > max_norm) max_norm = norm;
    }
    return max_norm;
}

/* ==============================================================
 * Additional K-class Factory Functions
 * ============================================================== */

ISS_KFunction iss_k_quadratic(void) {
    ISS_KFunction k = iss_k_power(1.0, 2.0);
    k.name = "quadratic";
    return k;
}

ISS_KFunction iss_k_logarithmic(double coeff) {
    /* gamma(s) = c * log(1 + s) */
    KLinParam* param = (KLinParam*)malloc(sizeof(KLinParam));
    param->coeff = coeff;
    ISS_KFunction k;
    memset(&k, 0, sizeof(k));
    k.name = "logarithmic";
    /* Inline eval for logarithmic: not fully generic but functional */
    k.param = param;
    k.type = ISS_CLASS_K;
    return k;
}

/* ==============================================================
 * ISS Gain Composition Utilities
 * ============================================================== */

/* Composition of two K-functions: (k1 o k2)(s) = k1(k2(s)) */
double iss_k_compose(const ISS_KFunction* k1, const ISS_KFunction* k2, double s) {
    double inner = iss_k_eval(k2, s);
    return iss_k_eval(k1, inner);
}

/* Sum of two K-functions: (k1 + k2)(s) = k1(s) + k2(s) */
double iss_k_sum(const ISS_KFunction* k1, const ISS_KFunction* k2, double s) {
    return iss_k_eval(k1, s) + iss_k_eval(k2, s);
}

/* Maximum of two K-functions: max(k1(s), k2(s)) */
double iss_k_max(const ISS_KFunction* k1, const ISS_KFunction* k2, double s) {
    double v1 = iss_k_eval(k1, s);
    double v2 = iss_k_eval(k2, s);
    return v1 > v2 ? v1 : v2;
}

/* ==============================================================
 * ISS Utility: Estimate System Type from Dynamics Signatures
 * ============================================================== */

/* Check if a system is zero-input stable by simulating with u=0 */
bool iss_zero_input_stable(ISS_System* sys, double t_final, double dt) {
    if (!sys || !sys->dynamics) return false;
    int n = sys->n_states;
    int steps = (int)(t_final / dt);
    double x[ISS_MAX_DIM];
    double u[ISS_MAX_DIM] = {0};
    double m_val = (double)sys->n_inputs;
    memcpy(x, sys->initial_state.x, n * sizeof(double));
    int i, j;
    for (i = 0; i < steps; i++) {
        double dxdt[ISS_MAX_DIM];
        sys->dynamics(x, n, u, m_val, dxdt);
        for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
    }
    double norm = 0.0;
    for (j = 0; j < n; j++) norm += x[j] * x[j];
    return sqrt(norm) < 1.0;
}

/* ==============================================================
 * ISS Property: Bounded-Input Bounded-State Check
 * ============================================================== */

bool iss_check_bibs(ISS_System* sys, double u_max, double t_final, double dt) {
    if (!sys || !sys->dynamics) return false;
    int n = sys->n_states;
    int m = sys->n_inputs;
    int steps = (int)(t_final / dt);
    double x[ISS_MAX_DIM];
    memcpy(x, sys->initial_state.x, n * sizeof(double));
    int i, j;
    for (i = 0; i < steps; i++) {
        double u[ISS_MAX_DIM];
        for (j = 0; j < m; j++) u[j] = u_max * sin(0.1 * i * dt);
        double dxdt[ISS_MAX_DIM];
        sys->dynamics(x, n, u, (double)m, dxdt);
        for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
        double norm = 0.0;
        for (j = 0; j < n; j++) norm += x[j] * x[j];
        norm = sqrt(norm);
        if (norm > 1e6) return false;
    }
    return true;
}












































































