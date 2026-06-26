#ifndef LDM_CORE_H
#define LDM_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * ldm_core.h - Lyapunov Direct Method Core Types
 *
 * Lyapunov (1892) "The General Problem of Stability of Motion"
 * The direct method determines stability without solving ODEs.
 *
 * Key theorem: If there exists V(x) such that:
 *   V(0)=0, V(x)>0 for x!=0 (positive definite)
 *   dV/dt = grad(V)*f(x) <= 0 (negative semi-definite)
 * Then origin is stable.
 * If dV/dt < 0, origin is asymptotically stable.
 * If V is radially unbounded, stability is global.
 *
 * References:
 *   Lyapunov (1892) Dissertation, Kharkov
 *   Khalil (2002) Nonlinear Systems, 3rd ed.
 *   Sastry (1999) Nonlinear Systems: Analysis, Stability, Control
 * ============================================================== */

typedef enum { LDM_STABLE=0, LDM_ASYM_STABLE=1, LDM_EXP_STABLE=2,
               LDM_GLOBAL_ASYM=3, LDM_GLOBAL_EXP=4, LDM_UNSTABLE=5,
               LDM_INCONCLUSIVE=6 } LDM_StabilityResult;

typedef enum { LDM_QUADRATIC=0, LDM_POLYNOMIAL=2, LDM_PIECEWISE=3,
               LDM_RADIAL=4, LDM_LOGARITHMIC=5, LDM_SUM_OF_SQUARES=6 } LDM_FunctionType;

typedef struct {
    double value;
    double derivative;
    double* grad;
    int n_dims;
    LDM_FunctionType type;
} LDM_LyapunovEval;

typedef struct {
    double* P;
    int n;
    bool is_pd;
    double min_eigenvalue;
    double max_eigenvalue;
    double condition_number;
} LDM_QuadraticForm;

typedef struct {
    double* A; int n;
    double* x;
    double t;
    bool is_linear;
    void (*dynamics)(double t, const double* x, int n, double* dx);
    void* context;
} LDM_System;

typedef struct {
    LDM_StabilityResult result;
    LDM_QuadraticForm* P;
    double convergence_rate;
    double* trajectory;
    int traj_len;
    double final_V;
    double final_dV;
} LDM_StabilityReport;

/* Core API — system lifecycle */
LDM_System* ldm_system_create(int n, void (*dynamics)(double,const double*,int,double*), void* ctx);
LDM_System* ldm_system_create_linear(int n, const double* A);
void ldm_system_set_state(LDM_System* sys, const double* x0);
void ldm_system_free(LDM_System* sys);
void ldm_system_step(LDM_System* sys, double dt);
void ldm_system_step_euler(LDM_System* sys, double dt);

/* Lyapunov evaluation */
LDM_LyapunovEval ldm_eval_quadratic(const LDM_System* sys, const double* P, int n);
LDM_LyapunovEval ldm_eval_general(const LDM_System* sys, double (*V)(const double*,int), double (*dV)(const double*,int,const double*,int), int n);

/* Quadratic forms */
double ldm_quadratic_form(const double* P, const double* x, int n);
double ldm_quadratic_form_sym(const double* P, const double* x, int n);

/* Definiteness tests (any n, via Cholesky) */
bool ldm_is_positive_definite(const double* P, int n);
bool ldm_is_positive_semidefinite(const double* P, int n);
bool ldm_is_negative_definite(const double* P, int n);
bool ldm_is_negative_semidefinite(const double* P, int n);

/* Quadratic form lifecycle */
LDM_QuadraticForm* ldm_quadratic_create(int n);
void ldm_quadratic_free(LDM_QuadraticForm* qf);

/* Eigenvalue bounds */
void ldm_quadratic_eigen_bounds(const double* P, int n, double* lmin, double* lmax);
void ldm_gershgorin_bounds(const double* P, int n, double* lmin, double* lmax);

/* Vector operations */
double ldm_vector_norm(const double* x, int n);
double ldm_vector_norm_inf(const double* x, int n);
double ldm_vector_dot(const double* a, const double* b, int n);

/* Matrix utilities */
double ldm_matrix_norm_frobenius(const double* A, int n);
void ldm_matrix_add(const double* A, const double* B, double* C, int n);
void ldm_matrix_scale(double* A, double s, int n);
double ldm_lyapunov_residual(const double* A, const double* P, const double* Q, int n);
int ldm_matrix_invert(double* A, int n);
double ldm_matrix_determinant(const double* A, int n);
double ldm_matrix_condition(const double* A, int n);

#endif
