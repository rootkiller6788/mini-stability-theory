#ifndef FTS_CORE_H
#define FTS_CORE_H

#include <stdbool.h>
#include <stddef.h>

/* ============================================================
 * Finite-Time Stability Core Definitions
 *
 * Numerical simulation framework for FTS analysis.
 * Key refs: Bhat & Bernstein (2000), Polyakov (2012),
 *           Moulay & Perruquetti (2006)
 * ============================================================ */

#define FTS_EPSILON    1e-12
#define FTS_MAX_DIM    10
#define FTS_MAX_HIST   100000
#define FTS_DEFAULT_DT 0.01
#define FTS_INF        1e99

typedef struct { double x[FTS_MAX_DIM]; int dim; } FTSState;
typedef struct { double re, im; } FTSComplex;

typedef struct {
    FTSState* pts; int n, cap;
} FTSTrajectory;

typedef enum {
    FTS_ASYMPTOTIC = 0, FTS_EXPONENTIAL = 1,
    FTS_FINITE_TIME = 2, FTS_FIXED_TIME = 3,
    FTS_PRACTICAL = 4
} FTSConvergenceType;

typedef enum {
    FTS_NORM_L2 = 0, FTS_NORM_L1 = 1,
    FTS_NORM_LINF = 2, FTS_NORM_WEIGHTED = 3
} FTSNormType;

typedef struct { double t, error; int step_idx; } FTSErrorPoint;

typedef struct FTSSystem {
    char* name; int dim; FTSState state; double t, dt;
    double* params; int np;
    void (*rhs)(struct FTSSystem*, FTSState, FTSState*);
    FTSTrajectory traj;
} FTSSystem;

/* Lifecycle */
FTSSystem* fts_create(int dim, int np, double dt);
void       fts_free(FTSSystem* sys);

/* State access */
void   fts_set_state(FTSSystem* sys, double* x);
double fts_get_state(FTSSystem* sys, int idx);
void   fts_copy_state(FTSState* dst, const FTSState* src);
void   fts_set_param(FTSSystem* sys, int idx, double v);
double fts_get_param(FTSSystem* sys, int idx);

/* Integration */
void fts_step_rk4(FTSSystem* sys);
void fts_step_euler(FTSSystem* sys);
void fts_step_heun(FTSSystem* sys);
void fts_integrate(FTSSystem* sys, double dur, int rec);

/* Trajectory */
void fts_record(FTSSystem* sys);
void fts_clear_trajectory(FTSSystem* sys);
int  fts_trajectory_count(FTSSystem* sys);

/* Norm and distance */
double fts_norm(FTSState s);
double fts_norm_type(FTSState s, FTSNormType type);
double fts_distance(FTSState a, FTSState b, int dim);
double fts_squared_norm(FTSState s);
double fts_dot_product(FTSState a, FTSState b);

/* Vector operations */
void fts_vector_add(FTSState* dst, FTSState a, FTSState b, double s);
void fts_vector_scale(FTSState* dst, FTSState src, double s);
void fts_vector_sub(FTSState* dst, FTSState a, FTSState b);
void fts_vector_copy(FTSState* dst, FTSState src);

/* System analysis */
void   fts_print(FTSSystem* sys);
int    fts_dimension(FTSSystem* sys);
double fts_current_time(FTSSystem* sys);
bool   fts_is_valid(FTSSystem* sys);

/* Factory functions */
FTSSystem* fts_create_double_integrator(double dt);
FTSSystem* fts_create_scalar_fts(double alpha, double dt);
FTSSystem* fts_create_harmonic_oscillator(double omega, double dt);
FTSSystem* fts_create_van_der_pol(double mu, double dt);
FTSSystem* fts_create_lorenz(double sigma, double rho, double beta, double dt);
FTSSystem* fts_create_brockett_integrator(double dt);

/* Energy */
double fts_energy_quadratic(FTSState x, double* P, int dim);
double fts_lyap_default(FTSState x, void* params);

/* Finite escape */
bool   fts_has_finite_escape(FTSSystem* sys, double Tmax);
double fts_divergence(FTSSystem* sys, FTSState x);
int    fts_escape_time(FTSSystem* sys, double thresh, double Tmax, double* te);

/* Utility */
void   fts_sleep_ms(int ms);
double fts_clamp(double x, double lo, double hi);
double fts_sign(double x);
double fts_saturate(double x, double limit);

#endif /* FTS_CORE_H */
