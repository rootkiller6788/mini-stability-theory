#ifndef INTERCON_CORE_H
#define INTERCON_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* ==============================================================
 * Stability of Interconnected Systems - Core Types & Operations
 *
 * Based on: Zames (1966), Willems (1972), Moylan & Hill (1978),
 *           Siljak (1978), Vidyasagar (1980)
 *
 * intercon_core.h: subsystem types, interconnection topology,
 *                  composite system construction, core operations.
 *
 * Companion source: intercon_core.c
 *
 * Dependent headers:
 *   intercon_small_gain.h    - Small-gain theorem
 *   intercon_passivity.h     - Passivity-based stability
 *   intercon_lyapunov.h      - Lyapunov methods for composites
 *   intercon_decentralized.h - Decentralized stability
 * ============================================================== */

#define IC_PI          3.14159265358979323846
#define IC_EPSILON     1e-12
#define IC_MAX_STATES  16
#define IC_MAX_SUBSYS  8

/* Subsystem dimensions */
#define IC_MAX_DIM     8

/* ==============================================================
 * Interconnection Topology Types
 * ============================================================== */
typedef enum {
    IC_SERIES        = 0,  /* S1 -> S2 (cascade) */
    IC_PARALLEL      = 1,  /* S1 || S2 (additive) */
    IC_FEEDBACK_POS  = 2,  /* positive feedback loop */
    IC_FEEDBACK_NEG  = 3,  /* negative feedback loop */
    IC_GENERAL       = 4   /* arbitrary interconnection matrix */
} ICTopology;

/* ==============================================================
 * Subsystem State-Space Model
 *
 * dx/dt = A*x + B*u
 * y     = C*x + D*u
 *
 * For discrete-time: x[k+1] = A*x[k] + B*u[k]
 * ============================================================== */
typedef struct {
    int     n;                  /* state dimension */
    int     m;                  /* input dimension */
    int     p;                  /* output dimension */
    double  A[IC_MAX_DIM][IC_MAX_DIM];
    double  B[IC_MAX_DIM][IC_MAX_DIM];
    double  C[IC_MAX_DIM][IC_MAX_DIM];
    double  D[IC_MAX_DIM][IC_MAX_DIM];
    double  x[IC_MAX_DIM];      /* current state */
    bool    is_discrete;
} ICSubsystem;

/* ==============================================================
 * Interconnection Matrix
 *
 * For N subsystems, the interconnection is specified by
 * an NxN block matrix H where H_ij describes how the
 * output of subsystem j feeds into the input of subsystem i.
 *
 * u_i = sum_j H_ij * y_j + r_i  (r_i = external input)
 * ============================================================== */
typedef struct {
    int     n_subsys;
    double  H[IC_MAX_SUBSYS][IC_MAX_SUBSYS];
    bool    active[IC_MAX_SUBSYS];
} ICInterconnection;

/* ==============================================================
 * Interconnected System
 *
 * Composed of multiple ICSubsystem blocks connected via an
 * ICInterconnection matrix. The composite system state is
 * the stacked state vectors of all subsystems.
 * ============================================================== */
typedef struct {
    int             n_subsys;
    ICSubsystem*    subsys[IC_MAX_SUBSYS];
    ICInterconnection* intercon;
    double          composite_x[IC_MAX_STATES];
    int             total_states;
    double          t;
} ICSystem;

/* ==============================================================
 * Stability Result
 * ============================================================== */
typedef enum {
    IC_UNKNOWN       = 0,
    IC_ASYM_STABLE   = 1,   /* asymptotically stable */
    IC_EXP_STABLE    = 2,   /* exponentially stable */
    IC_MARG_STABLE   = 3,   /* marginally stable */
    IC_UNSTABLE      = 4,   /* unstable */
    IC_BIBO_STABLE   = 5    /* bounded-input bounded-output stable */
} ICStability;

/* ==============================================================
 * Subsystem Lifecycle
 * ============================================================== */

ICSubsystem* ic_subsys_create(int n, int m, int p, bool discrete);
void         ic_subsys_free(ICSubsystem* s);
void         ic_subsys_set_A(ICSubsystem* s, const double* A_flat);
void         ic_subsys_set_B(ICSubsystem* s, const double* B_flat);
void         ic_subsys_set_C(ICSubsystem* s, const double* C_flat);
void         ic_subsys_set_D(ICSubsystem* s, const double* D_flat);
void         ic_subsys_set_state(ICSubsystem* s, const double* x0);
void         ic_subsys_get_state(const ICSubsystem* s, double* x);
void         ic_subsys_step(ICSubsystem* s, const double* u, double dt);
void         ic_subsys_output(const ICSubsystem* s, const double* u, double* y);
void         ic_subsys_print(const ICSubsystem* s);

/* Subsystem transfer function evaluation */
void         ic_subsys_freq_response(const ICSubsystem* s, double omega,
                                      double* mag, double* phase);

/* ==============================================================
 * Interconnection Matrix Operations
 * ============================================================== */

ICInterconnection* ic_intercon_create(int n_subsys);
void               ic_intercon_free(ICInterconnection* ic);
void               ic_intercon_set(ICInterconnection* ic, int from, int to,
                                    double weight);
double             ic_intercon_get(const ICInterconnection* ic, int from, int to);
void               ic_intercon_spectral_radius(const ICInterconnection* ic,
                                                double* radius);
bool               ic_intercon_is_connected(const ICInterconnection* ic);
void               ic_intercon_print(const ICInterconnection* ic);

/* ==============================================================
 * Interconnected System Lifecycle
 * ============================================================== */

ICSystem* ic_system_create(void);
void      ic_system_free(ICSystem* sys);
int       ic_system_add_subsys(ICSystem* sys, ICSubsystem* sub);
bool      ic_system_connect(ICSystem* sys, ICTopology topo, int idx1, int idx2,
                              double gain);
void      ic_system_step(ICSystem* sys, const double* r, double dt);
void      ic_system_get_output(const ICSystem* sys, const double* r, double* y);
void      ic_system_get_composite_state(const ICSystem* sys, double* x);
void      ic_system_print(const ICSystem* sys);

/* ==============================================================
 * Composite Matrix Operations
 * ============================================================== */

/* Build composite A matrix for interconnected system */
void ic_build_composite_A(const ICSystem* sys, double* A_comp, int max_dim);

/* Check eigenvalues of composite A for stability */
ICStability ic_composite_eigenvalue_stability(const ICSystem* sys);

/* ==============================================================
 * Utility Functions
 * ============================================================== */

void  ic_matrix_mult(const double* A, const double* B, double* C,
                      int m, int n, int p);
void  ic_matrix_add(const double* A, const double* B, double* C,
                     int m, int n);
void  ic_matrix_scale(const double* A, double alpha, double* B,
                       int m, int n);
double ic_matrix_norm_inf(const double* A, int m, int n);
double ic_matrix_norm_2(const double* A, int m, int n);
int   ic_eigenvalues_real(double* A, int n, double* real_parts,
                            double* imag_parts, int max_eig);
void  ic_solve_lyapunov(double* A, double* Q, double* P, int n);

#endif /* INTERCON_CORE_H */
double ic_eigenvalue_bound_gershgorin(const double* A, int n);
bool   ic_is_hurwitz(const double* A, int n);
void   ic_matrix_exp_pade(const double* A, int n, double t, double* expAt);
bool   ic_is_observable_pair(const double* A, const double* C, int n, int p);
bool   ic_is_controllable_pair(const double* A, const double* B, int n, int m);

double ic_lyapunov_condition_number(const double* A, int n);
double ic_distance_to_instability(const double* A, int n);
double ic_hankel_norm_bound(const double* A, const double* B, const double* C, int n);
double ic_h2_norm_approx(const double* A, const double* B, const double* C, int n);
void   ic_system_reset(ICSystem* sys);

void   ic_cross_gramian(const double* A, const double* B, const double* C, int n, double* X);
double ic_balanced_realization(const double* A, int n, double* T);
double ic_impulse_energy(const ICSubsystem* s, double T_end, double dt);
bool   ic_is_minimum_phase(const ICSubsystem* s);
double ic_stability_radius(const double* A, int n);

/* Extended Kalman-Yakubovich-Popov lemma variant */
double ic_kyp_lemma_check(const double* A, const double* B, const double* C, int n);
