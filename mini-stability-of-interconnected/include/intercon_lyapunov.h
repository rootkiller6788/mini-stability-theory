#ifndef INTERCON_LYAPUNOV_H
#define INTERCON_LYAPUNOV_H
#include "intercon_core.h"

/* ==============================================================
 * Lyapunov Methods for Interconnected Systems
 *
 * Composite Lyapunov functions: V(x) = sum_i d_i * V_i(x_i)
 * Vector Lyapunov functions:   V(x) = [V_1(x_1), ..., V_N(x_N)]^T
 * M-matrix stability conditions
 * Gershgorin diagonal dominance
 * ============================================================== */

typedef struct {
    double*   V;         /* vector Lyapunov function values */
    int       n;         /* number of subsystems */
    double*   d;         /* weighting coefficients */
    double    V_composite; /* scalar composite value */
    double    decay_rate;
    bool      is_stable;
} CompositeLyapunov;

typedef struct {
    double    M[IC_MAX_SUBSYS][IC_MAX_SUBSYS]; /* interconnection matrix */
    int       n;
    bool      is_m_matrix;
    double    stability_margin;
} MMatrixResult;

/* ==============================================================
 * Composite Lyapunov Function
 * ============================================================== */

/* Construct composite Lyapunov V(x) = sum d_i * x_i^T * P_i * x_i */
CompositeLyapunov* ic_composite_lyapunov(const ICSystem* sys,
                                           const double* weights);
void               ic_composite_lyapunov_free(CompositeLyapunov* cl);
double             ic_composite_lyapunov_eval(const CompositeLyapunov* cl,
                                                const double* x);
double             ic_composite_lyapunov_derivative(const CompositeLyapunov* cl,
                                                      const ICSystem* sys);
void               ic_composite_lyapunov_print(const CompositeLyapunov* cl);

/* ==============================================================
 * Vector Lyapunov Function (Bellman 1962, Siljak 1978)
 *
 * V_dot <= M * V (componentwise)
 * If M is an M-matrix, the origin is globally asymptotically stable.
 * ============================================================== */

/* Build comparison system M from subsystem Lyapunov functions */
MMatrixResult* ic_build_M_matrix(const ICSystem* sys);
void           ic_M_matrix_free(MMatrixResult* mr);

/* Check if matrix is an M-matrix (all off-diagonals <= 0, all principal minors > 0) */
bool           ic_is_M_matrix(const double* M_flat, int n);

/* Check if M-matrix condition implies composite stability */
bool           ic_M_matrix_stability(const MMatrixResult* mr);
void           ic_M_matrix_print(const MMatrixResult* mr);

/* ==============================================================
 * Gershgorin Circle Theorem for Stability
 *
 * All eigenvalues of A lie within the union of Gershgorin discs:
 *   D_i = { z : |z - a_ii| <= R_i } where R_i = sum_{j!=i} |a_ij|
 *
 * For stability: all discs must lie in left half-plane.
 * For diagonal dominance: |a_ii| > R_i for all i
 * ============================================================== */

typedef struct {
    double centers[IC_MAX_STATES];
    double radii[IC_MAX_STATES];
    int    n;
    bool   is_diag_dominant;
    bool   all_left_half;
    double max_real_part;   /* rightmost possible eigenvalue */
} GershgorinResult;

GershgorinResult* ic_gershgorin_test(const double* A, int n);
void              ic_gershgorin_free(GershgorinResult* gr);
bool              ic_gershgorin_stable(const GershgorinResult* gr);
void              ic_gershgorin_print(const GershgorinResult* gr);

/* ==============================================================
 * Diagonal Dominance & Stability
 * ============================================================== */

/* Check row diagonal dominance */
bool ic_is_row_diag_dominant(const double* A, int n);

/* Check column diagonal dominance */
bool ic_is_col_diag_dominant(const double* A, int n);

/* Generalized diagonal dominance: exists D>0 s.t. D*A*D^{-1} is diag dominant */
bool ic_is_generalized_diag_dominant(const double* A, int n);

/* ==============================================================
 * Comparison Principle (Siljak 1978)
 *
 * If V_dot <= W(t, V) and U_dot = W(t, U), U(0) = V(0),
 * then V(t) <= U(t) for all t >= 0.
 *
 * This allows bounding the composite system by a simpler comparison system.
 * ============================================================== */

void ic_comparison_system_step(const double* M, int n, const double* V_in,
                                 double* V_out, double dt);
void ic_comparison_system_simulate(const double* M, int n, const double* V0,
                                     double duration, double dt, double* V_final);

/* ==============================================================
 * Connective Stability (Siljak 1978)
 *
 * A system is connectively stable if it remains stable under
 * arbitrary disconnection/reconnection of subsystems.
 *
 * This is equivalent to the existence of a common Lyapunov function
 * for all possible interconnection topologies.
 * ============================================================== */

bool ic_connective_stability_test(const ICSystem* sys);

/* ==============================================================
 * Stability Margin Computation
 * ============================================================== */

/* Compute stability margin via structured singular value (mu) approximation */
double ic_stability_margin_mu(const ICSystem* sys);

/* Compute gain margin for feedback interconnection */
double ic_gain_margin(const ICSubsystem* s1, const ICSubsystem* s2);

/* Compute phase margin for feedback interconnection */
double ic_phase_margin(const ICSubsystem* s1, const ICSubsystem* s2);

#endif /* INTERCON_LYAPUNOV_H */
double ic_estimate_roa(const ICSystem* sys);
bool   ic_absolute_stability_sector(const ICSubsystem* s, double k);
double ic_max_interconnection_gain(const ICSystem* sys);

double ic_m_matrix_decay_rate(const MMatrixResult* mr);
bool   ic_circle_time_varying(const ICSubsystem* s, double k1, double k2);
double ic_min_damping_ratio(const ICSystem* sys);
