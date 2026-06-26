#ifndef INTERCON_DECENTRALIZED_H
#define INTERCON_DECENTRALIZED_H
#include "intercon_core.h"

/* ==============================================================
 * Decentralized Stability (Siljak 1978, Sandell 1978)
 *
 * Decentralized control: each subsystem has a local controller
 * with limited information. Stability must be guaranteed despite
 * interconnections treated as perturbations.
 * ============================================================== */

typedef struct {
    double K[IC_MAX_DIM][IC_MAX_DIM]; /* local gain matrix */
    int    n;                          /* state dimension */
    int    m;                          /* input dimension */
    bool   is_stabilizing;
} DecentralizedController;

typedef struct {
    int    n_subsys;
    DecentralizedController* controllers[IC_MAX_SUBSYS];
    double interaction_bound;       /* max norm of interconnection */
    double stability_margin;
    bool   is_stable;
} DecentralizedStability;

/* ==============================================================
 * Decentralized Controller Design
 * ============================================================== */

DecentralizedController* ic_decentralized_controller_create(int n, int m);
void  ic_decentralized_controller_free(DecentralizedController* dc);
void  ic_decentralized_controller_set_gain(DecentralizedController* dc,
                                             const double* K_flat);
bool  ic_decentralized_stabilize_check(const DecentralizedController* dc,
                                         const ICSubsystem* s);

/* ==============================================================
 * Decentralized Stability Analysis
 * ============================================================== */

DecentralizedStability* ic_decentralized_stability_test(const ICSystem* sys);
void    ic_decentralized_stability_free(DecentralizedStability* ds);
void    ic_decentralized_stability_print(const DecentralizedStability* ds);

/* ==============================================================
 * Interaction Bound Computation
 *
 * For decentralized stability, we need:
 *   min_i lambda_min(Q_i) > sum_j ||A_ij||
 *
 * where Q_i is from the Lyapunov equation for subsystem i.
 * The RHS is the "interaction bound" that must be smaller than
 * the stability margin of each isolated subsystem.
 * ============================================================== */

double ic_interaction_norm_bound(const ICSystem* sys);

/* ==============================================================
 * Overlapping Decomposition (Siljak 1978)
 *
 * Decompose system with overlapping states into disjoint
 * subsystems via inclusion principle.
 * ============================================================== */

typedef struct {
    double* V;       /* expansion matrix */
    double* U;       /* contraction matrix */
    int     n_orig;  /* original dimension */
    int     n_exp;   /* expanded dimension */
} OverlappingDecomp;

OverlappingDecomp* ic_overlapping_decompose(const double* A, int n);
void               ic_overlapping_decomp_free(OverlappingDecomp* od);
void               ic_overlapping_decomp_print(const OverlappingDecomp* od);

/* ==============================================================
 * Decentralized Fixed Modes
 *
 * A fixed mode is an eigenvalue of the composite system that
 * cannot be moved by any decentralized feedback.
 * ============================================================== */

int  ic_fixed_modes_count(const ICSystem* sys);
bool ic_is_decentralized_controllable(const ICSystem* sys);
bool ic_is_decentralized_observable(const ICSystem* sys);

/* ==============================================================
 * Sparse Decentralized H-infinity Control
 * ============================================================== */

typedef struct {
    double K_diag[IC_MAX_SUBSYS][IC_MAX_DIM][IC_MAX_DIM];
    int    n_subsys;
    double gamma;    /* H-infinity performance level */
    bool   is_feasible;
} SparseHinfController;

SparseHinfController* ic_sparse_hinf_design(const ICSystem* sys,
                                               double gamma_target);
void   ic_sparse_hinf_free(SparseHinfController* sc);
void   ic_sparse_hinf_print(const SparseHinfController* sc);

#endif /* INTERCON_DECENTRALIZED_H */
bool   ic_decentralized_lqr_design(const ICSubsystem* s, double Qw, double Rw, double* K_out);
double ic_relative_gain_array(const ICSystem* sys, int i, int j);
bool   ic_decentralized_integral_controllability(const ICSystem* sys);

bool   ic_fixed_modes_stable(const ICSystem* sys);
double ic_decentralized_performance_loss(const ICSystem* sys);
