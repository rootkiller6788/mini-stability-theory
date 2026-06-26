#ifndef FTS_SLIDING_H
#define FTS_SLIDING_H

#include "fts_core.h"

/* ============================================================
 * Sliding Mode Control for Finite-Time Convergence
 *
 * Key concepts:
 * - Terminal sliding surface: s = x1 + beta*sig(x1)^alpha
 *   where sig(x)^a = |x|^a * sign(x)
 * - Reaching phase: finite-time convergence to s=0
 * - Sliding phase: finite-time convergence to origin
 * - Super-twisting: 2nd-order SMC with continuous control
 *
 * Ref: Utkin et al. (2009), Levant (2003), Feng et al. (2002)
 * ============================================================ */

typedef struct {
    double* surface;    /* Surface coefficients */
    double gain;        /* Sliding gain */
    double exponent;    /* Fractional power exponent */
    int order;          /* Sliding surface order */
    bool is_terminal;   /* True for terminal sliding surface */
    int dim;            /* Dimension of surface */
} FTSSlidingSurface;

typedef struct {
    FTSSlidingSurface sigma;  /* Sliding surface */
    double* control;          /* Control input history */
    double* state_norm;       /* Norm history */
    int dim_ctrl;             /* Control dimension */
    int n_steps;              /* Simulation steps */
    double reaching_time;     /* Time to reach s=0 */
    double sliding_time;      /* Time spent in sliding mode */
    double settling_time;     /* Time to reach equilibrium */
    bool converged;           /* Convergence flag */
    double chattering;        /* Chattering measure */
    double energy;            /* Control energy */
} FTSTerminalSMC;

/* --- Surface Management --- */
FTSSlidingSurface* fts_smc_surface_create(int dim, double gain,
                                           double exponent, int order);
void fts_smc_surface_free(FTSSlidingSurface* s);
double fts_smc_eval_surface(FTSSlidingSurface* s, FTSState x);
double fts_smc_surface_derivative(FTSSlidingSurface* s, FTSState x,
                                   FTSState dx);

/* --- SMC Controller --- */
FTSTerminalSMC* fts_smc_create(FTSSlidingSurface* s, int dim_ctrl);
void fts_smc_free(FTSTerminalSMC* smc);
int  fts_smc_simulate(FTSTerminalSMC* smc, FTSSystem* sys,
                      FTSState x0, double T, double dt);

/* --- Control Law --- */
double fts_smc_control(FTSTerminalSMC* smc, FTSState x);
double fts_smc_super_twisting_control(FTSState x, double lambda, double alpha);
double fts_smc_terminal_control(FTSState x, double k, double alpha, double beta);
double fts_smc_nonsingular_tsm_control(FTSState x, double k,
                                        double alpha, double beta);

/* --- Bounds --- */
double fts_smc_reaching_time_bound(FTSSlidingSurface* s,
                                    FTSState x0, double gain);
double fts_smc_sliding_time_bound(FTSSlidingSurface* s, FTSState x0);
double fts_smc_settling_time_estimate(FTSTerminalSMC* smc);

/* --- Analysis --- */
bool   fts_smc_is_converged(FTSTerminalSMC* smc, double tol);
bool   fts_smc_is_nonsingular(FTSSlidingSurface* s, FTSState x);
double fts_smc_control_effort(FTSTerminalSMC* smc);
double fts_smc_chattering_measure(FTSTerminalSMC* smc);
double fts_smc_energy(FTSTerminalSMC* smc);

/* --- Output --- */
void fts_smc_print(FTSTerminalSMC* smc);
void fts_smc_print_surface(FTSSlidingSurface* s);

/* --- Adaptive SMC --- */
void fts_smc_adapt_gain(FTSSlidingSurface* s, FTSState x, double s_val,
                         double adaptation_rate);

/* --- Observer-Based SMC --- */
double fts_smc_observer_gain(double lambda, double L);
void fts_smc_levants_differentiator(double* x_hat, double* x,
                                     double dt, int order, double L);

#endif /* FTS_SLIDING_H */
