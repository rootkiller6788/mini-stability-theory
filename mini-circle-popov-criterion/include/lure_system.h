#ifndef LURE_SYSTEM_H
#define LURE_SYSTEM_H

/*
 * lure_system.h — Lur'e System Construction and Analysis
 *
 * A Lur'e system is the fundamental structure of absolute stability theory.
 * It consists of a linear time-invariant (LTI) system G(s) in the forward
 * path, with a memoryless nonlinearity phi(sigma) in the feedback path.
 *
 * Structure:
 *   Linear part: dx/dt = A*x + B*u,  y = C*x + D*u
 *   Feedback:    u = -phi(y)  (negative feedback)
 *
 *   where phi(.) belongs to a sector [alpha, beta]:
 *     alpha * sigma^2  <=  sigma * phi(sigma)  <=  beta * sigma^2
 *
 * The problem of absolute stability: determine if the origin is globally
 * asymptotically stable for ALL nonlinearities in the given sector.
 *
 * This module provides:
 *   1. Lur'e system construction and lifecycle.
 *   2. Unified stability test (circle and Popov criteria).
 *   3. Loop transformation between equivalent sectors.
 *   4. Gain and phase margin estimation via sector bounding.
 *   5. Multiplier-based tests for specific nonlinearity classes.
 *
 * References:
 *   - Lur'e & Postnikov, "On the theory of stability", 1944
 *   - Aizerman & Gantmacher, "Absolute Stability of Regulator Systems", 1964
 *   - Khalil, Nonlinear Systems, Chapter 7
 */

#include "kyp_lemma.h"

/* ---- Lur'e System ---- */
typedef struct {
    StateSpace  *linear;       /* linear part G(s)                 */
    SectorBounds nonlinear;   /* sector of phi                     */
    double       eta;         /* Popov multiplier (for Popov test) */
    int          config;      /* configuration flags               */
} LureSystem;

/* Configuration flags for LureSystem.config */
#define LURE_CONFIG_STANDARD     0   /* u = -phi(y), negative fb   */
#define LURE_CONFIG_POSITIVE_FB  1   /* u = +phi(y), positive fb   */
#define LURE_CONFIG_MULTIPLIER   2   /* use multiplier in forward   */

/* ---- Stability Test Method ---- */
#define LURE_METHOD_CIRCLE   0   /* Circle criterion               */
#define LURE_METHOD_POPOV    1   /* Popov criterion                */
#define LURE_METHOD_KYP      2   /* KYP lemma LMI approach         */
#define LURE_METHOD_ALL      3   /* Try all methods, return best   */

/* ---- Stability Margin ---- */
typedef struct {
    double gain_margin;      /* gain margin (linearized)           */
    double phase_margin;     /* phase margin in degrees            */
    double sector_margin;    /* margin in sector parameter space   */
    double frequency;        /* critical frequency                 */
} StabilityMargin;

/* ---- Lur'e System Lifecycle ---- */

/* Create a Lur'e system from a state-space linear part and sector bounds.
 * The state-space object is NOT copied — caller retains ownership.
 * For a deep copy, use lure_deep_copy(). */
LureSystem* lure_create(StateSpace* ss, SectorBounds sb);

/* Create a deep copy of a Lur'e system.
 * The linear part's state-space is also deep-copied. */
LureSystem* lure_deep_copy(const LureSystem* ls);

/* Free a Lur'e system. Does NOT free the linear state-space
 * unless the deep_copy flag is set. */
void lure_free(LureSystem* ls);

/* ---- Absolute Stability Tests ---- */

/* Run an absolute stability test using the specified method.
 *
 * Parameters:
 *   ls     - Lur'e system to test
 *   method - LURE_METHOD_CIRCLE, _POPOV, _KYP, or _ALL
 *   margin - output stability margin (>0 means stable)
 *   eta    - Popov multiplier (only for POPOV method; can be NULL)
 *
 * Returns: 1 = absolutely stable, 0 = inconclusive, -1 = error.
 *
 * For LURE_METHOD_ALL, tries circle first (fastest), then Popov,
 * then KYP. Returns the best result. */
int lure_absolute_stability(const LureSystem* ls, int method,
                             double *margin, double *eta);

/* Circle criterion test on a Lur'e system.
 * Constructs G(s), generates Nyquist plot, checks critical disk. */
int lure_circle_test(const LureSystem* ls, double *margin);

/* Popov criterion test on a Lur'e system.
 * Constructs G(s), generates Popov plot, searches for optimal eta. */
int lure_popov_test(const LureSystem* ls, double *margin, double *eta_opt);

/* KYP-based stability test.
 * Converts the Lur'e system to state-space, solves KYP LMI. */
int lure_kyp_test(const LureSystem* ls, double *margin);

/* ---- Margin Computation ---- */

/* Compute gain and phase margins for a Lur'e system.
 * The margins are computed via the linearized system (phi = gain * sigma)
 * and provide a baseline for comparison. */
StabilityMargin* lure_gain_margin(const LureSystem* ls);

/* Compute the robust stability margin by finding the largest sector
 * for which absolute stability can still be certified. */
StabilityMargin* lure_robust_margin(const LureSystem* ls, int method);

/* Free a stability margin structure. */
void lure_margin_free(StabilityMargin* sm);

/* ---- Sector Utilities ---- */

/* Get the lower sector bound (alpha) from a Lur'e system. */
double lure_sector_lower_bound(const LureSystem* ls);

/* Get the upper sector bound (beta) from a Lur'e system. */
double lure_sector_upper_bound(const LureSystem* ls);

/* Check if a given nonlinear function value satisfies the sector condition.
 * Returns true if alpha*sigma^2 <= sigma*phi <= beta*sigma^2. */
bool lure_sector_check(const LureSystem* ls, double sigma, double phi);

/* ---- Loop Transformation ---- */

/* Apply loop transformation to the Lur'e system.
 * Transforms the sector from [alpha, beta] to [0, k] where k = beta - alpha.
 * The transformed linear part is: G'(s) = G(s) / (1 + alpha*G(s)).
 *
 * Returns a new LureSystem with the transformed linear part.
 * Caller must free with lure_free(). */
LureSystem* lure_loop_transform(const LureSystem* ls);

/* ---- Multiplier Methods ---- */

/* Apply the Zames-Falb multiplier to the Lur'e system.
 * This generalizes both circle and Popov criteria using a class of
 * convolution multipliers. Returns the augmented Lur'e system. */
LureSystem* lure_zames_falb_multiplier(const LureSystem* ls,
                                         double *multiplier_gain);

/* ---- System Properties ---- */

/* Compute the H-infinity norm of the linear part.
 * Useful for small-gain theorem comparisons. */
double lure_hinf_norm(const LureSystem* ls);

/* Check if the linear part is stable (all eigenvalues of A in LHP). */
bool lure_linear_is_stable(const LureSystem* ls);

/* ---- Constants ---- */
#define LURE_MAX_SECTOR   1000.0  /* maximum sector bound considered */
#define LURE_NUM_FREQS     500    /* frequency grid for tests       */

#endif /* LURE_SYSTEM_H */
