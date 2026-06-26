#ifndef ABS_AIZERMAN_H
#define ABS_AIZERMAN_H

/*
 * abs_aizerman.h — Aizerman Conjecture and Counterexamples
 *
 * Aizerman's conjecture (1949):
 *   If the linearized system dx/dt = (A - k*b*c^T)*x is stable
 *   for all k in [alpha, beta], then the nonlinear Lur'e system
 *   with any phi in sector [alpha, beta] is globally asymptotically stable.
 *
 * This conjecture is FALSE in general (counterexamples by Pliss 1958,
 * Fitts 1966), but holds for:
 *   - First and second order systems
 *   - Systems with A symmetric (by circle criterion)
 *   - Specific structural conditions
 *
 * The Aizerman problem motivates the search for frequency-domain
 * criteria (circle, Popov) that provide sufficient conditions.
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "abs_core.h"

/* ── Aizerman Test Result ─────────────────────────────────────────── */

typedef enum {
    ABS_AIZ_HOLDS_PROVEN    = 0,  /* conjecture holds (provably) */
    ABS_AIZ_HOLDS_EMPIRICAL = 1,  /* conjecture holds numerically */
    ABS_AIZ_FAILS           = 2,  /* conjecture fails (counterexample exists) */
    ABS_AIZ_UNDECIDED       = 3   /* cannot determine */
} AbsAizermanStatus;

typedef struct {
    AbsAizermanStatus status;
    double  k_min;           /* lower Hurwitz bound */
    double  k_max;           /* upper Hurwitz bound */
    double  margin;          /* stability margin of worst-case linearization */
    double  gamma_min;       /* minimum H-infinity gain margin */
    bool    is_symmetric;    /* true if A is symmetric */
    int     order;           /* system order */
    char    reason[256];     /* human-readable explanation */
} AbsAizermanResult;

/* ── Aizerman Conjecture Testing ──────────────────────────────────── */

/** Test Aizerman's conjecture for a given Lur'e system.
 *
 *  1. Compute the Hurwitz sector [k_min, k_max] for A - k*b*c^T.
 *  2. Check if sector [alpha, beta] is within [k_min, k_max].
 *  3. Apply known results (order 1-2 -> holds; A symmetric -> holds).
 *  4. Check known counterexample patterns.
 *  5. Apply frequency-domain criteria as supplementary evidence. */
AbsAizermanResult* abs_aizerman_test(const AbsLureSystem *sys);

/** Free Aizerman result. */
void abs_aizerman_result_free(AbsAizermanResult *res);

/* ── Hurwitz Sector Computation ───────────────────────────────────── */

/** Compute the Hurwitz sector for parameterized matrix A - k*B*C.
 *
 *  For each k, the eigenvalues of A - k*B*C must have negative real parts.
 *  The set of k where this holds is an interval (k_min, k_max) for
 *  rank-1 perturbations B*C = b*c^T.
 *
 *  Uses the characteristic polynomial and Routh-Hurwitz criterion
 *  for small n, or root-locus computation for larger n. */
bool abs_aizerman_hurwitz_sector(const double *A, const double *b,
                                  const double *c, int n,
                                  double *k_min, double *k_max);

/** Check if a specific k makes A - k*b*c^T Hurwitz stable.
 *  Verifies all eigenvalues have negative real parts. */
bool abs_aizerman_is_stable_k(const double *A, const double *b,
                               const double *c, int n, double k);

/** Compute the root locus: eigenvalues of A - k*b*c^T for k in [k1, k2].
 *  Stores real/imag parts in pre-allocated arrays (num_k × n doubles each).
 *  Returns the critical k where an eigenvalue crosses the imaginary axis. */
double abs_aizerman_root_locus(const double *A, const double *b,
                                const double *c, int n,
                                double k1, double k2, int num_k,
                                double *eig_re, double *eig_im);

/* ── Counterexample Construction ──────────────────────────────────── */

/** Build the Pliss counterexample (3rd order, demonstrates conjecture
 *  fails for n >= 3 with appropriate parameters).
 *  Fills A[n*n], b[n], c[n] with the canonical counterexample. */
void abs_aizerman_pliss_counterexample(double *A, double *b, double *c);

/** Build the Fitts counterexample (4th order).
 *  More dramatic failure of the conjecture. */
void abs_aizerman_fitts_counterexample(double *A, double *b, double *c);

/** Check if a given system matches known counterexample patterns.
 *  Returns true if this system is structurally a counterexample. */
bool abs_aizerman_is_counterexample(const double *A, const double *b,
                                     const double *c, int n);

/* ── Aizerman vs Circle/Popov Comparison ──────────────────────────── */

/** Compare Aizerman Hurwitz sector with circle criterion sector.
 *  Returns the ratio: (sector from circle criterion) / (Hurwitz sector).
 *  Ratio < 1 means circle criterion is more conservative. */
double abs_aizerman_compare_circle(const AbsLureSystem *sys);

/** Compare Aizerman with Popov criterion sector.
 *  Returns ratio for time-invariant nonlinearities. */
double abs_aizerman_compare_popov(const AbsLureSystem *sys);

/** Compute the conservatism gap: difference between Hurwitz sector
 *  bound and best frequency-domain criterion bound. */
void abs_aizerman_conservatism_gap(const AbsLureSystem *sys,
                                    double *circle_gap, double *popov_gap);

/* ── Structural Sufficient Conditions ─────────────────────────────── */

/** Check if the system satisfies conditions where Aizerman holds:
 *   1. n = 1 or n = 2 (always holds)
 *   2. A is symmetric (circle criterion proves it for alpha >= 0)
 *   3. The system is positive (Metzler matrix)
 *   4. The transfer function is minimum phase */
bool abs_aizerman_trivial_cases(const AbsLureSystem *sys);

/** Print Aizerman test result. */
void abs_aizerman_print_result(const AbsAizermanResult *res);

#endif /* ABS_AIZERMAN_H */
