#ifndef ABS_POPOV_H
#define ABS_POPOV_H

/*
 * abs_popov.h — Popov Criterion for absolute stability
 *
 * The Popov criterion is a frequency-domain sufficient condition
 * for absolute stability of Lur'e systems with TIME-INVARIANT
 * nonlinearities phi(y) in sector [0, k].
 *
 * Criterion: exists q >= 0 such that for all w >= 0:
 *   Re[(1 + j*w*q) * G(jw)] + 1/k > 0
 *
 * Geometrically: the modified Nyquist plot
 *   X(w) = Re[G(jw)],  Y(w) = w * Im[G(jw)]
 * must lie strictly to the right of the Popov line
 * passing through (-1/k, 0) with slope 1/q.
 *
 * The Popov criterion is generally less conservative than the
 * circle criterion for time-invariant nonlinearities.
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "abs_core.h"

/* ── Popov Plot Point ─────────────────────────────────────────────── */

typedef struct {
    double X;     /* Re[G(jw)] */
    double Y;     /* w * Im[G(jw)] */
    double w;     /* frequency */
} AbsPopovPoint;

/* ── Popov Line ───────────────────────────────────────────────────── */

typedef struct {
    double slope;     /* 1/q, the Popov line slope */
    double intercept; /* Popov line X-intercept = -1/k */
    double q;         /* Popov multiplier >= 0 */
    double k;         /* sector upper bound */
} AbsPopovLine;

/** Compute the Popov line for given sector bound k.
 *  Line passes through (-1/k, 0) with variable slope. */
AbsPopovLine abs_popov_line_make(double k);

/** Check if a Popov point (X, Y) lies strictly to the right
 *  of the Popov line with given slope (1/q). */
bool abs_popov_point_right_of_line(AbsPopovPoint pt, AbsPopovLine line);

/* ── Popov Criterion Results ──────────────────────────────────────── */

typedef enum {
    ABS_POPOV_PASS       = 0,  /* criterion satisfied */
    ABS_POPOV_FAIL       = 1,  /* no q found satisfying condition */
    ABS_POPOV_INFEASIBLE = 2,  /* sector bound too large */
    ABS_POPOV_DEGENERATE = 3   /* trivial case */
} AbsPopovStatus;

typedef struct {
    AbsPopovStatus status;
    double optimal_q;       /* q that gives maximum margin */
    double optimal_slope;   /* 1/q for optimal q */
    double min_margin;      /* minimum of Re[G] + q*w*Im[G] + 1/k */
    double critical_freq;   /* frequency at minimum margin */
    int    n_freqs_tested;
    double *popov_X;        /* [n_freqs] Popov X coordinates */
    double *popov_Y;        /* [n_freqs] Popov Y coordinates */
    double *freqs;          /* [n_freqs] frequency points */
    int    n_freqs;
} AbsPopovResult;

/** Run the Popov criterion test on a Lur'e system.
 *
 *  For sector [0, k], searches for q >= 0 such that the Popov
 *  frequency condition holds. Computes Popov plot and finds
 *  the optimal multiplier via tangent-line method.
 *
 *  @param sys       Lur'e system (sector must be [0, k] or [alpha, beta])
 *  @param n_freqs   number of frequency points (0 for default=200)
 *  @param w_min     minimum frequency (0 for auto)
 *  @param w_max     maximum frequency (0 for auto)
 *  @return allocated result */
AbsPopovResult* abs_popov_test(const AbsLureSystem *sys,
                                int n_freqs, double w_min, double w_max);

/** Free Popov result. */
void abs_popov_result_free(AbsPopovResult *res);

/** Quick Popov check without detailed sweep.
 *  Tests a few strategic q values and frequency points. */
AbsPopovStatus abs_popov_quick_check(const AbsLureSystem *sys);

/* ── Popov Criterion Variants ─────────────────────────────────────── */

/** Generalized Popov criterion for sector [alpha, beta]:
 *   Uses loop transformation to map to [0, k] form,
 *   then applies standard Popov criterion.
 *   Returns k_max, the maximum sector bound confirmed. */
double abs_popov_generalized(const AbsLureSystem *sys,
                              int n_freqs, double w_min, double w_max);

/** Popov criterion via Lyapunov function construction.
 *  Builds V(x) = x^T*P*x + q*integral_0^{c^T*x} phi(sigma) dsigma
 *  Solves the LMI to verify V_dot < 0.
 *  Returns true if a valid P, q pair is found. */
bool abs_popov_lyapunov_synthesis(const AbsLureSystem *sys,
                                   double *P_out, double *q_out);

/** Multiplier method: generalized Popov with dynamic multipliers.
 *  Uses Zames-Falb multipliers or Popov multiplier M(s) = 1 + q*s.
 *  Returns the maximum stability margin. */
double abs_popov_multiplier_method(const double *A, const double *b,
                                    const double *c, int n,
                                    double k, double *optimal_q);

/** Frequency-domain Popov inequality check.
 *  For a given q, evaluates:
 *    Re[(1 + j*w*q)*G(jw)] + 1/k  at n_freqs points.
 *  Returns the minimum value (positive = stable). */
double abs_popov_eval_at_q(const AbsLureSystem *sys, double q,
                            int n_freqs, double w_min, double w_max);

/** Find optimal Popov multiplier q* that maximizes stability margin.
 *  Uses golden-section search on q >= 0.
 *  Returns optimal q, stores margin in *margin_out. */
double abs_popov_find_optimal_q(const AbsLureSystem *sys,
                                 int n_freqs, double w_min, double w_max,
                                 double *margin_out);

/* ── Utilities ────────────────────────────────────────────────────── */

/** Generate Popov plot data: fills X[n_freqs], Y[n_freqs], freqs[n_freqs]. */
void abs_popov_plot_data(const AbsLureSystem *sys,
                          double *X, double *Y, double *freqs,
                          int n_freqs, double w_min, double w_max);

/** Print Popov result */
void abs_popov_print_result(const AbsPopovResult *res);

/** Compute the maximum k for which Popov criterion guarantees stability.
 *  Uses bisection on k. */
double abs_popov_max_stable_k(const AbsLureSystem *sys,
                               int n_freqs, double w_min, double w_max);

#endif /* ABS_POPOV_H */
