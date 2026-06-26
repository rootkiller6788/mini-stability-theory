#ifndef ABS_CIRCLE_H
#define ABS_CIRCLE_H

/*
 * abs_circle.h — Circle Criterion for absolute stability
 *
 * For a Lur'e system with phi in sector [alpha, beta], the circle
 * criterion states:
 *
 * Case 0 < alpha < beta:
 *   The Nyquist plot of G(jw) must not intersect the disk D(alpha,beta)
 *   with center at -(alpha+beta)/(2*alpha*beta) on the real axis and
 *   radius (beta-alpha)/(2*alpha*beta), and must not encircle it.
 *
 * Case 0 = alpha < beta (infinite sector):
 *   The Nyquist plot of G(jw) must lie strictly to the right of the
 *   vertical line Re(s) = -1/beta. Equivalent to:
 *     Re[G(jw)] > -1/beta   for all w.
 *
 * Case alpha < 0 < beta:
 *   The Nyquist plot must lie strictly inside the disk.
 *
 * The circle criterion is a sufficient condition for absolute stability.
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "abs_core.h"

/* ── Nyquist Disk ─────────────────────────────────────────────────── */

typedef struct {
    double center_real;  /* real part of disk center (imag = 0) */
    double radius;       /* disk radius */
    bool   is_left_half; /* true if disk is left-half-plane (alpha*beta>0) */
    bool   is_degenerate; /* true if alpha=beta (trivial linear case) */
} AbsCircleDisk;

/** Compute the Nyquist disk from sector bounds.
 *  Handles all cases: 0<alpha<beta, 0=alpha<beta, alpha<0<beta. */
AbsCircleDisk abs_circle_compute_disk(double alpha, double beta);

/** Check if a complex point z lies inside or on the Nyquist disk. */
bool abs_circle_point_in_disk(AbsComplex z, AbsCircleDisk disk);

/* ── Circle Criterion Verification ────────────────────────────────── */

/** Result of circle criterion test. */
typedef struct {
    bool   is_stable;        /* true if circle criterion confirms stability */
    double min_margin;       /* minimum distance from Nyquist to disk boundary
                                (positive = safely outside) */
    double critical_freq;    /* frequency where margin is minimal */
    int    num_freqs_tested; /* number of frequency points evaluated */
    double *nyquist_re;      /* [num_freqs] real parts of G(jw) */
    double *nyquist_im;      /* [num_freqs] imag parts of G(jw) */
    double *freqs;           /* [num_freqs] frequency points */
    int    n_freqs;          /* length of above arrays */
} AbsCircleResult;

/** Run the circle criterion test on a Lur'e system.
 *
 *  Evaluates G(jw) at logarithmically-spaced frequencies from
 *  w_min to w_max (n_freqs points), and checks the Nyquist
 *  disk non-intersection condition.
 *
 *  @param sys       Lur'e system description
 *  @param n_freqs   number of frequency points (0 for default=200)
 *  @param w_min     minimum frequency (0 for auto = 1e-3)
 *  @param w_max     maximum frequency (0 for auto = 1e3)
 *  @return allocated result, caller must abs_circle_result_free() */
AbsCircleResult* abs_circle_test(const AbsLureSystem *sys,
                                  int n_freqs, double w_min, double w_max);

/** Free circle criterion result. */
void abs_circle_result_free(AbsCircleResult *res);

/** Quick pass/fail check without detailed frequency sweep.
 *  Checks DC gain (w=0) and infinite frequency gain, plus
 *  a few strategic frequency points. */
bool abs_circle_quick_check(const AbsLureSystem *sys);

/* ── Alternative Formulations ─────────────────────────────────────── */

/** Circle criterion via the equivalent frequency-domain inequality:
 *    Re[ (1 + beta*G(jw)) / (1 + alpha*G(jw)) ] > 0   for all w
 *  Evaluates at n_freqs points, returns minimum value. */
double abs_circle_freq_inequality(const AbsLureSystem *sys,
                                   int n_freqs, double w_min, double w_max);

/** Off-axis circle criterion (for complex sectors).
 *  Handles the case where the disk may be offset from the real axis.
 *  Returns true if absolute stability holds. */
bool abs_circle_off_axis_test(const AbsLureSystem *sys,
                               double gamma_real, double gamma_imag);

/** Multi-loop circle criterion for MIMO Lur'e systems.
 *  Uses the multivariable generalization:
 *    (I + beta*G(jw))^H * (I + alpha*G(jw)) > 0.
 *  Simplified for diagonal sector nonlinearities.
 *  Returns true if condition holds at all tested frequencies. */
bool abs_circle_multiloop_check(const double *A, int n,
                                 const double *B, int m,
                                 const double *C, int p,
                                 double alpha, double beta,
                                 int n_freqs);

/* ── Utilities ────────────────────────────────────────────────────── */

/** Print the circle criterion result to stdout. */
void abs_circle_print_result(const AbsCircleResult *res);

/** Generate Nyquist plot data for a given transfer function.
 *  Fills pre-allocated arrays re[n_freqs], im[n_freqs]. */
void abs_circle_nyquist_sweep(const AbsLureSystem *sys,
                               double *re, double *im,
                               double *freqs, int n_freqs,
                               double w_min, double w_max);

/** Find critical frequency where Nyquist is closest to critical point. */
double abs_circle_find_critical_freq(const AbsLureSystem *sys,
                                      double w_min, double w_max);

#endif /* ABS_CIRCLE_H */
