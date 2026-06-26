#ifndef POPOV_CRITERION_H
#define POPOV_CRITERION_H

/*
 * popov_criterion.h — Popov Criterion for Absolute Stability
 *
 * The Popov criterion is a frequency-domain test for absolute stability of
 * Lur'e systems with a TIME-INVARIANT nonlinearity in sector [0, k] (or
 * shifted to [alpha, beta] via loop transformation).
 *
 * Unlike the circle criterion (which works for time-varying nonlinearities),
 * the Popov criterion exploits time-invariance through the Popov multiplier
 * (1 + eta*s), leading to a less conservative test.
 *
 * Popov Condition:
 *   There exists eta >= 0 such that for all w >= 0:
 *     Re[(1 + j*w*eta) * G(jw)] + 1/k > 0
 *
 * Geometrically (Popov Plot):
 *   Plot Re[G(jw)] on x-axis, w*Im[G(jw)] on y-axis.
 *   If there exists a line through (-1/k, 0) with slope 1/eta such that
 *   the entire Popov plot lies strictly to the right of this line,
 *   then the system is absolutely stable.
 *
 * References:
 *   - Popov, Automatika i Telemekhanika, 1961
 *   - Khalil, Nonlinear Systems (3rd ed.), Section 7.3
 *   - Haddad & Chellaboina, Nonlinear Dynamical Systems, 2008
 */

#include "circle_criterion.h"

/* ---- Popov Plot Data ---- */
/* Popov plot maps: x = Re[G(jw)], y = w * Im[G(jw)] */
typedef struct {
    double *omega_re;   /* x-axis: Re[G(jw)]              */
    double *omega_im;   /* y-axis: w * Im[G(jw)]          */
    int     n;          /* number of frequency points     */
} PopovData;

/* Allocate Popov data storage for n frequency points.
 * Must be freed with popov_data_free(). */
PopovData* popov_data_create(int n);

/* Free Popov data and its internal arrays. */
void popov_data_free(PopovData* pd);

/* Compute the Popov plot: Re[G(jw)] and w*Im[G(jw)] for
 * frequencies from wmin to wmax, linearly spaced. */
void popov_compute(const TransferFunction* tf, PopovData* pd,
                   double wmin, double wmax);

/* ---- Popov Criterion Core Functions ---- */

/* Check the Popov criterion for a given transfer function and sector,
 * with a specified multiplier parameter eta.
 *
 * Parameters:
 *   tf     - transfer function G(s) of the linear part
 *   sb     - sector bounds [alpha, beta] for the nonlinearity
 *   eta    - Popov multiplier parameter (eta >= 0)
 *   margin - output: distance margin (positive => stable)
 *
 * Returns 1 if Popov condition holds, 0 otherwise.
 *
 * The sector [alpha, beta] is internally shifted to [0, k] via
 * loop transformation before applying the Popov test.
 */
int popov_criterion_check(const TransferFunction* tf, SectorBounds sb,
                           double eta, double *margin);

/* Find the optimal Popov multiplier eta that maximizes the stability
 * margin. Searches over a grid of eta values.
 *
 * Parameters:
 *   tf        - transfer function
 *   sb        - sector bounds
 *   eta_opt   - output: optimal eta value found
 *   margin_opt- output: stability margin at optimal eta
 */
void popov_find_optimal_slope(const TransferFunction* tf,
                               SectorBounds sb,
                               double *eta_opt, double *margin_opt);

/* ---- Popov Line Test ---- */
/* Given a Popov plot and a line y = slope * x + intercept, test if
 * the entire plot lies to the left of the line (for the Popov condition:
 * plot must lie to the RIGHT, so we check y < slope*x + intercept for all
 * points when testing violation).
 *
 * For the Popov criterion: line passes through (-1/k, 0).
 *   intercept = slope * (1/k)  => intercept = slope/k
 * But in canonical form: line is y = (1/eta)*(x + 1/k)
 *   slope = 1/eta, intercept = 1/(eta*k)
 *
 * Returns true if all points satisfy the line condition.
 * On violation, *violation = max(line_y - data_y). */
bool popov_line_test(const PopovData* pd, double slope, double intercept,
                     double *violation);

/* ---- Multiplier Search ---- */
/* Grid search for the optimal Popov multiplier eta.
 * Tests eta values from eta_min to eta_max in n_steps.
 * Returns the eta giving the largest stability margin. */
double popov_grid_search_eta(const TransferFunction* tf,
                              SectorBounds sb,
                              double eta_min, double eta_max,
                              int n_steps, double *best_margin);

/* ---- Sector Transformation ---- */
/* Apply loop transformation to shift a general sector [alpha, beta]
 * to the canonical sector [0, k].
 *
 * Given: G(s) with phi in [alpha, beta]
 * Transformed: G'(s) = G(s) / (1 + alpha*G(s))
 *              k = beta - alpha
 *
 * Returns the transformed transfer function (caller must free).
 */
TransferFunction* popov_loop_transform(const TransferFunction* tf,
                                        SectorBounds sb, double *k);

/* ---- Comparison ---- */
/* Compare the circle criterion and Popov criterion on the same system.
 * Returns: 0 = neither, 1 = circle only, 2 = popov only, 3 = both pass.
 * Demonstrates that Popov is less conservative for time-invariant phi. */
int popov_compare_criteria(const TransferFunction* tf, SectorBounds sb,
                            double *circle_margin, double *popov_margin,
                            double *eta_found);

/* ---- Constants ---- */
#define POPOV_MAX_FREQ      1000.0  /* max frequency for Popov sweep    */
#define POPOV_N_POINTS       500    /* default frequency resolution     */
#define POPOV_ETA_MAX        100.0  /* max multiplier search range      */
#define POPOV_ETA_STEPS       50    /* grid search steps for eta        */

#endif /* POPOV_CRITERION_H */
