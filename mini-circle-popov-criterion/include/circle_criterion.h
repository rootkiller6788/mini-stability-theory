#ifndef CIRCLE_CRITERION_H
#define CIRCLE_CRITERION_H

/*
 * circle_criterion.h — Circle Criterion for Absolute Stability
 *
 * The Circle Criterion provides a sufficient condition for absolute stability
 * of a Lur'e system with a (possibly time-varying) nonlinearity in sector
 * [alpha, beta]. The criterion requires that the Nyquist plot of the linear
 * part G(jw) does not intersect or encircle a critical disk D(alpha, beta).
 *
 * Critical Disk:
 *   Center c = -0.5 * (1/alpha + 1/beta)
 *   Radius r =  0.5 * |1/alpha - 1/beta|
 *
 * If the Nyquist plot of G(jw) stays strictly outside D(alpha, beta) and
 * winds around it zero times (if the open-loop is stable), then the closed-loop
 * Lur'e system is absolutely stable for ALL nonlinearities in [alpha, beta].
 *
 * References:
 *   - Khalil, Nonlinear Systems (3rd ed.), Section 7.1
 *   - Zames, IEEE TAC, 1966
 *   - Vidyasagar, Nonlinear Systems Analysis, 1993
 */

#include <stdbool.h>
#include <stddef.h>

/* ---- Transfer Function Representation ---- */
/* Polynomial transfer function: G(s) = num(s)/den(s)
 * num[0] + num[1]*s + num[2]*s^2 + ...
 * den[0] + den[1]*s + den[2]*s^2 + ...  (den[0] typically 1.0)
 */
typedef struct {
    double *num;
    double *den;
    int     n_num;    /* number of numerator coefficients   */
    int     n_den;    /* number of denominator coefficients */
} TransferFunction;

/* Create a transfer function with given polynomial degrees.
 * Allocates and zero-initializes coefficient arrays.
 * den[0] is set to 1.0 by default. Caller must fill remaining coeffs. */
TransferFunction* tf_create(int n_num, int n_den);

/* Free a transfer function and its coefficient arrays. */
void tf_free(TransferFunction* tf);

/* Evaluate |G(j*omega)| (magnitude) at a given frequency. */
double tf_eval(const TransferFunction* tf, double omega);

/* Evaluate G(j*omega) = real + j*imag at a given frequency.
 * Uses direct polynomial evaluation with s = j*omega:
 *   (j*omega)^k = omega^k * (j^k)
 *   j^0=1, j^1=j, j^2=-1, j^3=-j, j^4=1, ...
 */
void tf_eval_complex(const TransferFunction* tf, double omega,
                     double *re, double *im);

/* ---- Nyquist Plot Data ---- */
typedef struct {
    double *real;   /* Re[G(jw)] at each frequency point */
    double *imag;   /* Im[G(jw)] at each frequency point */
    int     n;      /* number of frequency points        */
} NyquistData;

/* Allocate Nyquist data storage for n frequency points. */
NyquistData* nyquist_create(int n);

/* Free Nyquist data. */
void nyquist_free(NyquistData* nd);

/* Compute Nyquist plot data for frequencies from wmin to wmax.
 * Evenly spaced on linear frequency axis. */
void nyquist_compute(const TransferFunction* tf, NyquistData* nd,
                     double wmin, double wmax);

/* Log-spaced Nyquist for better low-frequency resolution. */
void nyquist_compute_log(const TransferFunction* tf, NyquistData* nd,
                          double wmin, double wmax);

/* Full Nyquist contour (positive and negative frequencies). */
void nyquist_compute_full(const TransferFunction* tf, NyquistData* nd,
                           double wmin, double wmax);

/* ---- Sector Bounds ---- */
typedef struct SectorBounds {
    double alpha;           /* lower sector bound                  */
    double beta;            /* upper sector bound (beta > alpha)   */
    int    is_time_varying; /* 1 = time-varying, 0 = time-invariant */
    int    is_slope_bounded;/* 1 = slope-restricted nonlinearity   */
    double slope_min;       /* minimum slope (if slope-restricted) */
    double slope_max;       /* maximum slope (if slope-restricted) */
} SectorBounds;

/* ---- Circle Criterion Core Functions ---- */

/* Compute the center of the critical disk on the real axis.
 * Formula: c = -0.5 * (1/alpha + 1/beta)
 * For alpha=0: c = -infinity (half-plane test)
 * For alpha=-beta: c = 0
 */
double circle_criterion_center(SectorBounds sb);

/* Compute the radius of the critical disk.
 * Formula: r = 0.5 * |1/alpha - 1/beta|
 * For alpha=0: r = infinity (half-plane)
 * For alpha=beta: r = 0 (point criterion)
 */
double circle_criterion_radius(SectorBounds sb);

/* Check the circle criterion for a given transfer function and sector.
 * Returns 1 if absolutely stable (criterion satisfied), 0 otherwise.
 * On return, *margin = r - min_distance (positive => stable).
 *
 * Algorithm:
 *   1. Compute critical disk center and radius.
 *   2. Generate Nyquist plot of G(jw).
 *   3. Check if Nyquist enters the disk: min_distance < r => violation.
 *   4. Check for encirclements (winding number != 0).
 */
int circle_criterion_check(const TransferFunction* tf, SectorBounds sb,
                            double *margin);

/* Check if the Nyquist plot encircles a given disk center.
 * Uses winding number computation via cumulative angle change.
 * center_x is the disk center (real axis), center_y is typically 0.
 * radius is provided for future use (currently only winding is checked).
 */
bool circle_criterion_encircles(const NyquistData* nd,
                                 double center_x, double center_y);

/* Generate (x,y) data pairs for plotting the Nyquist curve.
 * Convenience wrapper that fills arrays with Re[G] and Im[G]. */
void circle_criterion_plot_data(const TransferFunction* tf,
                                 SectorBounds sb,
                                 double *x, double *y, int n);

/* Compute the minimum distance from Nyquist plot to critical disk center.
 * Returns the minimum Euclidean distance across all frequency points. */
double circle_criterion_min_distance(const NyquistData* nd,
                                      double center_x, double center_y);

/* ---- Disk Intersection Test ---- */
/* Check if any point on the Nyquist plot lies inside the critical disk.
 * Returns true if the Nyquist curve enters or touches the disk. */
bool circle_criterion_in_disk(const NyquistData* nd,
                               double center_x, double radius);

/* ---- Special Cases ---- */
/* Handle the infinite-sector case (alpha = beta): reduces to linear
 * Nyquist criterion with gain margin interpretation. */
int circle_criterion_linear_case(const TransferFunction* tf,
                                  double gain, double *margin);

/* Handle the half-plane case (alpha = 0): requires Re[G(jw)] > -1/beta
 * for all frequencies. */
int circle_criterion_halfplane_case(const TransferFunction* tf,
                                     double beta, double *margin);

/* ---- Constants ---- */
#define CIRCLE_MAX_FREQ   100.0   /* maximum frequency for Nyquist sweep  */
#define CIRCLE_N_POINTS   1000    /* default number of frequency points   */
#define CIRCLE_EPS        1e-12   /* numerical tolerance                  */

/* ---- Crossover Frequency ---- */
/* Find the frequency where |G(jw)| = 1 (0 dB gain crossover). */
double circle_criterion_crossover_freq(const TransferFunction* tf);

#endif /* CIRCLE_CRITERION_H */
