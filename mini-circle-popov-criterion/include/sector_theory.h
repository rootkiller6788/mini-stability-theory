#ifndef SECTOR_THEORY_H
#define SECTOR_THEORY_H

/*
 * sector_theory.h — Sector-Bounded Nonlinearities Theory
 *
 * Sector conditions are the mathematical foundation of absolute stability
 * theory. A memoryless nonlinearity phi : R -> R satisfies the sector
 * condition [alpha, beta] if:
 *
 *   alpha * sigma^2  <=  sigma * phi(sigma)  <=  beta * sigma^2
 *
 * This module provides:
 *   1. Sector condition verification for arbitrary functions.
 *   2. Sector transformations: loop shifting, multiplier scaling.
 *   3. Sector arithmetic: intersection, union, scaling.
 *   4. Passivity <-> sector [0, inf) equivalence.
 *   5. Common nonlinearity characterizations (saturation, dead-zone, etc).
 *   6. Sector-based robust stability margins.
 *
 * References:
 *   - Zames, IEEE TAC, 1966 (Sector conditions)
 *   - Desoer & Vidyasagar, Feedback Systems, 1975
 *   - Safonov, Stability and Robustness, 1980
 */

#include <stdbool.h>
#include "circle_criterion.h"

/* ---- Nonlinearity Function Types ---- */
typedef double (*NonlinearityFunc)(double sigma);

typedef struct {
    NonlinearityFunc phi;
    NonlinearityFunc dphi;
    SectorBounds     sector;
    const char      *name;
} Nonlinearity;

/* ---- Nonlinearity Lifecycle ---- */
Nonlinearity* nonlinearity_create(NonlinearityFunc phi, NonlinearityFunc dphi,
                                   const char *name);
void nonlinearity_free(Nonlinearity* nl);

/* ---- Sector Condition Verification ---- */
bool sector_check_point(SectorBounds sb, double sigma, double phi);
bool sector_verify_function(NonlinearityFunc phi, SectorBounds sb,
                             double sigma_min, double sigma_max, int n);
SectorBounds sector_tightest_bounds(NonlinearityFunc phi,
                                     double sigma_min, double sigma_max, int n);

/* ---- Sector Transformations ---- */
SectorBounds sector_shift(SectorBounds sb, double c);
SectorBounds sector_scale(SectorBounds sb, double k);
SectorBounds sector_invert(SectorBounds sb);
SectorBounds sector_normalize(SectorBounds sb, double *k_out);

/* ---- Sector Arithmetic ---- */
SectorBounds sector_intersection(SectorBounds a, SectorBounds b);
SectorBounds sector_union(SectorBounds a, SectorBounds b);
bool sector_contains(SectorBounds outer, SectorBounds inner);

/* ---- Standard Nonlinearities ---- */
double nonlinearity_saturation(double sigma, double limit);
double nonlinearity_saturation_deriv(double sigma, double limit);
double nonlinearity_deadzone(double sigma, double delta);
double nonlinearity_deadzone_deriv(double sigma, double delta);
double nonlinearity_relay(double sigma);
double nonlinearity_quantizer(double sigma, double q_level);
double nonlinearity_cubic(double sigma, double c);

/* ---- Passivity and Sector Equivalence ---- */
bool sector_is_passive(SectorBounds sb);
bool sector_is_finite_gain(SectorBounds sb, double *gain_bound);

/* ---- Combined Nonlinearity Tests ---- */
SectorBounds sector_parallel(SectorBounds s1, SectorBounds s2);
SectorBounds sector_feedback(SectorBounds forward, SectorBounds feedback);

/* ---- Robust Stability via Sectors ---- */
double sector_max_allowed_gain(const TransferFunction* tf, int method, double *margin_at_max);

/* ---- Constants ---- */
#define SECTOR_SAMPLES   200
#define SECTOR_SIGMA_MAX 10.0
#define SECTOR_EPS       1e-10

#endif /* SECTOR_THEORY_H */

/* ---- Additional Nonlinearity Types ---- */
double nonlinearity_hysteresis_relay(double sigma, double hysteresis, double *prev_output);
double nonlinearity_sigmoid(double sigma, double gain);
double nonlinearity_polynomial(double sigma, const double *coeffs, int degree);

/* ---- Lyapunov Function Construction ---- */
double sector_lyapunov_value(const double *P, int n, const double *x,
                              NonlinearityFunc phi, const double *C, double alpha);

/* ---- Maximal Sector via Binary Search ---- */
double sector_max_beta_binary_search(const TransferFunction* tf,
                                      double alpha, int method);

/* ---- Conjecture Tests ---- */
int sector_test_aizerman(const TransferFunction* tf, SectorBounds sb);
int sector_test_kalman(const TransferFunction* tf, SectorBounds sb);

/* ---- Mixed Sector + Slope Restriction ---- */
bool sector_check_slope_restricted(SectorBounds sb, double sigma,
                                    double phi, double dphi_dsigma);

/* ---- Sector Sensitivity Analysis ---- */
double sector_sensitivity(const TransferFunction* tf, SectorBounds sb, double delta);
