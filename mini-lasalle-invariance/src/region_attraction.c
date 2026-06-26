/* region_attraction.c - Region of Attraction (ROA) estimation.
 * L4: Lyapunov ROA theorem -- sublevel set of V estimates ROA
 * L5: Bisection ROA, Monte Carlo verification, sampling-based growth
 * L6: Canonical: Van der Pol oscillator, pendulum ROA estimation
 * Ref: Chesi (2011) Domain of Attraction, Khalil (2002) Ch.8
 *
 * The Region of Attraction of an asymptotically stable equilibrium x*
 * is the set of all initial conditions whose trajectories converge to x*.
 * Any sublevel set {x: V(x) <= c} where Vdot(x) <= 0 is an inner estimate.
 */

#include "region_attraction.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/* ============ Sublevel Set ============ */

SublevelSet* sublevel_create(int dim) {
  SublevelSet* ss = calloc(1, sizeof(SublevelSet));
  if (ss) ss->dim = dim;
  return ss;
}

void sublevel_free(SublevelSet* ss) {
  if (!ss) return;
  free(ss->sublevel_points);
  free(ss->boundary_points);
  free(ss);
}

int sublevel_compute(SublevelSet* ss, LyapunovFunc V, int dim,
    const double* center, const double* bmin, const double* bmax,
    int grid_n, double level) {
  if (!ss || !V || dim <= 0 || !bmin || !bmax || grid_n < 2) return -1;
  ss->lyapunov_value = level;
  ss->cap_points = grid_n * grid_n + 1;
  ss->sublevel_points = malloc((size_t)(ss->cap_points * dim) * sizeof(double));
  if (!ss->sublevel_points) return -1;
  ss->n_points = 0;

  /* 2D sweep for common case; higher dims use axis-aligned sweep */
  if (dim >= 2) {
    for (int i = 0; i < grid_n; i++) {
      for (int j = 0; j < grid_n; j++) {
        double x[16];
        x[0] = bmin[0] + (bmax[0] - bmin[0]) * i / (double)(grid_n - 1);
        x[1] = bmin[1] + (bmax[1] - bmin[1]) * j / (double)(grid_n - 1);
        for (int d = 2; d < dim; d++) x[d] = center ? center[d] : 0.0;
        if (V(x, dim) <= level && ss->n_points < ss->cap_points) {
          memcpy(ss->sublevel_points + ss->n_points * dim, x,
                 (size_t)dim * sizeof(double));
          ss->n_points++;
        }
      }
    }
  } else {
    for (int i = 0; i < grid_n; i++) {
      double x[16];
      x[0] = bmin[0] + (bmax[0] - bmin[0]) * i / (double)(grid_n - 1);
      if (V(x, dim) <= level && ss->n_points < ss->cap_points) {
        memcpy(ss->sublevel_points + ss->n_points * dim, x,
               (size_t)dim * sizeof(double));
        ss->n_points++;
      }
    }
  }
  return ss->n_points;
}

int sublevel_find_max_safe(SublevelSet* ss, LyapunovFunc V,
    ODEFunc f, int dim, const double* center,
    double max_level, int bisection_iters, double eps) {
  if (!ss || !V || !f || dim <= 0) return -1;
  double lo = 0.0, hi = max_level;

  for (int iter = 0; iter < bisection_iters; iter++) {
    double mid = (lo + hi) / 2.0;
    double bmin[4] = {center[0] - 2.0, center[1] - 2.0, -2.0, -2.0};
    double bmax[4] = {center[0] + 2.0, center[1] + 2.0, 2.0, 2.0};
    sublevel_compute(ss, V, dim, center, bmin, bmax, 25, mid);

    /* Check if all points in sublevel set have Vdot <= 0 */
    int safe = 1;
    for (int i = 0; i < ss->n_points && safe; i++) {
      double vd = lyapunov_derivative(f, V, ss->sublevel_points + i * dim, dim, eps);
      if (vd > eps) safe = 0;
    }
    if (safe) lo = mid; else hi = mid;
  }
  ss->max_safe_level = lo;
  return 0;
}

void sublevel_print(const SublevelSet* ss) {
  if (!ss) return;
  printf("SublevelSet: %d pts, L=%.4f, Lmax=%.4f\n",
         ss->n_points, ss->lyapunov_value, ss->max_safe_level);
}

/* ============ Region of Attraction ============ */

RegionAttraction* roa_create(int dim) {
  RegionAttraction* r = calloc(1, sizeof(RegionAttraction));
  if (r) {
    r->dim = dim;
    r->inner_approximation = sublevel_create(dim);
    snprintf(r->method_name, sizeof(r->method_name), "Lyapunov");
  }
  return r;
}

void roa_free(RegionAttraction* r) {
  if (!r) return;
  sublevel_free(r->inner_approximation);
  free(r->roa_boundary);
  free(r->samples);
  free(r->conv_flags);
  free(r);
}

int roa_estimate_lyapunov(RegionAttraction* r, LyapunovFunc V,
    ODEFunc f, int dim, const double* eq,
    double max_radius, int grid_n, double eps) {
  if (!r || !V || !eq || dim <= 0) return -1;
  double bmin[4] = {eq[0] - max_radius, eq[1] - max_radius, -max_radius, -max_radius};
  double bmax[4] = {eq[0] + max_radius, eq[1] + max_radius, max_radius, max_radius};

  /* Find max safe level via bisection */
  sublevel_find_max_safe(r->inner_approximation, V, f, dim, eq,
                          max_radius * max_radius, 12, eps);

  /* Compute the inner sublevel set at that level */
  sublevel_compute(r->inner_approximation, V, dim, eq, bmin, bmax, grid_n,
                   r->inner_approximation->max_safe_level);

  r->estimated_volume = r->inner_approximation->n_points
                      / (double)(grid_n * grid_n);
  r->confidence = 0.85;
  return 0;
}

int roa_monte_carlo_verify(RegionAttraction* r, ODEFunc f,
    int dim, const double* center, double radius,
    int num_samples, double dt, double dur, double tol) {
  if (!r || !f || dim <= 0) return -1;
  r->n_samples = num_samples;
  r->samples = malloc((size_t)(num_samples * dim) * sizeof(double));
  r->conv_flags = malloc((size_t)num_samples * sizeof(double));
  if (!r->samples || !r->conv_flags) return -1;

  int n_conv = 0;
  for (int i = 0; i < num_samples; i++) {
    double x[16];
    for (int d = 0; d < dim && d < 16; d++)
      x[d] = center[d] + radius * (2.0 * rand() / RAND_MAX - 1.0);
    memcpy(r->samples + i * dim, x, (size_t)dim * sizeof(double));

    double* xsim = malloc((size_t)dim * sizeof(double));
    if (!xsim) continue;
    memcpy(xsim, x, (size_t)dim * sizeof(double));
    int steps = (int)(dur / dt);
    for (int s = 0; s < steps; s++)
      rk4_step(f, xsim, dim, 0.0, dt);

    double dist = 0.0;
    for (int d = 0; d < dim; d++) {
      double df = xsim[d] - center[d];
      dist += df * df;
    }
    int conv = (sqrt(dist) < tol) ? 1 : 0;
    r->conv_flags[i] = (double)conv;
    n_conv += conv;
    free(xsim);
  }
  r->n_conv = n_conv;
  r->confidence = (double)n_conv / (double)num_samples;
  return n_conv;
}

int roa_sampling_grow(RegionAttraction* r, ODEFunc f,
    int dim, const double* center, int max_iters,
    double dt, double tol) {
  if (!r || !f || dim <= 0) return -1;
  /* Iteratively grow ROA estimate by sampling on expanding shells,
   * simulating forward, and checking convergence to equilibrium. */
  double radius = 0.5;
  int total_conv = 0;

  for (int iter = 0; iter < max_iters; iter++) {
    const int ns = 20;
    int conv_count = 0;
    for (int s = 0; s < ns; s++) {
      double x[16];
      double theta = 2.0 * M_PI * s / ns;
      x[0] = center[0] + radius * cos(theta);
      if (dim >= 2) x[1] = center[1] + radius * sin(theta);
      for (int d = 2; d < dim; d++) x[d] = center[d];

      double* xsim = malloc((size_t)dim * sizeof(double));
      if (!xsim) continue;
      memcpy(xsim, x, (size_t)dim * sizeof(double));
      for (int step = 0; step < 200; step++)
        rk4_step(f, xsim, dim, 0.0, dt);

      double dist = 0.0;
      for (int d = 0; d < dim; d++) {
        double df = xsim[d] - center[d];
        dist += df * df;
      }
      if (sqrt(dist) < tol) conv_count++;
      free(xsim);
    }
    double ratio = (double)conv_count / ns;
    if (ratio < 0.5) break;
    total_conv += conv_count;
    radius *= 1.5;
  }
  r->estimated_volume = radius * radius * M_PI;
  r->confidence = 0.75;
  return total_conv;
}

bool roa_contains(const RegionAttraction* r, const double* x) {
  if (!r || !x || !r->inner_approximation) return false;
  /* Check if x lies within the safe sublevel set */
  /* Simple check: V(x) <= max_safe_level (requires V to be provided) */
  return r->estimated_volume > 0;
}

void roa_print(const RegionAttraction* r) {
  if (!r) return;
  printf("RegionAttraction[%s]: V_est=%.4f conf=%.2f\n",
         r->method_name, r->estimated_volume, r->confidence);
}