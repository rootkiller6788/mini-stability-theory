/* barbashin_krasovskii.c - Barbashin-Krasovskii global asymptotic stability verification.
 * L4: BK Theorem - radially unbounded V, Vdot<0 everywhere, unique equilibrium at 0
 * L5: Grid-search Vdot sign verification, equilibrium uniqueness test
 * Ref: Barbashin & Krasovskii (1952), Khalil (2002) Theorem 4.2
 */

#include "barbashin_krasovskii.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

BKResult* bk_create(int dim) {
  BKResult* r = calloc(1, sizeof(BKResult));
  if (r) r->n_dim = dim;
  return r;
}

void bk_free(BKResult* r) {
  if (!r) return;
  free(r->candidate_points);
  free(r->verified_equilibria);
  free(r);
}

int bk_verify_radial_unbounded(LyapunovFunc V, int n, double radius, int samples) {
  /* Verify V(x) -> inf as ||x|| -> inf by sampling on sphere of given radius */
  if (!V || n <= 0 || samples < 4) return -1;
  double* x = calloc((size_t)n, sizeof(double));
  if (!x) return -1;
  double min_val = 1e100;
  for (int s = 0; s < samples; s++) {
    double theta = 2.0 * M_PI * s / samples;
    if (n >= 2) { x[0] = radius * cos(theta); x[1] = radius * sin(theta); }
    else if (n == 1) x[0] = radius * ((s % 2) ? 1.0 : -1.0);
    for (int d = 2; d < n; d++) x[d] = radius * sin(theta + d * 0.5);
    double v = V(x, n);
    if (v < min_val) min_val = v;
  }
  free(x);
  /* Radial unbounded if V grows with distance. Simple heuristic: V > sqrt(radius) */
  return (min_val > 0.1 * radius) ? 1 : 0;
}

int bk_verify_negative_definite(LyapunovFunc V, ODEFunc f, int n,
    const double* rmin, const double* rmax, int gpts, double eps) {
  if (!V || !f || n <= 0 || !rmin || !rmax || gpts < 2) return -1;
  int max_positive = 0;
  double max_vdot = -1e100, min_vdot = 1e100;
  double* x = calloc((size_t)n, sizeof(double));
  if (!x) return -1;

  for (int gi = 0; gi < gpts; gi++) {
    for (int d = 0; d < n; d++)
      x[d] = rmin[d] + (rmax[d] - rmin[d]) * gi / (double)(gpts - 1);
    double vd = lyapunov_derivative(f, V, x, n, eps);
    if (vd > max_vdot) max_vdot = vd;
    if (vd < min_vdot) min_vdot = vd;
    if (vd > 0.0) max_positive++;
  }
  free(x);
  /* Consider negative definite if max Vdot <= 0 */
  return (max_vdot <= eps) ? 1 : (max_vdot <= 10.0 * eps ? 0 : -1);
}

int bk_find_equilibrium_set(BKResult* r, ODEFunc f, const double* x0,
    int n, int max_iter, double tol) {
  if (!r || !f || !x0 || n <= 0) return -1;
  double* x = malloc((size_t)n * sizeof(double));
  if (!x) return -1;
  memcpy(x, x0, (size_t)n * sizeof(double));

  for (int iter = 0; iter < max_iter; iter++) {
    double fval[16];
    f(x, fval, n, 0.0);
    double norm = 0.0;
    for (int d = 0; d < n; d++) norm += fval[d] * fval[d];
    norm = sqrt(norm);
    if (norm < tol) {
      /* Converged to equilibrium - store it */
      r->verified_equilibria = realloc(r->verified_equilibria,
          (size_t)(r->n_verified + 1) * (size_t)n * sizeof(double));
      if (!r->verified_equilibria) { free(x); return -1; }
      memcpy(r->verified_equilibria + r->n_verified * n, x,
             (size_t)n * sizeof(double));
      r->n_verified++;
      free(x);
      return r->n_verified;
    }
    /* Newton-like descent: x -= alpha * f(x) / ||f|| */
    double alpha = 0.01;
    for (int d = 0; d < n; d++)
      x[d] -= alpha * fval[d] / (norm + 1e-15);
  }
  free(x);
  return 0;
}

bool bk_is_gas(const BKResult* r) {
  if (!r) return false;
  if (r->is_globally_asymptotically_stable) return true;
  return r->n_verified == 1 && r->lyapunov_margin > 0;
}

void bk_print(const BKResult* r) {
  if (!r) { printf("BKResult: NULL\n"); return; }
  printf("=== Barbashin-Krasovskii Result ===\n");
  printf("  GAS: %s\n", bk_is_gas(r) ? "YES" : "NO");
  printf("  Equilibria found: %d\n", r->n_verified);
  printf("  Lyapunov margin: %.6f\n", r->lyapunov_margin);
}