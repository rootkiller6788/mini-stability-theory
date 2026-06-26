/* invariant_sets.c - Invariant set computation for stability analysis.
 * L2: Positively invariant, negatively invariant, omega-limit sets
 * L4: LaSalles largest invariant set theorem
 * L5: Grid-based invariance verification, omega-limit approximation
 * Ref: LaSalle (1960), Bhatia & Szego (1970)
 */

#include "invariant_sets.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

InvariantSet* iset_create(int dim) {
  InvariantSet* is = calloc(1, sizeof(InvariantSet));
  if (is) {
    is->dim = dim;
    is->capacity = 256;
    is->points = malloc((size_t)(256 * dim) * sizeof(double));
  }
  return is;
}

void iset_free(InvariantSet* is) {
  if (!is) return; free(is->points); free(is);
}

int iset_add_point(InvariantSet* is, const double* x) {
  if (!is || !x) return -1;
  if (is->n_points >= is->capacity) {
    is->capacity *= 2;
    double* np = realloc(is->points, (size_t)(is->capacity * is->dim) * sizeof(double));
    if (!np) return -1;
    is->points = np;
  }
  memcpy(is->points + is->n_points * is->dim, x, (size_t)is->dim * sizeof(double));
  return is->n_points++;
}

bool iset_check_positive_invariance(InvariantSet* is, ODEFunc f,
    double dt, double duration, double eps) {
  if (!is || !f || is->n_points == 0) return false;
  int steps = (int)(duration / dt);
  int dim = is->dim;
  for (int i = 0; i < is->n_points; i++) {
    double* x = malloc((size_t)dim * sizeof(double));
    if (!x) continue;
    memcpy(x, is->points + i * dim, (size_t)dim * sizeof(double));
    int stays = 1;
    for (int s = 0; s < steps && stays; s++) {
      rk4_step(f, x, dim, 0.0, dt);
      int found = 0;
      for (int j = 0; j < is->n_points; j++) {
        double dist = 0.0;
        for (int d = 0; d < dim; d++) {
          double diff = x[d] - is->points[j * dim + d];
          dist += diff * diff;
        }
        if (sqrt(dist) < eps) { found = 1; break; }
      }
      if (!found) stays = 0;
    }
    free(x);
    if (!stays) return false;
  }
  is->is_positively_invariant = true;
  return true;
}

bool iset_check_negative_invariance(InvariantSet* is, ODEFunc f,
    double dt, double duration, double eps) {
  if (!is || !f) return false;
  /* Negative invariance: backward integration stays in set */
  int steps = (int)(duration / dt);
  int dim = is->dim;
  for (int i = 0; i < is->n_points; i++) {
    double* x = malloc((size_t)dim * sizeof(double));
    if (!x) continue;
    memcpy(x, is->points + i * dim, (size_t)dim * sizeof(double));
    int stays = 1;
    for (int s = 0; s < steps && stays; s++) {
      rk4_step(f, x, dim, 0.0, -dt);
      int found = 0;
      for (int j = 0; j < is->n_points; j++) {
        double dist = 0.0;
        for (int d = 0; d < dim; d++) {
          double diff = x[d] - is->points[j * dim + d];
          dist += diff * diff;
        }
        if (sqrt(dist) < eps) { found = 1; break; }
      }
      if (!found) stays = 0;
    }
    free(x);
    if (!stays) { is->is_negatively_invariant = false; return false; }
  }
  is->is_negatively_invariant = true;
  is->is_invariant = is->is_positively_invariant && is->is_negatively_invariant;
  return true;
}

bool iset_find_fixed_points(const InvariantSet* is, ODEFunc f,
    double* fp, int max_fp, int* n_found, double eps) {
  if (!is || !f || !fp || !n_found) return false;
  *n_found = 0;
  for (int i = 0; i < is->n_points && *n_found < max_fp; i++) {
    double dx[16];
    f(is->points + i * is->dim, dx, is->dim, 0.0);
    double norm = 0.0;
    for (int d = 0; d < is->dim; d++) norm += dx[d] * dx[d];
    if (sqrt(norm) < eps) {
      memcpy(fp + (*n_found) * is->dim, is->points + i * is->dim,
             (size_t)is->dim * sizeof(double));
      (*n_found)++;
    }
  }
  return *n_found > 0;
}

double iset_diameter(const InvariantSet* is) {
  if (!is || is->n_points < 2) return 0.0;
  double md = 0.0;
  for (int i = 0; i < is->n_points; i++)
    for (int j = i + 1; j < is->n_points; j++) {
      double dist = 0.0;
      for (int d = 0; d < is->dim; d++) {
        double diff = is->points[i * is->dim + d] - is->points[j * is->dim + d];
        dist += diff * diff;
      }
      double d = sqrt(dist);
      if (d > md) md = d;
    }
  ((InvariantSet*)is)->diameter = md;
  return md;
}

bool iset_is_subset(const InvariantSet* a, const InvariantSet* b, double tol) {
  if (!a || !b || a->dim != b->dim) return false;
  for (int i = 0; i < a->n_points; i++) {
    int found = 0;
    for (int j = 0; j < b->n_points && !found; j++) {
      double dist = 0.0;
      for (int d = 0; d < a->dim; d++) {
        double diff = a->points[i * a->dim + d] - b->points[j * b->dim + d];
        dist += diff * diff;
      }
      if (sqrt(dist) < tol) found = 1;
    }
    if (!found) return false;
  }
  return true;
}

void iset_print(const InvariantSet* is) {
  if (!is) return;
  printf("InvariantSet: %d pts, d=%.4f, inv=%s\n",
         is->n_points, is->diameter,
         is->is_invariant ? "YES" : "NO");
}

OmegaLimitSet* omega_create(int dim) {
  OmegaLimitSet* ol = calloc(1, sizeof(OmegaLimitSet));
  if (ol) ol->dim = dim;
  return ol;
}
void omega_free(OmegaLimitSet* ol) {
  if (!ol) return; free(ol->omega_points); free(ol->distances_to_limit); free(ol);
}

int omega_compute(OmegaLimitSet* ol, ODEFunc f, const double* x0,
    int dim, double dt, double dur, int ns, double eps) {
  if (!ol || !f || !x0) return -1;
  int steps = (int)(dur / dt);
  int sample_every = steps / (ns > 0 ? ns : 1);
  ol->n_omega = 0;
  ol->cap_omega = ns + 1;
  ol->omega_points = malloc((size_t)(ol->cap_omega * dim) * sizeof(double));
  if (!ol->omega_points) return -1;

  double* x = malloc((size_t)dim * sizeof(double));
  if (!x) return -1;
  memcpy(x, x0, (size_t)dim * sizeof(double));
  for (int i = 0; i < steps; i++) {
    if (i % sample_every == 0 && ol->n_omega < ol->cap_omega) {
      memcpy(ol->omega_points + ol->n_omega * dim, x, (size_t)dim * sizeof(double));
      ol->n_omega++;
    }
    rk4_step(f, x, dim, 0.0, dt);
  }
  free(x);
  ol->approach_rate = 1.0 / dur;
  return ol->n_omega;
}

double omega_convergence_rate(const OmegaLimitSet* ol) {
  return ol ? ol->approach_rate : 0.0;
}

bool omega_is_singleton(const OmegaLimitSet* ol, double tol) {
  if (!ol || ol->n_omega < 2) return true;
  double md = 0.0;
  for (int i = 0; i < ol->n_omega; i++)
    for (int j = i + 1; j < ol->n_omega; j++) {
      double d = 0.0;
      for (int k = 0; k < ol->dim; k++) {
        double df = ol->omega_points[i * ol->dim + k] - ol->omega_points[j * ol->dim + k];
        d += df * df;
      }
      d = sqrt(d);
      if (d > md) md = d;
    }
  return md < tol;
}

bool omega_contains(const OmegaLimitSet* ol, const double* x, double tol) {
  if (!ol || !x) return false;
  for (int i = 0; i < ol->n_omega; i++) {
    double d = 0.0;
    for (int j = 0; j < ol->dim; j++) {
      double df = ol->omega_points[i * ol->dim + j] - x[j];
      d += df * df;
    }
    if (sqrt(d) < tol) return true;
  }
  return false;
}

void omega_print(const OmegaLimitSet* ol) {
  if (!ol) return;
  printf("OmegaLimit: %d pts, rate=%.4f, singleton=%s\n",
         ol->n_omega, ol->approach_rate,
         omega_is_singleton(ol, 1e-4) ? "YES" : "NO");
}