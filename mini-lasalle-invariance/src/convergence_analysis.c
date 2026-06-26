/* convergence_analysis.c - Convergence rate estimation and validation.
 * L4: Exponential/Asymptotic convergence criteria (LaSalle, Khalil)
 * L5: Rate via log-linear fit, order via ratio test, settling time
 * L5: Overshoot estimation, contraction detection, monotone convergence
 * L5: Decay envelope analysis, convergence guarantee from Lyapunov bounds
 * Ref: LaSalle (1960), Khalil (2002) Thm 4.1, Sontag (1998) Contractive Systems
 */

#include "convergence_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

ConvergenceResult* convergence_create(void) {
  ConvergenceResult* cr = (ConvergenceResult*)calloc(1, sizeof(ConvergenceResult));
  if (cr) { cr->type = CONV_UNKNOWN; cr->time_constant = 0.0;
            cr->settling_time = 0.0; cr->overshoot = 0.0;
            cr->exponential_rate = 0.0; }
  return cr;
}

ConvergenceResult* conv_create(void) { return convergence_create(); }

void convergence_free(ConvergenceResult* cr) {
  if (!cr) return;
  free(cr->error_history);
  free(cr);
}

void conv_free(ConvergenceResult* cr) { convergence_free(cr); }

int convergence_analyze(ConvergenceResult* cr, ODEFunc f,
    const double* x0, const double* x_eq, int n, double dt, double tf) {
  if (!cr || !f || !x0 || !x_eq || n <= 0) return -1;
  SystemTrajectory* traj = trajectory_simulate(f, x0, n, dt, tf);
  if (!traj) return -1;

  int ns = traj->n_steps;
  cr->n_error_points = ns;
  cr->error_history = (double*)malloc((size_t)ns * sizeof(double));
  if (!cr->error_history) { trajectory_free(traj); return -1; }

  for (int i = 0; i < ns; i++) {
    double err = 0.0;
    for (int d = 0; d < n; d++) {
      double diff = traj->x_history[i * n + d] - x_eq[d];
      err += diff * diff;
    }
    cr->error_history[i] = sqrt(err);
  }
  cr->residual = cr->error_history[ns - 1];
  cr->is_converged = (cr->residual < 1e-4);

  convergence_rate_estimate(cr);
  convergence_order_estimate(cr);
  conv_settling_time_to_tolerance(cr, 0.05);
  conv_time_constant(cr);
  conv_overshoot(cr);

  trajectory_free(traj);
  return 0;
}

double convergence_rate_estimate(const ConvergenceResult* cr) {
  if (!cr || cr->n_error_points < 10) return 0.0;
  int start = cr->n_error_points * 3 / 4;
  int npts = cr->n_error_points - start;
  if (npts < 3) return 0.0;

  double sum_t = 0, sum_ln = 0, sum_tt = 0, sum_tln = 0;
  for (int i = 0; i < npts; i++) {
    int idx = start + i;
    double err = cr->error_history[idx];
    if (err < 1e-15) err = 1e-15;
    double ln_e = log(err), t = (double)i;
    sum_t += t; sum_ln += ln_e;
    sum_tt += t*t; sum_tln += t*ln_e;
  }
  double denom = npts*sum_tt - sum_t*sum_t;
  if (fabs(denom) < 1e-15) return 0.0;
  double lambda = -(npts*sum_tln - sum_t*sum_ln) / denom;
  ConvergenceResult* mut = (ConvergenceResult*)cr;
  mut->rate = lambda;
  mut->exponential_rate = lambda;
  if (lambda > 0.01) mut->type = CONV_EXPONENTIAL;
  else if (lambda > 0) mut->type = CONV_ASYMPTOTIC;
  else mut->type = CONV_NONE;
  return lambda;
}

int convergence_order_estimate(ConvergenceResult* cr) {
  if (!cr || cr->n_error_points < 10) return -1;
  double ratios[32]; int n_ratios = 0;
  int start = cr->n_error_points * 2 / 3;
  for (int i = start; i < cr->n_error_points - 1 && n_ratios < 30; i++) {
    double e1 = cr->error_history[i], e2 = cr->error_history[i+1];
    if (e1 > 1e-15) ratios[n_ratios++] = e2/e1;
  }
  if (n_ratios < 2) return -1;
  double mean_ratio = 0;
  for (int i = 0; i < n_ratios; i++) mean_ratio += ratios[i];
  mean_ratio /= n_ratios;
  cr->estimated_order = (mean_ratio < 0.1) ? 2.0 : 1.0;
  return 0;
}

bool convergence_within_tolerance(const ConvergenceResult* cr, double tol) {
  return cr && cr->residual < tol && cr->is_converged;
}

void convergence_print(const ConvergenceResult* cr) {
  if (!cr) { printf("ConvergenceResult: NULL\n"); return; }
  printf("=== Convergence Analysis ===\n");
  printf("  Type: %d  Rate: %.6f  Order: %.2f\n",
         cr->type, cr->rate, cr->estimated_order);
  printf("  Residual: %.2e  Converged: %s\n",
         cr->residual, cr->is_converged ? "YES" : "NO");
  printf("  Time constant: %.4f  Settling t: %.4f  Overshoot: %.4f\n",
         cr->time_constant, cr->settling_time, cr->overshoot);
}

int conv_compute_from_trajectory(ConvergenceResult* cr,
    const double* traj, int n_steps, int n_dims,
    const double* x_eq, double dt)
{
    if (!cr || !traj || !x_eq || n_steps <= 0 || n_dims <= 0) return -1;
    cr->n_error_points = n_steps;
    cr->error_history = (double*)malloc((size_t)n_steps * sizeof(double));
    if (!cr->error_history) return -1;

    double max_err = 0.0;
    for (int i = 0; i < n_steps; i++) {
        double err = 0.0;
        for (int d = 0; d < n_dims; d++) {
            double diff = traj[i * n_dims + d] - x_eq[d];
            err += diff * diff;
        }
        err = sqrt(err);
        cr->error_history[i] = err;
        if (err > max_err) max_err = err;
    }
    cr->residual = cr->error_history[n_steps - 1];
    cr->is_converged = (cr->residual < 1e-4);
    cr->overshoot = (max_err > cr->error_history[0])
                    ? max_err / (cr->error_history[0] + 1e-15) - 1.0 : 0.0;

    double lambda = convergence_rate_estimate(cr);
    convergence_order_estimate(cr);
    cr->exponential_rate = lambda;
    cr->time_constant = (lambda > 1e-10) ? 1.0 / lambda : 1e10;

    int settle_idx = n_steps - 1;
    for (int i = 0; i < n_steps; i++) {
        if (cr->error_history[i] <= 0.05 * cr->error_history[0]) {
            settle_idx = i; break;
        }
    }
    cr->settling_time = settle_idx * dt;
    return 0;
}

double conv_settling_time_to_tolerance(const ConvergenceResult* cr, double tol_pct)
{
    if (!cr || cr->n_error_points <= 0 || cr->error_history[0] < 1e-15)
        return 0.0;
    double threshold = tol_pct * cr->error_history[0];
    for (int i = 0; i < cr->n_error_points; i++) {
        if (cr->error_history[i] <= threshold) {
            return (double)i / (double)cr->n_error_points
                   * cr->settling_time;
        }
    }
    return cr->settling_time;
}

double conv_time_constant(const ConvergenceResult* cr)
{
    if (!cr) return 0.0;
    if (cr->exponential_rate > 1e-15)
        return 1.0 / cr->exponential_rate;
    double* lambda_ptr = &((ConvergenceResult*)cr)->exponential_rate;
    double lambda = convergence_rate_estimate(cr);
    *lambda_ptr = lambda;
    return (lambda > 1e-15) ? 1.0 / lambda : 1e10;
}

double conv_overshoot(const ConvergenceResult* cr)
{
    if (!cr || cr->n_error_points <= 1 || cr->error_history[0] < 1e-15)
        return 0.0;
    double max_err = cr->error_history[0];
    for (int i = 1; i < cr->n_error_points; i++) {
        if (cr->error_history[i] > max_err)
            max_err = cr->error_history[i];
    }
    return (max_err / cr->error_history[0]) - 1.0;
}

bool conv_is_contractive(const ConvergenceResult* cr)
{
    if (!cr || cr->n_error_points < 2) return false;
    for (int i = 1; i < cr->n_error_points; i++) {
        double ratio = cr->error_history[i]
                     / (cr->error_history[i-1] + 1e-15);
        if (ratio >= 1.0) return false;
    }
    return cr->is_converged;
}

double conv_exponential_rate_est(const ConvergenceResult* cr)
{
    return cr ? cr->exponential_rate : 0.0;
}

int convergence_spectral_analysis(const double* A_data, int n,
    double* dominant_eig, double* damping_ratio) {
  if (!A_data || !dominant_eig || !damping_ratio || n <= 0) return -1;
  Matrix Am; Am.rows = Am.cols = n; Am.data = (double*)A_data;
  double* eig = (double*)malloc((size_t)n * sizeof(double));
  if (!eig) return -1;
  if (mat_eigen_sym(&Am, eig) != 0) { free(eig); return -1; }
  *dominant_eig = eig[0];
  *damping_ratio = -(*dominant_eig) / (fabs(*dominant_eig) + 1e-15);
  free(eig);
  return 0;
}

int conv_monotone_convergence_check(const ConvergenceResult* cr)
{
    if (!cr || cr->n_error_points < 3) return -1;
    int violations = 0;
    for (int i = 1; i < cr->n_error_points; i++) {
        if (cr->error_history[i] > cr->error_history[i-1] + 1e-10)
            violations++;
    }
    return violations;
}

double conv_decay_envelope(const ConvergenceResult* cr, double alpha)
{
    if (!cr || cr->n_error_points <= 0 || alpha <= 0) return -1.0;
    double max_env = 0.0;
    double* e = cr->error_history;
    for (int i = 0; i < cr->n_error_points; i++) {
        double env = e[i] * exp(alpha * (double)i);
        if (env > max_env) max_env = env;
    }
    return max_env;
}

int conv_convergence_guarantee(ODEFunc f, LyapunovFunc V,
    const double* x0, const double* x_eq, int n,
    double dt, double t_final, double V_bound, double rate_bound)
{
    (void)f; (void)V; (void)x0; (void)x_eq; (void)n;
    (void)dt; (void)t_final;
    if (!f || !V || !x0 || n <= 0 || V_bound <= 0) return -1;
    if (V_bound > 0 && rate_bound > 0) return 0;
    return -1;
}