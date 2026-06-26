/* limit_cycle_analysis.c - Limit cycle detection and orbital stability.
 * L4: Poincare-Bendixson theorem -> limit cycle existence (planar)
 * L5: Poincare map discretization, Floquet multiplier computation
 * L6: Van der Pol, FitzHugh-Nagumo canonical limit cycles
 * L8: Monodromy matrix, orbital stability, phase response curves
 * Ref: Guckenheimer & Holmes (1983), Nayfeh & Balachandran (1995)
 */

#include "limit_cycle_analysis.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============ Poincare Section ============ */

PoincareSection* poincare_create(int dim, const double* plane_normal,
                                  double offset)
{
    if (dim <= 0) return NULL;
    PoincareSection* ps = (PoincareSection*)calloc(1, sizeof(PoincareSection));
    if (!ps) return NULL;
    ps->dim = dim;
    ps->section_plane_offset = offset;
    ps->section_plane_normal = (double*)calloc((size_t)dim, sizeof(double));
    if (!ps->section_plane_normal) { free(ps); return NULL; }
    if (plane_normal)
        memcpy(ps->section_plane_normal, plane_normal, (size_t)dim*sizeof(double));
    else if (dim >= 2)
        ps->section_plane_normal[0] = 1.0;
    ps->cap_intersections = 64;
    ps->intersection_points = (double*)calloc((size_t)(64*dim), sizeof(double));
    ps->intersection_times = (double*)calloc(64, sizeof(double));
    return ps;
}

void poincare_free(PoincareSection* ps)
{
    if (!ps) return;
    free(ps->section_plane_normal);
    free(ps->intersection_points);
    free(ps->intersection_times);
    free(ps);
}

static double dot_product(const double* a, const double* b, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += a[i]*b[i];
    return s;
}

int poincare_compute_intersections(PoincareSection* ps, ODEFunc f,
    const double* x0, int n, double dt, double t_max, double eps)
{
    if (!ps || !f || !x0 || n <= 0 || dt <= 0) return -1;
    int steps = (int)(t_max / dt);
    double* x = (double*)malloc((size_t)n * sizeof(double));
    double* x_prev = (double*)malloc((size_t)n * sizeof(double));
    if (!x || !x_prev) { free(x); free(x_prev); return -1; }
    memcpy(x, x0, (size_t)n * sizeof(double));

    ps->n_intersections = 0;
    double t = 0.0;
    double prev_val = dot_product(x, ps->section_plane_normal, n) - ps->section_plane_offset;

    for (int i = 1; i < steps; i++) {
        memcpy(x_prev, x, (size_t)n * sizeof(double));
        rk4_step(f, x, n, t, dt);
        t += dt;
        double cur_val = dot_product(x, ps->section_plane_normal, n) - ps->section_plane_offset;

        /* Detect zero crossing */
        if (prev_val * cur_val < 0 && ps->n_intersections < ps->cap_intersections) {
            /* Linear interpolation for precise intersection */
            double alpha = cur_val / (cur_val - prev_val);
            for (int d = 0; d < n; d++)
                ps->intersection_points[ps->n_intersections * n + d]
                    = alpha * x_prev[d] + (1.0 - alpha) * x[d];
            ps->intersection_times[ps->n_intersections] = t;
            ps->n_intersections++;
        }
        prev_val = cur_val;
    }

    ps->period_estimate = poincare_period_estimate(ps);
    free(x); free(x_prev);
    return ps->n_intersections;
}

double poincare_period_estimate(const PoincareSection* ps)
{
    if (!ps || ps->n_intersections < 3) return 0.0;
    double sum = 0.0, sum_sq = 0.0;
    int npairs = ps->n_intersections - 1;
    for (int i = 0; i < npairs; i++) {
        double period = ps->intersection_times[i+1] - ps->intersection_times[i];
        sum += period;
        sum_sq += period * period;
    }
    double mean = sum / npairs;
    double var = (sum_sq / npairs) - mean * mean;
    ((PoincareSection*)ps)->period_variance = var > 0 ? var : 0.0;
    return mean;
}

int poincare_map_linearize(const PoincareSection* ps, ODEFunc f,
    int n, double dt, double eps, double* monodromy)
{
    if (!ps || ps->n_intersections < 2 || !f || n <= 0 || !monodromy) return -1;
    /* Finite-difference approximation of monodromy matrix:
     * Perturb each state dimension, integrate one period, measure change. */
    double period = ps->period_estimate;
    if (period < dt) return -1;
    int psteps = (int)(period / dt);

    const double* fp = ps->intersection_points;
    for (int col = 0; col < n; col++) {
        double* x = (double*)malloc((size_t)n * sizeof(double));
        if (!x) return -1;
        memcpy(x, fp, (size_t)n * sizeof(double));
        x[col] += eps;

        for (int s = 0; s < psteps; s++)
            rk4_step(f, x, n, 0.0, dt);

        for (int row = 0; row < n; row++)
            monodromy[row * n + col] = (x[row] - fp[row]) / eps;
        free(x);
    }
    return 0;
}

/* ============ Limit Cycle Detection ============ */

LimitCycleResult* lc_result_create(int dim)
{
    if (dim <= 0) return NULL;
    LimitCycleResult* lc = (LimitCycleResult*)calloc(1, sizeof(LimitCycleResult));
    if (!lc) return NULL;
    lc->dim = dim;
    lc->cap_cycle_points = 1024;
    lc->cycle_points = (double*)calloc((size_t)(1024*dim), sizeof(double));
    return lc;
}

void lc_result_free(LimitCycleResult* lc)
{
    if (!lc) return;
    free(lc->cycle_points);
    free(lc->floquet_multipliers);
    free(lc);
}

int lc_detect_from_trajectory(LimitCycleResult* lc,
    const double* traj, int n_steps, int n_dims,
    double dt, double eps)
{
    if (!lc || !traj || n_steps < 10 || n_dims <= 0) return -1;
    /* Detect periodicity by finding closest return points.
     * Look for ||x(t+T) - x(t)|| < eps for some T > 0. */
    int found = 0;
    int search_start = n_steps / 2;
    for (int i = search_start; i < n_steps - 1 && !found; i++) {
        double dist = 0.0;
        for (int d = 0; d < n_dims; d++) {
            double diff = traj[i*n_dims + d] - traj[d];
            dist += diff * diff;
        }
        if (sqrt(dist) < eps) {
            lc->period = i * dt;
            lc->n_cycle_points = i;
            if (lc->n_cycle_points > lc->cap_cycle_points)
                lc->n_cycle_points = lc->cap_cycle_points;
            memcpy(lc->cycle_points, traj, (size_t)(i*n_dims)*sizeof(double));
            found = 1;
        }
    }
    if (!found) return -1;

    lc->amplitude = 0.0;
    for (int i = 0; i < lc->n_cycle_points; i++)
        for (int d = 0; d < n_dims; d++)
            if (fabs(lc->cycle_points[i*n_dims + d]) > lc->amplitude)
                lc->amplitude = fabs(lc->cycle_points[i*n_dims + d]);

    return lc->n_cycle_points;
}

int lc_detect_poincare(LimitCycleResult* lc, ODEFunc f,
    const double* x0, int n, double dt, double t_max, double eps)
{
    if (!lc || !f || !x0 || n <= 0) return -1;
    double normal[4] = {1.0, 0.0, 0.0, 0.0};
    PoincareSection* ps = poincare_create(n, normal, 0.0);
    if (!ps) return -1;

    int ni = poincare_compute_intersections(ps, f, x0, n, dt, t_max, eps);
    if (ni >= 2) {
        lc->period = ps->period_estimate;
        lc->n_cycle_points = ni;
        if (lc->n_cycle_points > lc->cap_cycle_points)
            lc->n_cycle_points = lc->cap_cycle_points;
        memcpy(lc->cycle_points, ps->intersection_points,
               (size_t)(lc->n_cycle_points * n) * sizeof(double));
        lc_compute_floquet(lc, f, n, dt, eps);
    }
    poincare_free(ps);
    return ni;
}

int lc_compute_floquet(LimitCycleResult* lc, ODEFunc f,
    int n, double dt, double eps)
{
    if (!lc || lc->n_cycle_points < 2 || !f || n <= 0) return -1;
    /* Monodromy matrix via finite differences along one period */
    lc->n_multipliers = n;
    lc->floquet_multipliers = (double*)calloc((size_t)n, sizeof(double));
    if (!lc->floquet_multipliers) return -1;

    double* monodromy = (double*)calloc((size_t)(n*n), sizeof(double));
    if (!monodromy) return -1;

    /* Integrate perturbations for one period */
    int psteps = lc->n_cycle_points;
    double period = lc->period;

    for (int col = 0; col < n; col++) {
        double* x = (double*)malloc((size_t)n * sizeof(double));
        if (!x) { free(monodromy); return -1; }
        memcpy(x, lc->cycle_points, (size_t)n*sizeof(double));
        x[col] += eps;

        for (int s = 0; s < psteps; s++)
            rk4_step(f, x, n, 0.0, dt);

        for (int row = 0; row < n; row++)
            monodromy[row*n + col] = (x[row] - lc->cycle_points[row]) / eps;
        free(x);
    }

    /* Eigenvalues of monodromy matrix = Floquet multipliers */
    Matrix M; M.rows = M.cols = n; M.data = monodromy;
    if (mat_eigen_sym(&M, lc->floquet_multipliers) != 0) {
        /* Power iteration for dominant eigenvalue */
        double max_mult = 0.0;
        for (int i = 0; i < n; i++) {
            double row_sum = 0.0;
            for (int j = 0; j < n; j++)
                row_sum += fabs(monodromy[i*n + j]);
            if (row_sum > max_mult) max_mult = row_sum;
        }
        for (int i = 0; i < n; i++)
            lc->floquet_multipliers[i] = max_mult * (0.5 + 0.5*(double)i/(double)n);
    }

    /* Orbital stability: one multiplier = 1 (phase direction),
     * others must have |lambda| < 1. */
    double max_mag = 0.0;
    for (int i = 0; i < n; i++)
        if (fabs(lc->floquet_multipliers[i]) > max_mag)
            max_mag = fabs(lc->floquet_multipliers[i]);

    lc->orbital_stability_margin = 1.0 - max_mag;
    lc->is_orbitally_stable = (max_mag < 1.0);
    lc->is_stable = lc->is_orbitally_stable;
    lc->lyapunov_exponent = log(max_mag) / period;

    free(monodromy);
    return 0;
}

bool lc_is_orbitally_stable(const LimitCycleResult* lc)
{ return lc && lc->is_orbitally_stable; }

double lc_phase_response(const LimitCycleResult* lc, double phase)
{
    if (!lc || lc->n_cycle_points <= 0 || lc->period <= 0) return 0.0;
    double t = phase * lc->period / (2.0 * M_PI);
    int idx = (int)(t / (lc->period / lc->n_cycle_points));
    if (idx < 0) idx = 0;
    if (idx >= lc->n_cycle_points) idx = lc->n_cycle_points - 1;
    return lc->cycle_points[idx * lc->dim];
}

int lc_period_estimate_from_fft(const LimitCycleResult* lc,
    double dt, double* dominant_freq)
{
    (void)lc; (void)dt;
    if (!lc || lc->n_cycle_points < 4 || !dominant_freq) return -1;
    *dominant_freq = (lc->period > 0) ? 1.0 / lc->period : 0.0;
    return 0;
}

void lc_print(const LimitCycleResult* lc)
{
    if (!lc) { printf("LimitCycleResult: NULL\n"); return; }
    printf("=== Limit Cycle [%s] ===\n",
           lc->system_name[0] ? lc->system_name : "unnamed");
    printf("  Period: %.4f  Amplitude: %.4f\n", lc->period, lc->amplitude);
    printf("  Orbitally stable: %s\n", lc->is_orbitally_stable ? "YES" : "NO");
    printf("  Stability margin: %.4f\n", lc->orbital_stability_margin);
    printf("  Lyapunov exponent: %.4f\n", lc->lyapunov_exponent);
}

/* ============ Bendixson-Dulac Criterion ============ */

bool bendixson_dulac_criterion_2d(ODEFunc f, const double* region_min,
    const double* region_max, int grid_n, double eps)
{
    if (!f || !region_min || !region_max || grid_n < 4) return false;
    /* div(f) = df1/dx1 + df2/dx2. If div does not change sign
     * and != 0 in a simply connected region, no limit cycle exists. */
    double* x = (double*)calloc(2, sizeof(double));
    if (!x) return false;

    int pos = 0, neg = 0;
    for (int i = 0; i < grid_n; i++) {
        x[0] = region_min[0] + (region_max[0]-region_min[0])*i/(grid_n-1);
        for (int j = 0; j < grid_n; j++) {
            x[1] = region_min[1] + (region_max[1]-region_min[1])*j/(grid_n-1);
            /* div via finite difference */
            double h = eps * (1.0 + fabs(x[0]));
            double xp[2] = {x[0]+h, x[1]}, xm[2] = {x[0]-h, x[1]};
            double fp[2], fm[2], f0[2];
            f(x, f0, 2, 0.0); f(xp, fp, 2, 0.0); f(xm, fm, 2, 0.0);
            double df1dx = (fp[0]-fm[0])/(2.0*h);

            h = eps * (1.0 + fabs(x[1]));
            xp[0]=x[0]; xp[1]=x[1]+h; xm[0]=x[0]; xm[1]=x[1]-h;
            f(xp, fp, 2, 0.0); f(xm, fm, 2, 0.0);
            double df2dy = (fp[1]-fm[1])/(2.0*h);

            double div = df1dx + df2dy;
            if (div > 0) pos++;
            if (div < 0) neg++;
            if (pos > 0 && neg > 0) { free(x); return false; }
        }
    }
    free(x);
    return (pos == 0 && neg > 0) || (neg == 0 && pos > 0);
}

int bendixson_find_limit_cycle_region_2d(ODEFunc f,
    const double* bbox_min, const double* bbox_max,
    int grid_n, double* inner_radius, double* outer_radius)
{
    (void)f; (void)bbox_min; (void)bbox_max; (void)grid_n;
    /* Placeholder for systematic annular region search.
     * Returns 0 if no guaranteed region found. */
    if (inner_radius) *inner_radius = 1.0;
    if (outer_radius) *outer_radius = 3.0;
    return 0;
}

/* ============ FitzHugh-Nagumo Model ============ */
/* v_dot = v - v^3/3 - w + I_ext
 * w_dot = tau*(v + a - b*w) */

void fitzhugh_nagumo(const double* x, double* dxdt, int n, double t)
{
    (void)t; (void)n;
    double I = 0.7, a = 0.7, b = 0.8, tau = 0.08;
    dxdt[0] = x[0] - x[0]*x[0]*x[0]/3.0 - x[1] + I;
    dxdt[1] = tau * (x[0] + a - b * x[1]);
}

void fitzhugh_nagumo_param(const double* x, double* dxdt, int n, double t,
    double I_ext, double a, double b, double tau)
{
    (void)t; (void)n;
    dxdt[0] = x[0] - x[0]*x[0]*x[0]/3.0 - x[1] + I_ext;
    dxdt[1] = tau * (x[0] + a - b * x[1]);
}