#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fts_core.h"

/* ============================================================
 * fts_core.c -- Finite-Time Stability Core Implementation
 *
 * Simulation framework: system lifecycle, RK4/Euler/Heun
 * integration, trajectory recording, vector operations,
 * factory functions for standard dynamical systems.
 *
 * Key formulas:
 *   RK4: x_{n+1} = x_n + h/6*(k1 + 2k2 + 2k3 + k4)
 *   Euler: x_{n+1} = x_n + h*f(x_n)
 *   Heun: x_{n+1} = x_n + h/2*(f(x_n) + f(x_n + h*f(x_n)))
 *
 * Refs: Bhat & Bernstein (2000), Polyakov (2012)
 * ============================================================ */

/* ---------- Lifecycle ---------- */

FTSSystem* fts_create(int dim, int np, double dt) {
    if (dim <= 0 || dim > FTS_MAX_DIM) return NULL;
    FTSSystem* s = calloc(1, sizeof(FTSSystem));
    if (!s) return NULL;
    s->name = strdup("FTS");
    if (!s->name) { free(s); return NULL; }
    s->dim = dim; s->np = np;
    s->dt = (dt > 0) ? dt : FTS_DEFAULT_DT;
    s->t = 0.0;
    if (np > 0) {
        s->params = calloc(np, sizeof(double));
        if (!s->params) { free(s->name); free(s); return NULL; }
    }
    s->traj.cap = 1000; s->traj.n = 0;
    s->traj.pts = calloc(s->traj.cap, sizeof(FTSState));
    if (!s->traj.pts) {
        free(s->params); free(s->name); free(s); return NULL;
    }
    return s;
}

void fts_free(FTSSystem* sys) {
    if (!sys) return;
    free(sys->name); free(sys->params);
    free(sys->traj.pts); free(sys);
}

/* ---------- State Access ---------- */

void fts_set_state(FTSSystem* sys, double* x) {
    if (!sys || !x) return;
    memcpy(sys->state.x, x, sys->dim * sizeof(double));
}

double fts_get_state(FTSSystem* sys, int idx) {
    if (!sys || idx < 0 || idx >= sys->dim) return 0.0;
    return sys->state.x[idx];
}

void fts_copy_state(FTSState* dst, const FTSState* src) {
    if (!dst || !src) return;
    dst->dim = src->dim;
    memcpy(dst->x, src->x, src->dim * sizeof(double));
}

void fts_set_param(FTSSystem* sys, int idx, double v) {
    if (!sys || idx < 0 || idx >= sys->np) return;
    sys->params[idx] = v;
}

double fts_get_param(FTSSystem* sys, int idx) {
    if (!sys || idx < 0 || idx >= sys->np) return 0.0;
    return sys->params[idx];
}

/* ---------- Integration Methods ---------- */

void fts_step_euler(FTSSystem* sys) {
    /* Forward Euler: x_{n+1} = x_n + h*f(x_n)
     * First-order method. Local error O(h^2), global O(h).
     * Simple but requires small step sizes for accuracy. */
    if (!sys || !sys->rhs) return;
    int d = sys->dim; FTSState dx; dx.dim = d;
    sys->rhs(sys, sys->state, &dx);
    for (int i = 0; i < d; i++)
        sys->state.x[i] += sys->dt * dx.x[i];
    sys->t += sys->dt;
}

void fts_step_heun(FTSSystem* sys) {
    /* Heun's method (explicit trapezoidal, RK2):
     * k1 = f(x_n)
     * k2 = f(x_n + h*k1)
     * x_{n+1} = x_n + h/2*(k1 + k2)
     * Second-order method. Local error O(h^3), global O(h^2). */
    if (!sys || !sys->rhs) return;
    int d = sys->dim; FTSState x0 = sys->state;
    FTSState k1, k2, xtmp; k1.dim = k2.dim = xtmp.dim = d;
    sys->rhs(sys, x0, &k1);
    for (int i = 0; i < d; i++)
        xtmp.x[i] = x0.x[i] + sys->dt * k1.x[i];
    sys->rhs(sys, xtmp, &k2);
    for (int i = 0; i < d; i++)
        sys->state.x[i] = x0.x[i] + 0.5 * sys->dt * (k1.x[i] + k2.x[i]);
    sys->t += sys->dt;
}

void fts_step_rk4(FTSSystem* sys) {
    /* Classical 4th-order Runge-Kutta (RK4):
     * k1 = f(x_n)
     * k2 = f(x_n + h/2 * k1)
     * k3 = f(x_n + h/2 * k2)
     * k4 = f(x_n + h * k3)
     * x_{n+1} = x_n + h/6 * (k1 + 2k2 + 2k3 + k4)
     *
     * Local truncation error O(h^5), global error O(h^4).
     * Gold standard for non-stiff ODE simulation.
     * Four function evaluations per step. */
    if (!sys || !sys->rhs) return;
    int d = sys->dim; FTSState x0 = sys->state;
    FTSState k1, k2, k3, k4;
    FTSSystem ts = *sys; /* shallow copy for temporary state */
    k1.dim = k2.dim = k3.dim = k4.dim = d;

    sys->rhs(sys, x0, &k1);
    for (int i = 0; i < d; i++)
        ts.state.x[i] = x0.x[i] + 0.5 * sys->dt * k1.x[i];
    sys->rhs(&ts, ts.state, &k2);
    for (int i = 0; i < d; i++)
        ts.state.x[i] = x0.x[i] + 0.5 * sys->dt * k2.x[i];
    sys->rhs(&ts, ts.state, &k3);
    for (int i = 0; i < d; i++)
        ts.state.x[i] = x0.x[i] + sys->dt * k3.x[i];
    sys->rhs(&ts, ts.state, &k4);

    for (int i = 0; i < d; i++) {
        sys->state.x[i] = x0.x[i] + (sys->dt / 6.0) *
            (k1.x[i] + 2.0*k2.x[i] + 2.0*k3.x[i] + k4.x[i]);
    }
    sys->t += sys->dt;
}

void fts_integrate(FTSSystem* sys, double dur, int rec) {
    /* Integrate system for 'dur' time units using RK4.
     * If rec > 0, record state every 'rec' steps to trajectory. */
    if (!sys) return;
    int steps = (int)(dur / sys->dt);
    if (steps <= 0) steps = 1;
    for (int i = 0; i < steps; i++) {
        fts_step_rk4(sys);
        if (rec > 0 && (i % rec == 0)) {
            fts_record(sys);
        }
    }
}

/* ---------- Trajectory Management ---------- */

void fts_record(FTSSystem* sys) {
    if (!sys) return;
    if (sys->traj.n >= sys->traj.cap) {
        sys->traj.cap *= 2;
        if (sys->traj.cap > FTS_MAX_HIST)
            sys->traj.cap = FTS_MAX_HIST;
        FTSState* new_pts = realloc(sys->traj.pts,
            sys->traj.cap * sizeof(FTSState));
        if (!new_pts) return;
        sys->traj.pts = new_pts;
    }
    if (sys->traj.n < sys->traj.cap) {
        sys->traj.pts[sys->traj.n] = sys->state;
        sys->traj.n++;
    }
}

void fts_clear_trajectory(FTSSystem* sys) {
    if (!sys) return;
    sys->traj.n = 0;
}

int fts_trajectory_count(FTSSystem* sys) {
    return sys ? sys->traj.n : 0;
}

/* ---------- Norms and Distance ---------- */

double fts_norm(FTSState s) {
    /* Euclidean L2 norm: ||x||_2 = sqrt(sum(x_i^2)) */
    double sum = 0.0;
    int n = (s.dim < FTS_MAX_DIM) ? s.dim : FTS_MAX_DIM;
    for (int i = 0; i < n; i++) sum += s.x[i] * s.x[i];
    return sqrt(sum);
}

double fts_norm_type(FTSState s, FTSNormType type) {
    int n = (s.dim < FTS_MAX_DIM) ? s.dim : FTS_MAX_DIM;
    switch (type) {
    case FTS_NORM_L2: {
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += s.x[i] * s.x[i];
        return sqrt(sum);
    }
    case FTS_NORM_L1: {
        double sum = 0.0;
        for (int i = 0; i < n; i++) sum += fabs(s.x[i]);
        return sum;
    }
    case FTS_NORM_LINF: {
        double mx = 0.0;
        for (int i = 0; i < n; i++) {
            double a = fabs(s.x[i]);
            if (a > mx) mx = a;
        }
        return mx;
    }
    default: return fts_norm(s);
    }
}

double fts_squared_norm(FTSState s) {
    double sum = 0.0;
    int n = (s.dim < FTS_MAX_DIM) ? s.dim : FTS_MAX_DIM;
    for (int i = 0; i < n; i++) sum += s.x[i] * s.x[i];
    return sum;
}

double fts_dot_product(FTSState a, FTSState b) {
    int n = (a.dim < b.dim) ? a.dim : b.dim;
    if (n > FTS_MAX_DIM) n = FTS_MAX_DIM;
    double sum = 0.0;
    for (int i = 0; i < n; i++) sum += a.x[i] * b.x[i];
    return sum;
}

double fts_distance(FTSState a, FTSState b, int dim) {
    if (dim <= 0) return 0.0;
    if (dim > FTS_MAX_DIM) dim = FTS_MAX_DIM;
    double sum = 0.0;
    for (int i = 0; i < dim; i++) {
        double d = a.x[i] - b.x[i]; sum += d * d;
    }
    return sqrt(sum);
}

/* ---------- Vector Operations ---------- */

void fts_vector_add(FTSState* dst, FTSState a, FTSState b, double s) {
    if (!dst) return;
    int n = (a.dim < b.dim) ? a.dim : b.dim;
    if (n > FTS_MAX_DIM) n = FTS_MAX_DIM;
    for (int i = 0; i < n; i++)
        dst->x[i] = a.x[i] + s * b.x[i];
    dst->dim = n;
}

void fts_vector_scale(FTSState* dst, FTSState src, double s) {
    if (!dst) return;
    int n = (src.dim < FTS_MAX_DIM) ? src.dim : FTS_MAX_DIM;
    for (int i = 0; i < n; i++)
        dst->x[i] = s * src.x[i];
    dst->dim = n;
}

void fts_vector_sub(FTSState* dst, FTSState a, FTSState b) {
    if (!dst) return;
    int n = (a.dim < b.dim) ? a.dim : b.dim;
    if (n > FTS_MAX_DIM) n = FTS_MAX_DIM;
    for (int i = 0; i < n; i++)
        dst->x[i] = a.x[i] - b.x[i];
    dst->dim = n;
}

void fts_vector_copy(FTSState* dst, FTSState src) {
    if (!dst) return;
    dst->dim = src.dim;
    memcpy(dst->x, src.x, src.dim * sizeof(double));
}

/* ---------- System Analysis ---------- */

void fts_print(FTSSystem* sys) {
    if (!sys) return;
    printf("t=%.6f ", sys->t);
    for (int i = 0; i < sys->dim; i++)
        printf("x%d=%.6f ", i, sys->state.x[i]);
    printf("\n");
}

int fts_dimension(FTSSystem* sys) {
    return sys ? sys->dim : 0;
}

double fts_current_time(FTSSystem* sys) {
    return sys ? sys->t : 0.0;
}

bool fts_is_valid(FTSSystem* sys) {
    if (!sys) return false;
    if (sys->dim <= 0 || sys->dim > FTS_MAX_DIM) return false;
    if (sys->dt <= 0) return false;
    if (!sys->rhs) return false;
    return true;
}

/* ---------- Factory Functions ---------- */

static void di_rhs(FTSSystem* sys, FTSState x, FTSState* dx) {
    /* Double integrator: dx1/dt = x2, dx2/dt = u
     * params[0] = control input u */
    (void)sys;
    dx->x[0] = x.x[1];
    dx->x[1] = sys->params[0];
}

FTSSystem* fts_create_double_integrator(double dt) {
    FTSSystem* s = fts_create(2, 1, dt);
    if (!s) return NULL;
    s->params[0] = 0.0; s->rhs = di_rhs;
    s->state.x[0] = 1.0; s->state.x[1] = 0.0;
    return s;
}

static void sfts_rhs(FTSSystem* sys, FTSState x, FTSState* dx) {
    /* Scalar FTS: dx/dt = -k * |x|^alpha * sign(x)
     * params[0] = k (gain), params[1] = alpha (exponent)
     *
     * For 0 < alpha < 1: finite-time stable.
     * Settling time: T = |x0|^(1-alpha) / (k*(1-alpha))
     *
     * For alpha = 1: exponential convergence.
     * For alpha > 1: faster-than-exponential (not FTS in strict sense). */
    double sg = (x.x[0] > 0) ? 1.0 : ((x.x[0] < 0) ? -1.0 : 0.0);
    dx->x[0] = -sys->params[0] * pow(fabs(x.x[0]), sys->params[1]) * sg;
}

FTSSystem* fts_create_scalar_fts(double alpha, double dt) {
    FTSSystem* s = fts_create(1, 2, dt);
    if (!s) return NULL;
    s->params[0] = 1.0; s->params[1] = alpha;
    s->rhs = sfts_rhs; s->state.x[0] = 1.0;
    return s;
}

static void ho_rhs(FTSSystem* sys, FTSState x, FTSState* dx) {
    /* Harmonic oscillator: dx1/dt = x2, dx2/dt = -omega^2 * x1
     * params[0] = omega.
     * Energy: E = 0.5*(x2^2 + omega^2*x1^2), conserved. */
    double w2 = sys->params[0] * sys->params[0];
    dx->x[0] = x.x[1]; dx->x[1] = -w2 * x.x[0];
}

FTSSystem* fts_create_harmonic_oscillator(double omega, double dt) {
    FTSSystem* s = fts_create(2, 1, dt);
    if (!s) return NULL;
    s->params[0] = omega; s->rhs = ho_rhs;
    s->state.x[0] = 1.0; s->state.x[1] = 0.0;
    return s;
}

static void vdp_rhs(FTSSystem* sys, FTSState x, FTSState* dx) {
    /* Van der Pol oscillator:
     * dx1/dt = x2
     * dx2/dt = mu*(1-x1^2)*x2 - x1
     * params[0] = mu. Limit cycle for mu > 0. */
    double mu = sys->params[0];
    dx->x[0] = x.x[1];
    dx->x[1] = mu * (1.0 - x.x[0]*x.x[0]) * x.x[1] - x.x[0];
}

FTSSystem* fts_create_van_der_pol(double mu, double dt) {
    FTSSystem* s = fts_create(2, 1, dt);
    if (!s) return NULL;
    s->params[0] = mu; s->rhs = vdp_rhs;
    s->state.x[0] = 2.0; s->state.x[1] = 0.0;
    return s;
}

static void lorenz_rhs(FTSSystem* sys, FTSState x, FTSState* dx) {
    /* Lorenz system:
     * dx/dt = sigma*(y-x)
     * dy/dt = x*(rho-z) - y
     * dz/dt = x*y - beta*z
     * params[0]=sigma, params[1]=rho, params[2]=beta.
     * Classic parameters: sigma=10, rho=28, beta=8/3. */
    double s = sys->params[0], r = sys->params[1], b = sys->params[2];
    dx->x[0] = s * (x.x[1] - x.x[0]);
    dx->x[1] = x.x[0] * (r - x.x[2]) - x.x[1];
    dx->x[2] = x.x[0] * x.x[1] - b * x.x[2];
}

FTSSystem* fts_create_lorenz(double sigma, double rho, double beta,
                              double dt) {
    FTSSystem* s = fts_create(3, 3, dt);
    if (!s) return NULL;
    s->params[0] = sigma; s->params[1] = rho; s->params[2] = beta;
    s->rhs = lorenz_rhs;
    s->state.x[0] = 1.0; s->state.x[1] = 1.0; s->state.x[2] = 1.0;
    return s;
}

static void brockett_rhs(FTSSystem* sys, FTSState x, FTSState* dx) {
    /* Brockett integrator (nonholonomic system):
     * dx/dt = u1, dy/dt = u2, dz/dt = x*u2 - y*u1
     * params[0]=u1, params[1]=u2.
     * Classic example in nonholonomic control. */
    double u1 = sys->params[0], u2 = sys->params[1];
    dx->x[0] = u1; dx->x[1] = u2;
    dx->x[2] = x.x[0] * u2 - x.x[1] * u1;
}

FTSSystem* fts_create_brockett_integrator(double dt) {
    FTSSystem* s = fts_create(3, 2, dt);
    if (!s) return NULL;
    s->params[0] = 0.0; s->params[1] = 0.0;
    s->rhs = brockett_rhs;
    s->state.x[0] = 0.0; s->state.x[1] = 0.0; s->state.x[2] = 0.0;
    return s;
}

/* ---------- Energy Functions ---------- */

double fts_energy_quadratic(FTSState x, double* P, int dim) {
    /* V(x) = x^T * P * x, for symmetric positive definite P.
     * P stored row-major as dim*dim array.
     * Used for quadratic Lyapunov functions. */
    if (!P || dim <= 0) return 0.0;
    double energy = 0.0;
    for (int i = 0; i < dim; i++) {
        double row_sum = 0.0;
        for (int j = 0; j < dim; j++)
            row_sum += P[i*dim + j] * x.x[j];
        energy += x.x[i] * row_sum;
    }
    return energy;
}

double fts_lyap_default(FTSState x, void* params) {
    /* Default Lyapunov function: V(x) = 1/2 * ||x||^2
     * V_dot = x^T * f(x). If V_dot < 0, origin is stable.
     * params unused but kept for interface compatibility. */
    (void)params;
    return 0.5 * fts_squared_norm(x);
}

/* ---------- Finite Escape Detection ---------- */

bool fts_has_finite_escape(FTSSystem* sys, double Tmax) {
    /* Check for finite escape time:
     * ||x(t)|| -> infinity as t -> T_escape < infinity.
     * Practical detection: ||x|| exceeds 1e6 within Tmax. */
    if (!sys) return false;
    FTSSystem* cp = fts_create(sys->dim, sys->np, sys->dt);
    if (!cp) return false;
    fts_copy_state(&cp->state, &sys->state);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));

    bool escaped = false;
    const double thresh = 1e6;
    while (cp->t < Tmax) {
        fts_step_rk4(cp);
        double nrm = fts_norm(cp->state);
        if (nrm > thresh || isnan(nrm) || isinf(nrm)) {
            escaped = true; break;
        }
    }
    fts_free(cp);
    return escaped;
}

double fts_divergence(FTSSystem* sys, FTSState x) {
    /* Compute divergence of vector field f at point x:
     * div f = sum_i df_i/dx_i
     * Approximated by central finite differences with h = 1e-6.
     * div f < 0 indicates volume contraction. */
    if (!sys || !sys->rhs) return 0.0;
    int d = sys->dim;
    double h = 1e-6, div = 0.0;
    FTSState xp = x, xm = x, fp, fm;
    fp.dim = fm.dim = d;
    for (int i = 0; i < d; i++) {
        xp.x[i] = x.x[i] + h; xm.x[i] = x.x[i] - h;
        sys->rhs(sys, xp, &fp); sys->rhs(sys, xm, &fm);
        div += (fp.x[i] - fm.x[i]) / (2.0 * h);
        xp.x[i] = xm.x[i] = x.x[i];
    }
    return div;
}

int fts_escape_time(FTSSystem* sys, double threshold, double Tmax,
                     double* t_escape) {
    /* Find time when ||x(t)|| first exceeds threshold.
     * Returns: 0 = no escape, 1 = escaped (t_escape set), -1 = error. */
    if (!sys || !t_escape) return -1;
    FTSSystem* cp = fts_create(sys->dim, sys->np, sys->dt);
    if (!cp) return -1;
    fts_copy_state(&cp->state, &sys->state);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));
    int result = 0;
    while (cp->t < Tmax) {
        fts_step_rk4(cp);
        if (fts_norm(cp->state) > threshold) {
            *t_escape = cp->t; result = 1; break;
        }
    }
    fts_free(cp);
    return result;
}

/* ---------- Utility Functions ---------- */

double fts_clamp(double x, double lo, double hi) {
    if (x < lo) return lo; if (x > hi) return hi; return x;
}

double fts_sign(double x) {
    if (x > 0.0) return 1.0; if (x < 0.0) return -1.0; return 0.0;
}

double fts_saturate(double x, double limit) {
    return fts_clamp(x, -limit, limit);
}
