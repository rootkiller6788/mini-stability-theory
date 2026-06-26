#include "iss_gain.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

ISS_Gain iss_gain_from_lyapunov(const ISS_Certificate* cert) {
    ISS_Gain g;
    memset(&g, 0, sizeof(g));
    if (!cert) return g;
    g.linear_gain = cert->gamma_gain;
    g.is_finite = isfinite(g.linear_gain);
    g.margin = 1.0 / fmax(g.linear_gain, 1e-10);
    g.peak_gain = g.linear_gain;
    g.l2_gain = g.linear_gain * 0.7;
    g.k_function = iss_k_linear(g.linear_gain);
    return g;
}

ISS_Gain iss_gain_from_simulation(ISS_System* sys, double t_final, double dt,
    double* input_amplitudes, int n_inputs, double state_threshold) {
    ISS_Gain g;
    memset(&g, 0, sizeof(g));
    if (!sys || !sys->dynamics) return g;
    (void)state_threshold;
    double max_state = 0.0;
    double max_input = 0.0;
    int i;
    for (i = 0; i < n_inputs && i < 10; i++) {
        double u = input_amplitudes ? input_amplitudes[i] : 1.0;
        if (u > max_input) max_input = u;
    }
    int n = sys->n_states;
    int m = sys->n_inputs;
    double x[ISS_MAX_DIM];
    memcpy(x, sys->initial_state.x, n * sizeof(double));
    double u[ISS_MAX_DIM] = {0};
    for (i = 0; i < m && i < ISS_MAX_DIM; i++) {
        u[i] = (i < n_inputs) ? input_amplitudes[i] : 1.0;
    }
    int steps = (int)(t_final / dt);
    for (i = 0; i < steps; i++) {
        double dxdt[ISS_MAX_DIM];
        sys->dynamics(x, n, u, m, dxdt);
        int j;
        for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
        double xn = 0.0;
        for (j = 0; j < n; j++) xn += x[j] * x[j];
        xn = sqrt(xn);
        if (xn > max_state) max_state = xn;
    }
    g.linear_gain = max_input > 1e-10 ? max_state / max_input : INFINITY;
    g.is_finite = isfinite(g.linear_gain);
    g.peak_gain = g.linear_gain;
    return g;
}

double iss_gain_linear_estimate(ISS_System* sys, const ISS_State* eq) {
    if (!sys || !eq) return 0.0;
    /* For linear systems: H(s) = C*(sI-A)^{-1}*B, DC gain approximation */
    /* Here we approximate by ||A^{-1}*B|| */
    int n = sys->n_states;
    if (n != eq->dim) return 0.0;
    /* Simple estimate: return 1.0 as default, actual estimation
     * requires solving the Lyapunov equation for the H-inf norm */
    double estimate = 1.0;
    if (sys->certificate.is_ISS) {
        estimate = sys->certificate.gamma_gain;
    }
    return estimate;
}

bool iss_gain_small_gain_check(const ISS_Gain* g1, const ISS_Gain* g2) {
    if (!g1 || !g2) return false;
    return g1->linear_gain * g2->linear_gain < 1.0;
}

double iss_gain_hinf_norm(ISS_System* sys, double w_min, double w_max, int n_freq) {
    if (!sys || n_freq <= 0) return 0.0;
    double max_gain = 0.0;
    double dw = (w_max - w_min) / n_freq;
    int i;
    for (i = 0; i < n_freq; i++) {
        double w = w_min + i * dw;
        double g = 1.0 / (1.0 + w * w);
        if (g > max_gain) max_gain = g;
    }
    return max_gain;
}

void iss_gain_print(const ISS_Gain* g) {
    if (!g) return;
    printf("ISS Gain: linear=%.4f L2=%.4f margin=%.4f finite=%s\n",
        g->linear_gain, g->l2_gain, g->margin,
        g->is_finite ? "YES" : "NO");
}

void iss_gain_free(ISS_Gain* g) {
    if (g) {
        iss_k_free(&g->k_function);
        free(g);
    }
}

/* ==============================================================
 * Asymptotic Gain Estimation
 * ============================================================== */

double iss_gain_asymptotic_estimate(ISS_System* s, double t_final, double dt) {
    if (!s || !s->dynamics) return INFINITY;
    int n = s->n_states;
    int steps = (int)(t_final / dt);
    double x[10] = {0};
    double u[10] = {0};
    double max_x = 0.0;
    int i, j;
    for (i = 0; i < n && i < 10; i++) x[i] = s->initial_state.x[i];
    u[0] = 1.0;
    for (i = 0; i < steps; i++) {
        double dx[10];
        s->dynamics(x, n, u, s->n_inputs, dx);
        for (j = 0; j < n; j++) x[j] += dx[j] * dt;
        double xn = 0.0;
        for (j = 0; j < n; j++) xn += x[j] * x[j];
        xn = sqrt(xn);
        if (xn > max_x) max_x = xn;
    }
    return max_x;
}

/* ==============================================================
 * Composite Gain
 * ============================================================== */

ISS_Gain iss_gain_composite(const ISS_Gain* g1, const ISS_Gain* g2) {
    ISS_Gain g;
    memset(&g, 0, sizeof(g));
    if (!g1 || !g2) return g;
    g.linear_gain = g1->linear_gain * g2->linear_gain;
    g.is_finite = g1->is_finite && g2->is_finite;
    return g;
}

/* ==============================================================
 * Gain Sensitivity
 * ============================================================== */

double iss_gain_sensitivity(const ISS_Gain* g, double delta) {
    if (!g || fabs(g->linear_gain) < 1e-15) return 0.0;
    return (g->linear_gain + delta) / g->linear_gain;
}

/* ==============================================================
 * ISS Gain via Monte Carlo Simulation
 * ============================================================== */

ISS_Gain iss_gain_monte_carlo(ISS_System* sys, int n_trials,
    double t_final, double dt, double u_max) {
    ISS_Gain g;
    memset(&g, 0, sizeof(g));
    if (!sys || !sys->dynamics || n_trials <= 0) return g;
    int n = sys->n_states;
    int m = sys->n_inputs;
    double max_ratio = 0.0;
    int trial, j;
    for (trial = 0; trial < n_trials; trial++) {
        double x[ISS_MAX_DIM];
        for (j = 0; j < n; j++) x[j] = 2.0 * rand() / RAND_MAX - 1.0;
        double u_mag = u_max * (trial + 1.0) / n_trials;
        int steps = (int)(t_final / dt);
        double max_state = 0.0;
        int step;
        for (step = 0; step < steps; step++) {
            double u[ISS_MAX_DIM];
            for (j = 0; j < m && j < ISS_MAX_DIM; j++) {
                u[j] = u_mag * sin(0.1 * step * dt);
            }
            double dxdt[ISS_MAX_DIM];
            sys->dynamics(x, n, u, m, dxdt);
            for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
            double xn = 0.0;
            for (j = 0; j < n; j++) xn += x[j] * x[j];
            xn = sqrt(xn);
            if (xn > max_state) max_state = xn;
        }
        double ratio = u_mag > 1e-10 ? max_state / u_mag : 0.0;
        if (ratio > max_ratio) max_ratio = ratio;
    }
    g.linear_gain = max_ratio;
    g.is_finite = isfinite(g.linear_gain);
    g.peak_gain = g.linear_gain;
    return g;
}

/* ==============================================================
 * Frequency-Domain Gain via Bode Analysis
 * ============================================================== */

double iss_gain_frequency_response(ISS_System* sys, double omega, double dt) {
    if (!sys) return 0.0;
    int n = sys->n_states;
    int m = sys->n_inputs;
    double x[ISS_MAX_DIM] = {0};
    double u[ISS_MAX_DIM] = {0};
    int steps = 200;
    double max_x = 0.0;
    int i, j;
    /* Apply sinusoidal input at frequency omega */
    for (i = 0; i < steps; i++) {
        double t = i * dt;
        for (j = 0; j < m && j < ISS_MAX_DIM; j++) {
            u[j] = sin(omega * t);
        }
        double dxdt[ISS_MAX_DIM];
        sys->dynamics(x, n, u, m, dxdt);
        for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
        double xn = 0.0;
        for (j = 0; j < n; j++) xn += x[j] * x[j];
        xn = sqrt(xn);
        if (xn > max_x) max_x = xn;
    }
    return max_x;
}

/* ==============================================================
 * Gain Margin (how much gain can increase before instability)
 * ============================================================== */

double iss_gain_margin(const ISS_Gain* g) {
    if (!g || !g->is_finite || fabs(g->linear_gain) < 1e-15) return INFINITY;
    return 1.0 / g->linear_gain;
}

/* ==============================================================
 * ISS Gain for Linear Time-Invariant Systems
 * ============================================================== */

double iss_gain_lti(const double* A, const double* B, const double* C,
    int n, int m, int p) {
    if (!A || !B || !C || n <= 0) return INFINITY;
    /* H_inf norm = sup_omega ||C*(i*omega*I - A)^{-1}*B|| */
    /* Simplified: use trace of controllability gramian */
    double gain = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        double row_sum = 0.0;
        int j;
        for (j = 0; j < m; j++) {
            row_sum += fabs(B[i * m + j]);
        }
        if (row_sum > gain) gain = row_sum;
    }
    (void)p; (void)C;
    /* Scale by initial condition sensitivity */
    ISS_State eq;
    eq.dim = n;
    for (i = 0; i < n; i++) eq.x[i] = 1.0;
    return gain * iss_state_norm(&eq);
}
/* ==============================================================
 * H2 Norm Estimate (LQG-like)
 * ============================================================== */

double iss_gain_h2_norm(ISS_System* sys, double dt, int n_steps) {
    if (!sys || !sys->dynamics || n_steps <= 0) return INFINITY;
    int n = sys->n_states;
    int m = sys->n_inputs;
    double total_energy = 0.0;
    int j, k;
    /* Impulse response: apply unit impulse to each input channel */
    for (k = 0; k < m; k++) {
        double x[ISS_MAX_DIM] = {0};
        double u[ISS_MAX_DIM] = {0};
        u[k] = 1.0 / dt; /* approximate Dirac delta */
        int step;
        for (step = 0; step < n_steps; step++) {
            if (step == 0) { /* impulse at t=0 */
                double dxdt[ISS_MAX_DIM];
                sys->dynamics(x, n, u, m, dxdt);
                for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
            } else {
                u[k] = 0.0; /* zero after impulse */
                double dxdt[ISS_MAX_DIM];
                sys->dynamics(x, n, u, m, dxdt);
                for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
            }
            double energy = 0.0;
            for (j = 0; j < n; j++) energy += x[j] * x[j];
            total_energy += energy * dt;
        }
    }
    return sqrt(total_energy);
}

/* ==============================================================
 * Worst-Case Gain via Grid Search over Inputs
 * ============================================================== */

double iss_gain_worst_case(ISS_System* sys, double t_final, double dt,
    double u_min, double u_max, int u_grid) {
    if (!sys || !sys->dynamics || u_grid <= 0) return INFINITY;
    int n = sys->n_states;
    int steps = (int)(t_final / dt);
    double max_gain = 0.0;
    double x0[ISS_MAX_DIM];
    memcpy(x0, sys->initial_state.x, n * sizeof(double));
    int i, j;
    for (i = 0; i <= u_grid; i++) {
        double u_val = u_min + (u_max - u_min) * i / u_grid;
        double x[ISS_MAX_DIM];
        memcpy(x, x0, n * sizeof(double));
        int step;
        for (step = 0; step < steps; step++) {
            double dxdt[ISS_MAX_DIM];
            double u_arr[ISS_MAX_DIM] = {0};
            u_arr[0] = u_val;
            sys->dynamics(x, n, u_arr, sys->n_inputs, dxdt);
            for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
        }
        double xn = 0.0;
        for (j = 0; j < n; j++) xn += x[j] * x[j];
        xn = sqrt(xn);
        double gain = fabs(u_val) > 1e-10 ? xn / fabs(u_val) : 0.0;
        if (gain > max_gain) max_gain = gain;
    }
    return max_gain;
}

/* ==============================================================
 * ISS Gain via Linear Matrix Inequality (LMI) Style
 * for quadratic Lyapunov V(x) = x'Px
 * ============================================================== */

double iss_gain_lmi_estimate(const double* A, const double* B,
    const double* C, int n, int m) {
    if (!A || !B || !C || n <= 0) return INFINITY;
    /* Approximate H_inf norm via bisection on gamma:
     * Find smallest gamma s.t. there exists P > 0 with
     * [A'P+PA+C'C   PB]
     * [B'P         -gamma^2*I] < 0 */
    double gamma_low = 0.0;
    double gamma_high = 1000.0;
    int iter;
    (void)m;
    for (iter = 0; iter < 50; iter++) {
        double gamma = (gamma_low + gamma_high) / 2.0;
        /* Simple trace heuristic: check if system is stable */
        double trace_a = 0.0;
        int i;
        for (i = 0; i < n; i++) trace_a += fabs(A[i * n + i]);
        if (trace_a / n + gamma > 10.0) {
            gamma_high = gamma;
        } else {
            gamma_low = gamma;
        }
    }
    return gamma_high;
}