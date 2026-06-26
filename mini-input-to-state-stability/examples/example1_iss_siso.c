/*
 * Example 1: SISO (Single-Input Single-Output) ISS System
 *
 * System: dx/dt = -x + u  (linear, ISS with gain gamma=1)
 * Lyapunov: V(x) = 0.5*x^2
 *
 * ISS bound: |x(t)| <= beta(|x0|, t) + gamma(||u||_inf)
 * where beta(s,t) = s*exp(-t) and gamma(s) = s
 *
 * This example demonstrates:
 * 1. Creating an ISS system with dynamics and Lyapunov function
 * 2. Verifying ISS via random sampling
 * 3. Computing the ISS gain
 * 4. Simulating the system response with a bounded input
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "iss_core.h"
#include "iss_lyapunov.h"
#include "iss_gain.h"
#include "iss_verify.h"

static void siso_dynamics(const double* x, int n, const double* u, int m, double* d) {
    (void)n; (void)m;
    d[0] = -x[0] + u[0];
}

static double siso_lyap(const double* x, int n) {
    (void)n;
    return 0.5 * x[0] * x[0];
}

static void siso_lyap_grad(const double* x, int n, double* g) {
    (void)n;
    g[0] = x[0];
}

int main(void) {
    printf("=== Example 1: SISO ISS System ===\n\n");

    /* Create system */
    ISS_System* sys = iss_system_create("siso_iss", 1, 1);
    sys->dynamics = siso_dynamics;
    sys->lyapunov = siso_lyap;
    sys->lyapunov_grad = siso_lyap_grad;
    sys->initial_state.x[0] = 2.0; /* x(0) = 2 */
    sys->initial_state.dim = 1;

    printf("System: dx/dt = -x + u, x(0) = %.1f\n", sys->initial_state.x[0]);

    /* Verify ISS */
    ISS_VerifyConfig cfg = iss_verify_config_default();
    cfg.monte_carlo_samples = 500;
    ISS_Verification ver = iss_verify(sys, &cfg);

    printf("ISS Verification: %s\n", ver.is_ISS ? "PASSED" : "FAILED");
    printf("  Checks: %d/%d passed\n", ver.passed_checks, ver.total_checks);
    printf("  ISS Gain: %.4f\n", ver.gain.linear_gain);

    /* Compute ISS bound */
    ISS_KFunction gamma = iss_k_linear(ver.gain.linear_gain);
    ISS_KLFunction beta = iss_kl_exponential(1.0);

    double u_inf = 0.5; /* constant input */
    double x0_norm = iss_state_norm(&sys->initial_state);
    printf("\nISS Bound Analysis:\n");
    printf("  |x0| = %.2f, ||u||_inf = %.2f\n", x0_norm, u_inf);

    int i;
    for (i = 0; i <= 10; i++) {
        double t = (double)i * 0.5;
        double bound = iss_bound_full(&beta, &gamma, x0_norm, t, u_inf);
        printf("  t=%.1f: |x(t)| <= %.4f\n", t, bound);
    }

    /* Simulate trajectory */
    printf("\nTrajectory Simulation (u=0.5 constant):\n");
    double x = 2.0;
    double u_val = 0.5;
    for (i = 0; i <= 10; i++) {
        double t = (double)i * 0.5;
        printf("  t=%.1f: x=%.6f\n", t, x);
        /* Euler step */
        double dxdt;
        siso_dynamics(&x, 1, &u_val, 1, &dxdt);
        x += dxdt * 0.5;
    }

    iss_k_free(&gamma);
    free(beta.param);
    iss_system_free(sys);
    printf("\n=== Example 1 Complete ===\n");
    return 0;
}
