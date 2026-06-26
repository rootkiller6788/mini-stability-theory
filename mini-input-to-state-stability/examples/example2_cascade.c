/*
 * Example 2: Cascade ISS System
 *
 * System:  dx1/dt = -x1 + x2      (S1: ISS w.r.t. x2, gain gamma1=1)
 *          dx2/dt = -2*x2 + u      (S2: GAS with ISS gain gamma2=0.5)
 *
 * Cascade ISS property: If S2 is GAS and S1 is ISS w.r.t. x2,
 * then the cascade (x1, x2) is GAS.
 *
 * This example demonstrates:
 * 1. Creating interconnected ISS systems
 * 2. Verifying cascade ISS property
 * 3. Computing cascade ISS certificate
 * 4. Trajectory-based verification
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "iss_core.h"
#include "iss_interconnected.h"

static void s1_dynamics(const double* x, int n, const double* u, int m, double* d) {
    (void)n; (void)m;
    d[0] = -x[0] + u[0]; /* S1: dx1/dt = -x1 + x2 */
}

static double s1_lyap(const double* x, int n) {
    (void)n;
    return 0.5 * x[0] * x[0];
}

static void s1_lyap_grad(const double* x, int n, double* g) {
    (void)n;
    g[0] = x[0];
}

static void s2_dynamics(const double* x, int n, const double* u, int m, double* d) {
    (void)n; (void)m;
    d[0] = -2.0 * x[0] + u[0]; /* S2: dx2/dt = -2*x2 + u */
}

static double s2_lyap(const double* x, int n) {
    (void)n;
    return 0.5 * x[0] * x[0];
}

static void s2_lyap_grad(const double* x, int n, double* g) {
    (void)n;
    g[0] = x[0];
}

int main(void) {
    printf("=== Example 2: Cascade ISS System ===\n\n");

    /* Create subsystems */
    ISS_System* s1 = iss_system_create("driven", 1, 1);
    s1->dynamics = s1_dynamics;
    s1->lyapunov = s1_lyap;
    s1->lyapunov_grad = s1_lyap_grad;
    s1->initial_state.x[0] = 1.0;
    s1->initial_state.dim = 1;

    ISS_System* s2 = iss_system_create("driving", 1, 1);
    s2->dynamics = s2_dynamics;
    s2->lyapunov = s2_lyap;
    s2->lyapunov_grad = s2_lyap_grad;
    s2->initial_state.x[0] = 2.0;
    s2->initial_state.dim = 1;

    printf("S1: dx1/dt = -x1 + x2,  Lyapunov: V1 = 0.5*x1^2\n");
    printf("S2: dx2/dt = -2*x2 + u, Lyapunov: V2 = 0.5*x2^2\n\n");

    /* Build cascade */
    ISS_Cascade* cascade = iss_cascade_create(s1, s2);
    printf("Cascade ISS: %s\n", iss_cascade_is_ISS(cascade) ? "YES" : "NO");

    /* Verify cascade ISS with explicit gains */
    iss_cascade_verify_ISS(cascade, 1e-6);
    double gain = iss_cascade_gain(cascade);
    printf("Cascade gain: %.4f\n", gain);

    /* Certificate */
    ISS_Certificate cert = iss_cascade_certificate(cascade);
    printf("Certificate: alpha=%.4f sigma=%.4f gamma=%.4f\n",
        cert.alpha_coeff, cert.sigma_coeff, cert.gamma_gain);

    /* Trajectory simulation of cascade */
    printf("\nCascade Trajectory (u=0.1 constant):\n");
    printf("  t       x1(t)      x2(t)\n");
    printf("  ------  ---------  ---------\n");

    double x1 = 1.0, x2 = 2.0;
    double u_ext = 0.1;
    double dt = 0.1;
    int i;
    for (i = 0; i <= 20; i++) {
        double t = i * dt;
        printf("  %6.2f  %9.6f  %9.6f\n", t, x1, x2);

        /* S2 dynamics: dx2/dt = -2*x2 + u */
        double dx2 = -2.0 * x2 + u_ext;
        x2 += dx2 * dt;

        /* S1 dynamics: dx1/dt = -x1 + x2 */
        double dx1 = -x1 + x2;
        x1 += dx1 * dt;
    }

    /* Composite Lyapunov */
    printf("\nComposite Lyapunov: V = V1 + V2\n");
    printf("  V(0) = %.4f\n", s1_lyap(&x1, 1) + s2_lyap(&x2, 0));

    iss_cascade_free(cascade);
    iss_system_free(s1);
    iss_system_free(s2);
    printf("\n=== Example 2 Complete ===\n");
    return 0;
}
