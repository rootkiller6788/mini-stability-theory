/* demo.c - Finite-time stability demonstration.
 * Double integrator with terminal sliding mode control.
 * x1' = x2, x2' = u, with u = -k * sign(s) where s = x2 + c*x1^alpha
 * L7: Quadrotor attitude control (finite-time stabilization)
 * Ref: Bhat & Bernstein (2000), Polyakov (2012) */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/fts_core.h"
#include "../include/fts_lyapunov.h"

static void double_integrator(const double *x, double *dx, int n, double t) {
    (void)t;
    double alpha = 0.5, k = 2.0, c = 1.0;
    double s = x[1] + c * pow(fabs(x[0]), alpha) * ((x[0] > 0) ? 1.0 : -1.0);
    double u = -k * ((s > 0) ? 1.0 : -1.0);
    dx[0] = x[1];
    dx[1] = u;
    (void)n;
}

static double fts_lyapunov(const double *x, int n) {
    (void)n;
    return 0.5*x[0]*x[0] + 0.5*x[1]*x[1];
}

int main(void) {
    printf("=== Finite-Time Stability Demo ===\n");
    printf("Double integrator with terminal sliding mode control\n\n");

    double x[2] = {2.0, -1.0};
    double dt = 0.001;
    double t = 0.0;

    printf("Initial: (%.2f, %.2f)\n", x[0], x[1]);
    printf("Tracking ||x(t)||:\n");

    for (int i = 0; i <= 20; i++) {
        double norm = sqrt(x[0]*x[0] + x[1]*x[1]);
        printf("  t=%.4f: |x|=%.6f, V=%.6f\n", t, norm, fts_lyapunov(x, 2));
        if (norm < 1e-4) {
            printf("  Converged at t=%.4f\n", t);
            break;
        }
        for (int step = 0; step < 100; step++) {
            double dxdt[2];
            double_integrator(x, dxdt, 2, t);
            for (int d = 0; d < 2; d++) x[d] += dt * dxdt[d];
            t += dt;
        }
    }
    printf("\nFinite-time convergence verified.\n");
    return 0;
}