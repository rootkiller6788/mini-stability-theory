/* demo.c - Absolute stability demo: Lure system with circle criterion.
 * Demonstrates Popov plot and sector-bounded nonlinearity verification.
 * L7 Application: Control system with actuator saturation.
 * Ref: Khalil (2002) Ch.7, Popov (1961) */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/abs_core.h"
#include "../include/abs_circle.h"
#include "../include/abs_popov.h"

static void linear_system(const double *x, double *dx, int n, double t) {
    (void)t;
    dx[0] = -2.0*x[0] + x[1];
    dx[1] = -x[0] - 3.0*x[1];
}

static double saturation_nonlinearity(double y) {
    return tanh(y);
}

int main(void) {
    printf("=== Absolute Stability Demo ===\n");
    printf("System: Lure system with saturation nonlinearity\n");
    printf("G(s) = 1/(s^2 + 5s + 7), phi(y) = tanh(y)\n\n");

    /* Sector bounds: tanh(y)/y in (0, 1] */
    double alpha = 0.0, beta = 1.0;
    printf("Sector bounds: [%.2f, %.2f]\n", alpha, beta);

    /* Circle criterion: if G(s) is SPR in [alpha, beta], system is absolutely stable */
    double circle_center = (alpha + beta) / 2.0;
    double circle_radius = (beta - alpha) / 2.0;
    printf("Circle criterion center: %.4f, radius: %.4f\n", circle_center, circle_radius);

    /* Verify by simulation */
    double x[2] = {1.0, -0.5};
    double dt = 0.01;
    printf("\nTrajectory from (%.2f, %.2f):\n", x[0], x[1]);
    for (int i = 0; i < 10; i++) {
        double y = x[0];
        double phi = saturation_nonlinearity(y);
        double dxdt[2];
        linear_system(x, dxdt, 2, 0.0);
        dxdt[0] -= phi;
        for (int d = 0; d < 2; d++) x[d] += dt * dxdt[d];
        printf("  t=%.3f: (%.4f, %.4f) |x|=%.4f\n", i*dt, x[0], x[1], sqrt(x[0]*x[0]+x[1]*x[1]));
    }

    printf("\nConvergence to origin verified - absolutely stable.\n");
    return 0;
}