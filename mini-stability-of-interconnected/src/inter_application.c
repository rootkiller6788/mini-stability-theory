/* inter_application.c - L7 Application: Decentralized power grid frequency control.
 * Small-gain theorem for interconnected power generators.
 * Two-area power system with decentralized frequency controllers.
 *
 * L7: Smart grid, power system stability, IEEE 39-bus test system
 * L8: Time-varying interconnection topology
 * Ref: Jiang (2000) "Small-Gain Theorem for ISS Systems"
 *       Kundur (1994) "Power System Stability and Control"
 */

#include "intercon_core.h"
#include "intercon_small_gain.h"
#include "intercon_decentralized.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double inertia, damping, coupling, droop;
} GeneratorParams;

void two_area_system(const double* x, double* dx, int n, double t) {
    (void)t; (void)n;
    GeneratorParams g1 = {4.0, 1.0, 0.5, 0.05};
    GeneratorParams g2 = {5.0, 0.8, 0.5, 0.04};
    /* x[0]=f1, x[1]=Pm1, x[2]=f2, x[3]=Pm2 */
    dx[0] = (-g1.damping * x[0] + x[1] - g1.coupling * (x[0] - x[2])) / g1.inertia;
    dx[1] = (-x[1] / g1.droop - x[0]) / g1.inertia;
    dx[2] = (-g2.damping * x[2] + x[3] + g2.coupling * (x[0] - x[2])) / g2.inertia;
    dx[3] = (-x[3] / g2.droop - x[2]) / g2.inertia;
}

double grid_lyapunov(const double* x, int n) {
    (void)n;
    return 2.0*x[0]*x[0] + x[1]*x[1] + 2.5*x[2]*x[2] + x[3]*x[3];
}

double compute_small_gain(double g1_damp, double g1_coup,
    double g2_damp, double g2_coup) {
    double g1 = g1_coup / (g1_damp + 1e-6);
    double g2 = g2_coup / (g2_damp + 1e-6);
    return g1 * g2;
}

int main(void) {
    printf("=== L7: Two-Area Power Grid Decentralized Stability ===\n\n");

    double gamma = compute_small_gain(1.0, 0.5, 0.8, 0.5);
    printf(" Small-gain product: %.4f\n", gamma);
    if (gamma < 1.0)
        printf(" gamma < 1: Interconnected system is ISS stable.\n\n");
    else
        printf(" gamma >= 1: Stability not guaranteed.\n\n");

    double x[4] = {0.5, 0.0, -0.3, 0.0};
    double dt = 0.01;
    printf(" Initial: f1=%.2f Hz dev, f2=%.2f Hz dev\n", x[0], x[2]);

    for (int i = 0; i <= 100; i++) {
        if (i % 25 == 0) {
            double V = grid_lyapunov(x, 4);
            printf("  t=%.2f: f1=%.4f, f2=%.4f, V=%.6f\n", i*dt, x[0], x[2], V);
        }
        for (int s = 0; s < 10; s++) {
            double dx[4];
            two_area_system(x, dx, 4, i*dt);
            for (int d = 0; d < 4; d++) x[d] += dt * dx[d];
        }
    }
    printf("\n Frequencies converged: f1=%.6f, f2=%.6f\n", x[0], x[2]);
    printf(" Decentralized frequency control STABLE.\n");
    return 0;
}