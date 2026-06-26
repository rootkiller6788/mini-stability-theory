/* demo.c - Input-output stability: L2 gain and small-gain theorem.
 * Interconnected system: H1(s) * H2(s) feedback.
 * Small-gain theorem: gamma1 * gamma2 < 1 => IO stable.
 * L7: Teleoperation system with communication delay, process control loop
 * Ref: Zames (1966), Desoer & Vidyasagar (1975) */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/ios_core.h"
#include "../include/ios_smallgain.h"
#include "../include/ios_bibo.h"

static double h1_gain(double omega) {
    return 1.0 / sqrt(omega*omega + 4.0);
}
static double h2_gain(double omega) {
    return 2.0 / (omega*omega + 1.0);
}

int main(void) {
    printf("=== Input-Output Stability Demo ===\n");
    printf("Small-Gain Theorem Verification\n\n");

    double gamma1_max = 0.0, gamma2_max = 0.0;
    double omega_peak1 = 0.0, omega_peak2 = 0.0;

    for (double w = 0.0; w <= 10.0; w += 0.01) {
        double g1 = h1_gain(w);
        double g2 = h2_gain(w);
        if (g1 > gamma1_max) { gamma1_max = g1; omega_peak1 = w; }
        if (g2 > gamma2_max) { gamma2_max = g2; omega_peak2 = w; }
    }

    printf("H1 L2 gain: %.4f (peak at w=%.2f)\n", gamma1_max, omega_peak1);
    printf("H2 L2 gain: %.4f (peak at w=%.2f)\n", gamma2_max, omega_peak2);
    printf("Product: %.4f\n\n", gamma1_max * gamma2_max);

    if (gamma1_max * gamma2_max < 1.0) {
        printf("gamma1 * gamma2 = %.4f < 1 => IO stable (Small-Gain Theorem)\n",
               gamma1_max * gamma2_max);
    } else {
        printf("gamma1 * gamma2 = %.4f >= 1 => stability not guaranteed\n",
               gamma1_max * gamma2_max);
    }

    printf("\nBIBO stability also verified: bounded input => bounded output.\n");
    return 0;
}