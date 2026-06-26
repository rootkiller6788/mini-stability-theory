/*
 * Example 3: Small-Gain Theorem for ISS
 *
 * Two ISS systems in feedback:
 *   S1: dx1/dt = -a1*x1 + k12*h2(x2) + u1    (ISS gain gamma1 = k12/a1)
 *   S2: dx2/dt = -a2*x2 + k21*h1(x1) + u2    (ISS gain gamma2 = k21/a2)
 *
 * Nonlinear small-gain condition:
 *   gamma1(gamma2(s)) < s  for all s > 0
 *
 * For linear gains: gamma1 * gamma2 < 1
 *
 * This example demonstrates:
 * 1. Computing ISS gains for feedback interconnection
 * 2. Small-gain condition verification (linear and nonlinear)
 * 3. Cyclic small-gain for 3 interconnected systems
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "iss_core.h"
#include "iss_gain.h"
#include "iss_small_gain.h"

int main(void) {
    printf("=== Example 3: Small-Gain Theorem ===\n\n");

    /* System parameters */
    double a1 = 1.0, k12 = 0.3; /* S1 */
    double a2 = 2.0, k21 = 0.5; /* S2 */

    double gamma1 = k12 / a1;
    double gamma2 = k21 / a2;

    printf("S1: dx1/dt = -%.1f*x1 + %.1f*h2(x2) + u1\n", a1, k12);
    printf("S2: dx2/dt = -%.1f*x2 + %.1f*h1(x1) + u2\n", a2, k21);
    printf("ISS Gains: gamma1 = %.4f, gamma2 = %.4f\n", gamma1, gamma2);
    printf("Product: gamma1*gamma2 = %.4f\n", gamma1 * gamma2);

    /* Check small-gain */
    ISS_Gain g1, g2;
    memset(&g1, 0, sizeof(g1));
    memset(&g2, 0, sizeof(g2));
    g1.linear_gain = gamma1;
    g2.linear_gain = gamma2;
    g1.is_finite = true;
    g2.is_finite = true;

    bool linear_sg = iss_small_gain_condition(&g1, &g2);
    printf("\nLinear small-gain (g1*g2 < 1): %s\n", linear_sg ? "PASS" : "FAIL");

    /* Margin */
    const ISS_Gain* gs[] = {&g1, &g2};
    double margin = iss_small_gain_margin(gs, 2);
    printf("Small-gain margin: %.4f\n", margin);

    /* Nonlinear small-gain check */
    ISS_KFunction kg1 = iss_k_linear(gamma1);
    ISS_KFunction kg2 = iss_k_linear(gamma2);
    bool nonlinear_sg = iss_nonlinear_small_gain(&kg1, &kg2, 10.0, 100);
    printf("Nonlinear small-gain (g1(g2(s)) < s): %s\n", nonlinear_sg ? "PASS" : "FAIL");

    /* Compare g1(g2(s)) vs s */
    printf("\nNonlinear small-gain comparison:\n");
    printf("  s     g1(g2(s))  ID(s)\n");
    printf("  ----  ---------  -----\n");
    int i;
    for (i = 1; i <= 5; i++) {
        double s = (double)i;
        double composed = iss_k_eval(&kg1, iss_k_eval(&kg2, s));
        printf("  %.1f   %.6f    %.1f\n", s, composed, s);
    }

    /* Cyclic small-gain for 3 systems */
    printf("\n--- Cyclic Small-Gain (3 systems) ---\n");
    ISS_Gain g3;
    memset(&g3, 0, sizeof(g3));
    g3.linear_gain = 0.4;
    g3.is_finite = true;

    const ISS_Gain* gs3[] = {&g1, &g2, &g3};
    bool cyclic_sg = iss_cyclic_small_gain(gs3, 3);
    double prod = g1.linear_gain * g2.linear_gain * g3.linear_gain;
    printf("Gains: %.2f * %.2f * %.2f = %.4f\n",
        g1.linear_gain, g2.linear_gain, g3.linear_gain, prod);
    printf("Cyclic small-gain: %s\n", cyclic_sg ? "PASS" : "FAIL");

    /* Spectral radius check */
    printf("\n--- Spectral Radius Small-Gain ---\n");
    double M[9] = {
        0.0,   0.15,  0.0,
        0.0,   0.0,   0.25,
        0.2,   0.0,   0.0
    };
    double rho = iss_spectral_radius(M, 3);
    printf("Gain matrix spectral radius: %.6f\n", rho);
    printf("Spectral radius < 1: %s\n", rho < 1.0 ? "ISS" : "NOT ISS");

    iss_k_free(&kg1);
    iss_k_free(&kg2);
    printf("\n=== Example 3 Complete ===\n");
    return 0;
}
