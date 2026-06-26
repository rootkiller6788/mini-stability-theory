#include "ios_core.h"
#include <stdio.h>

/* Example 2: Small-Gain Theorem.
 * Demonstrates the Zames small-gain theorem:
 * If gamma1 * gamma2 < 1, the feedback interconnection is L2-stable.
 * Shows both scalar and system-level small-gain tests. */
int main(void) {
    printf("=== Example 2: Small-Gain Theorem ===\n\n");

    /* Scalar small-gain test */
    printf("Scalar Small-Gain Tests:\n");
    double pairs[][2] = {
        {0.3, 2.5},   /* 0.75 < 1 => stable */
        {0.5, 1.5},   /* 0.75 < 1 => stable */
        {0.8, 0.9},   /* 0.72 < 1 => stable */
        {0.9, 1.2},   /* 1.08 > 1 => fail */
        {1.5, 0.4}    /* 0.60 < 1 => stable */
    };
    for (int i = 0; i < 5; i++) {
        IOSSmallGainResult r = ios_small_gain_test(pairs[i][0], pairs[i][1]);
        printf("  g1=%.1f g2=%.1f product=%.4f: %s (margin=%.4f)\n",
               pairs[i][0], pairs[i][1], r.gamma,
               r.satisfied ? "STABLE" : "INCONCLUSIVE", r.margin);
    }

    /* System-level small-gain (transfer functions) */
    printf("\nSystem-Level Small-Gain:\n");
    IOSStateSpace* G = ios_ss_create(1, 1, 1);
    IOSStateSpace* H = ios_ss_create(1, 1, 1);
    if (G && H) {
        G->A[0] = -1.0; G->B[0] = 1.0; G->C[0] = 2.0;
        H->A[0] = -3.0; H->B[0] = 1.0; H->C[0] = 0.3;
        double gG = ios_L2_gain_ss(G);
        double gH = ios_L2_gain_ss(H);
        printf("  G: L2 gain = %.4f\n", gG);
        printf("  H: L2 gain = %.4f\n", gH);
        printf("  Product: %.4f -> %s\n", gG * gH,
               gG * gH < 1.0 ? "STABLE" : "INCONCLUSIVE");
    }
    ios_ss_free(G); ios_ss_free(H);

    return 0;
}
