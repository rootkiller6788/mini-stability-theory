#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "circle_criterion.h"
#include "popov_criterion.h"
#include "lure_system.h"

/* Example 1: Circle Criterion on a 1st-order Lur'e system.
 * G(s) = 1/(s+1), phi in sector [0, k].
 * Demonstrates the half-plane case of the circle criterion. */
int main(void) {
    printf("=== Example 1: Circle Criterion ===\n");

    TransferFunction* tf = tf_create(1, 2);
    tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;

    printf("G(s) = 1/(s+1)\n");

    /* Test sectors [0, k] for various k values */
    double k_values[] = {0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0};
    int nk = sizeof(k_values) / sizeof(k_values[0]);

    for (int i = 0; i < nk; i++) {
        double k = k_values[i];
        SectorBounds sb = {0.0, k, 0, 0, 0.0, 0.0};
        double margin;
        int stable = circle_criterion_check(tf, sb, &margin);
        printf("  Sector [0, %.1f]: %s (margin=%.4f)\n",
               k, stable ? "STABLE" : "INCONCLUSIVE", margin);
    }

    /* Show critical disk parameters for sector [0.5, 3.0] */
    SectorBounds sb2 = {0.5, 3.0, 0, 0, 0.0, 0.0};
    printf("\nCritical disk for [0.5, 3.0]:\n");
    printf("  Center = %.4f\n", circle_criterion_center(sb2));
    printf("  Radius = %.4f\n", circle_criterion_radius(sb2));
    double margin;
    int st = circle_criterion_check(tf, sb2, &margin);
    printf("  Stability: %s (margin=%.4f)\n",
           st ? "ABSOLUTELY STABLE" : "INCONCLUSIVE", margin);

    tf_free(tf);
    return 0;
}
