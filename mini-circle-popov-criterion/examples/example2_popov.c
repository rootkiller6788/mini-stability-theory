#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "circle_criterion.h"
#include "popov_criterion.h"
#include "lure_system.h"

/* Example 2: Popov Criterion with Optimal Multiplier.
 * G(s) = 1/(s^2 + s + 1), phi in sector [0, k].
 * Demonstrates that Popov can be less conservative than circle. */
int main(void) {
    printf("=== Example 2: Popov Criterion ===\n");

    TransferFunction* tf = tf_create(1, 3);
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0; tf->den[2] = 1.0;

    printf("G(s) = 1/(s^2 + s + 1)\n\n");

    printf("%-10s %-15s %-15s %-10s\n",
           "k", "Circle", "Popov", "PopovEta");
    printf("-----------------------------------------------\n");

    double k_values[] = {0.5, 1.0, 2.0, 3.0, 5.0, 10.0};
    for (int i = 0; i < 6; i++) {
        double k = k_values[i];
        SectorBounds sb = {0.0, k, 0, 0, 0.0, 0.0};

        double cm;
        int circle_ok = circle_criterion_check(tf, sb, &cm);

        double eta_opt, pm;
        popov_find_optimal_slope(tf, sb, &eta_opt, &pm);
        int popov_ok = (pm > 0.0) ? 1 : 0;

        printf("k=%-8.1f %-15s %-15s eta=%.2f\n",
               k,
               circle_ok ? "STABLE" : "inconclusive",
               popov_ok  ? "STABLE" : "inconclusive",
               eta_opt);
    }

    /* Compare criteria */
    printf("\n--- Comparison ---\n");
    double cm2, pm2, ef2;
    SectorBounds sb = {0.0, 3.0, 0, 0, 0.0, 0.0};
    int comp = popov_compare_criteria(tf, sb, &cm2, &pm2, &ef2);
    printf("Comparison result: Circle=%s Popov=%s (eta=%.2f)\n",
           (comp & 1) ? "OK" : "FAIL",
           (comp & 2) ? "OK" : "FAIL", ef2);

    tf_free(tf);
    return 0;
}
