#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "circle_criterion.h"
#include "popov_criterion.h"
#include "kyp_lemma.h"
#include "lure_system.h"
#include "sector_theory.h"

/* Example 3: Complete Criterion Comparison with State-Space.
 * Builds a Lur'e system from state-space, tests multiple methods. */
int main(void) {
    printf("=== Example 3: Full Criterion Comparison ===\n\n");

    /* Build state-space: G(s) = (s+2)/(s^2 + 0.5s + 1)
     * Controllable canonical form:
     * A = [0 1; -1 -0.5], B = [0; 1], C = [2 1], D = 0 */
    StateSpace* ss = ss_create(2);
    double A[4] = {0.0, 1.0, -1.0, -0.5};
    double B[2] = {0.0, 1.0};
    double C[2] = {2.0, 1.0};
    ss_set_matrices(ss, A, B, C, 0.0);

    printf("State-space system (n=2):\n");
    printf("  A = [%.2f %.2f; %.2f %.2f]\n",
           A[0], A[1], A[2], A[3]);
    printf("  B = [%.2f; %.2f]\n", B[0], B[1]);
    printf("  C = [%.2f %.2f]\n", C[0], C[1]);

    /* Convert to transfer function */
    TransferFunction* tf = tf_create(3, 3);
    kyp_to_transfer_function(ss, tf);
    printf("  G(s) = (%.2f s^2 + %.2f s + %.2f) / (%.2f s^2 + %.2f s + %.2f)\n",
           tf->num[2], tf->num[1], tf->num[0],
           tf->den[2], tf->den[1], tf->den[0]);

    /* Test stability for sector [0, k] with various k */
    double k_vals[] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0};
    printf("\n%s\n", "-------------------------------------------");
    printf("%-8s %-10s %-10s %-10s\n",
           "k", "Circle", "Popov", "KYP");
    printf("%s\n", "-------------------------------------------");

    for (int i = 0; i < 6; i++) {
        double k = k_vals[i];
        SectorBounds sb = {0.0, k, 0, 0, 0.0, 0.0};

        double cm, pm, eta;
        int c_ok = circle_criterion_check(tf, sb, &cm);
        popov_find_optimal_slope(tf, sb, &eta, &pm);
        int p_ok = (pm > 0.0) ? 1 : 0;

        KYPSolution* sol = kyp_solve(ss, 0.01);
        int k_ok = sol ? (sol->is_feasible ? 1 : 0) : 0;
        kyp_free(sol);

        printf("k=%-6.1f %-10s %-10s %-10s\n",
               k,
               c_ok ? "OK" : "FAIL",
               p_ok ? "OK" : "FAIL",
               k_ok ? "OK" : "FAIL");
    }

    /* Sector theory: find max allowed gain */
    printf("\n--- Max Sector via Circle Criterion ---\n");
    double max_k, mk_margin;
    max_k = sector_max_allowed_gain(tf, 0, &mk_margin);
    printf("Max k (circle): %.2f\n", max_k);

    /* Sector transformations */
    printf("\n--- Sector Transformations ---\n");
    SectorBounds orig = {1.0, 4.0, 0, 0, 0.0, 0.0};
    SectorBounds shifted = sector_shift(orig, -1.0);
    printf("Shift [1,4] by -1 => [%.1f, %.1f]\n",
           shifted.alpha, shifted.beta);

    SectorBounds scaled = sector_scale(orig, 0.5);
    printf("Scale [1,4] by 0.5 => [%.1f, %.1f]\n",
           scaled.alpha, scaled.beta);

    /* Cleanup */
    tf_free(tf);
    ss_free(ss);
    printf("\nDone.\n");
    return 0;
}
