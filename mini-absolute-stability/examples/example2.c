/* example2.c — Popov Criterion demonstration
 *
 * Demonstrates the Popov criterion on a 2nd-order Lur'e system
 * with time-invariant sector nonlinearity.
 *
 * The Popov criterion is sharper than the circle criterion for
 * time-invariant nonlinearities. We compare both criteria.
 *
 * System:  G(s) = (s+3)/((s+1)(s+2)(s+0.5))
 *          This is a 3rd-order minimum-phase system.
 */

#include "abs_core.h"
#include "abs_circle.h"
#include "abs_popov.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Example 2: Popov Criterion ===\n\n");

    /* 3rd-order system in companion form:
     * G(s) = (s+3) / (s^3 + 3.5*s^2 + 3.5*s + 1)
     *        = (s+3) / ((s+1)(s+2)(s+0.5))
     *
     * Controllable canonical form:
     * A = [0, 1, 0; 0, 0, 1; -1, -3.5, -3.5]
     * b = [0; 0; 1]
     * c = [3, 1, 0]
     */
    double A[] = { 0.0,  1.0,  0.0,
                    0.0,  0.0,  1.0,
                   -1.0, -3.5, -3.5};
    double b[] = {0.0, 0.0, 1.0};
    double c[] = {3.0, 1.0, 0.0};

    /* Sector: phi in [0, 10], time-invariant */
    double alpha = 0.0, beta = 10.0;

    AbsLureSystem *sys = abs_lure_create(A, b, c, 3, alpha, beta);
    if (!sys) {
        printf("Failed to create system.\n");
        return 1;
    }

    printf("System: G(s) = (s+3)/((s+1)(s+2)(s+0.5))\n");
    printf("Sector: phi in [0, %g] (time-invariant)\n\n", beta);

    /* ── Circle Criterion ── */
    printf("--- Circle Criterion ---\n");
    AbsCircleDisk disk = abs_circle_compute_disk(alpha, beta);
    printf("Disk: center = %g, radius = %g\n", disk.center_real, disk.radius);

    bool circle_pass = abs_circle_quick_check(sys);
    printf("Quick circle check: %s\n\n", circle_pass ? "PASS" : "FAIL");

    /* ── Popov Criterion ── */
    printf("--- Popov Criterion ---\n");
    AbsPopovResult *popov = abs_popov_test(sys, 150, 0.01, 100.0);
    if (!popov) {
        printf("Popov test failed.\n");
        abs_lure_free(sys);
        return 1;
    }

    abs_popov_print_result(popov);

    /* Show Popov plot points */
    printf("\nPopov plot (first 5 points):\n");
    for (int i = 0; i < 5 && i < popov->n_freqs; i++) {
        printf("  w=%.3f: X=%.4f, Y=%.4f\n",
               popov->freqs[i], popov->popov_X[i], popov->popov_Y[i]);
    }

    /* Compare with circle criterion conservatism */
    double circle_width = fabs(beta - alpha);
    double popov_max_k = abs_popov_max_stable_k(sys, 150, 0.01, 100.0);
    printf("\n--- Comparison ---\n");
    printf("Sector width tested:     %.2f\n", circle_width);
    printf("Max Popov-stable width:  %.2f\n", popov_max_k);

    /* Evaluate at specific q values */
    printf("\nPopov margin at various q:\n");
    double q_vals[] = {0.0, 0.5, 1.0, 2.0, 5.0};
    for (int i = 0; i < 5; i++) {
        double margin = abs_popov_eval_at_q(sys, q_vals[i], 150, 0.01, 100.0);
        printf("  q=%.1f: margin=%.4f\n", q_vals[i], margin);
    }

    abs_popov_result_free(popov);
    abs_lure_free(sys);

    printf("\nDone.\n");
    return 0;
}
