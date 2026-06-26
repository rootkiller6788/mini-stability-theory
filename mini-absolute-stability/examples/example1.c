/* example1.c — Circle Criterion demonstration
 *
 * Demonstrates the circle criterion on a 2nd-order Lur'e system
 * with sector [0, k] nonlinearity.
 *
 * System:  G(s) = 1 / ((s+1)(s+2))
 *          State-space: A = [-3, -2; 1, 0], b = [1; 0], c = [0, 1]
 *          Sector: phi in [0, 5]
 *
 * We compute the Nyquist plot and check the circle disk condition.
 */

#include "abs_core.h"
#include "abs_circle.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Example 1: Circle Criterion ===\n\n");

    /* Define a 2nd-order system with transfer function G(s) = 1/((s+1)(s+2)) */
    /* Controllable canonical form:
     *   ẋ₁ = -3x₁ - 2x₂ + u
     *   ẋ₂ =  x₁
     *   y  =  x₂
     */
    double A[] = {-3.0, -2.0,
                   1.0,  0.0};
    double b[] = {1.0, 0.0};
    double c[] = {0.0, 1.0};

    /* Sector: phi in [0, 5] */
    double alpha = 0.0, beta = 5.0;

    /* Create Lur'e system */
    AbsLureSystem *sys = abs_lure_create(A, b, c, 2, alpha, beta);
    if (!sys) {
        printf("Failed to create system.\n");
        return 1;
    }

    printf("System: G(s) = 1/((s+1)(s+2)), phi in [0, %g]\n", beta);

    /* Compute Nyquist disk */
    AbsCircleDisk disk = abs_circle_compute_disk(alpha, beta);
    printf("Nyquist disk: center = %g, radius = %g\n",
           disk.center_real, disk.radius);

    /* Run circle criterion test */
    AbsCircleResult *res = abs_circle_test(sys, 100, 0.01, 100.0);
    if (!res) {
        printf("Test failed.\n");
        abs_lure_free(sys);
        return 1;
    }

    abs_circle_print_result(res);

    /* Show sampled Nyquist values */
    printf("\nSampled Nyquist points (first 5):\n");
    for (int i = 0; i < 5 && i < res->n_freqs; i++) {
        printf("  w=%.3f: G(jw) = %+.4f %+.4fj\n",
               res->freqs[i], res->nyquist_re[i], res->nyquist_im[i]);
    }

    /* Also check the frequency inequality form */
    double freq_ineq = abs_circle_freq_inequality(sys, 100, 0.01, 100.0);
    printf("\nFrequency inequality min: %g (>0 => stable)\n", freq_ineq);

    /* Cleanup */
    abs_circle_result_free(res);
    abs_lure_free(sys);

    printf("\nDone.\n");
    return 0;
}
