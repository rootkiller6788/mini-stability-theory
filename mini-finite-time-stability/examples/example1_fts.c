#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fts_core.h"
#include "fts_lyapunov.h"
#include "fts_convergence.h"

/* ================================================================
 * Example 1: Finite-Time Convergence via Lyapunov Analysis
 *
 * Demonstrates:
 * 1. Scalar FTS system: dx/dt = -|x|^alpha * sign(x) with alpha=0.5
 * 2. Lyapunov function V(x) = 0.5*x^2
 * 3. Settling time bound: T <= V0^(1-alpha)/(c*(1-alpha))
 * 4. Numerical convergence verification
 * ================================================================ */

int main(void) {
    printf("=== Example 1: Finite-Time Stability via Lyapunov ===\n\n");

    /* Create scalar FTS system with alpha=1/2 (0.5) */
    double alpha = 0.5;
    double dt = 0.001;
    FTSSystem* sys = fts_create_scalar_fts(alpha, dt);
    if (!sys) { printf("Error creating system\n"); return 1; }

    /* Set initial condition x0 = 1.0 */
    double x_init[] = {1.0};
    fts_set_state(sys, x_init);

    printf("System: dx/dt = -|x|^%.1f * sign(x), x0 = %.2f\n", alpha, x_init[0]);
    printf("dt = %.4f\n\n", dt);

    /* Analytical settling time:
     * For dx/dt = -|x|^alpha * sign(x), the solution is:
     *   x(t) = sign(x0)*(|x0|^{1-alpha} - (1-alpha)*t)^{1/(1-alpha)}
     *   for t < |x0|^{1-alpha}/(1-alpha), and x(t)=0 thereafter.
     *
     * T_settling = |x0|^{1-alpha} / (1-alpha) = 1^{0.5}/0.5 = 2.0 */
    double T_analytical = pow(fabs(x_init[0]), 1.0 - alpha) / (1.0 - alpha);
    printf("Analytical settling time: T = |x0|^(1-%.1f)/(1-%.1f) = %.4f\n",
           alpha, alpha, T_analytical);

    /* Lyapunov analysis */
    FTSLyapunovFunc* lf = fts_lyap_create(fts_lyap_default, NULL,
                                           0.5, 1.0, 0.0);
    if (!lf) { printf("Error creating Lyapunov function\n"); fts_free(sys); return 1; }

    /* Lyapunov-based settling time bound:
     * V(x) = 0.5*x^2, V0 = 0.5*1^2 = 0.5
     * c = 1, alpha = 0.5
     * T_bound = V0^(1-alpha) / (c*(1-alpha)) = 0.5^0.5 / (1*0.5) = 0.7071/0.5 = 1.414 */
    FTSState x0; x0.dim = 1; x0.x[0] = 1.0;
    double T_bound = fts_lyap_settling_time_bound(lf, x0);
    printf("Lyapunov settling time bound: T <= %.4f\n", T_bound);

    /* Full Lyapunov analysis */
    printf("\n--- Full Lyapunov Analysis ---\n");
    FTSLyapunovAnalysis* la = fts_lyap_analyze(sys, x0, 10.0, dt, 100);
    if (la) {
        fts_lyap_print(la);
        printf("FTS verified: %s\n",
               fts_lyap_is_finite_time_stable(la) ? "YES" : "NO");
        fts_lyap_analysis_free(la);
    }

    /* Empirical settling time */
    printf("\n--- Empirical Settling Time ---\n");
    /* Reset system and simulate */
    fts_set_state(sys, x_init);
    double T_emp = fts_st_empirical(sys, x0, dt, 1e-3);
    printf("Empirical settling time (tol=1e-3): %.6f\n", T_emp);

    printf("\nComparison:\n");
    printf("  Analytical:          %.4f\n", T_analytical);
    printf("  Lyapunov bound:      %.4f\n", T_bound);
    printf("  Empirical:           %.4f\n", T_emp);

    fts_lyap_free(lf);
    fts_free(sys);
    printf("\nExample 1 completed.\n");
    return 0;
}
