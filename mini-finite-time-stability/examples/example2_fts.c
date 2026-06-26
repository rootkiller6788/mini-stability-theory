#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fts_core.h"
#include "fts_homogeneous.h"

/* ================================================================
 * Example 2: Homogeneity-Based Finite-Time Stability Detection
 *
 * Demonstrates:
 * 1. Homogeneous vector field definition
 * 2. Dilation operator
 * 3. Homogeneity verification
 * 4. Negative degree FTS detection
 * 5. Convergence testing with homogeneous norms
 * ================================================================ */

int main(void) {
    printf("=== Example 2: Homogeneity-Based FTS Detection ===\n\n");

    /* Set up scalar FTS system */
    FTSSystem* sys = fts_create_scalar_fts(0.5, 0.001);
    if (!sys) { printf("Error creating system\n"); return 1; }

    double x_init[] = {1.0};
    fts_set_state(sys, x_init);

    printf("System: dx/dt = -|x|^0.5 * sign(x)\n");

    /* --- Homogeneity Analysis --- */
    printf("\n--- Homogeneous Weights ---\n");

    /* For scalar system with dx/dt = -|x|^alpha*sign(x):
     * f(lambda^r*x) = -|lambda^r*x|^alpha*sign(lambda^r*x)
     *               = -lambda^{r*alpha}*|x|^alpha*sign(x)
     *               = lambda^{r*alpha} * f(x)
     *
     * For homogeneity of degree d: f(lambda^r*x) = lambda^{d+r}*f(x)
     * So: r*alpha = d + r => d = r*(alpha - 1)
     *
     * With r=1, alpha=0.5: d = 1*(0.5-1) = -0.5
     * Negative degree => FTS! */

    double w[] = {1.0};
    FTSHomogeneousWeights* hw = fts_hom_weights_create(1, w);
    if (!hw) { printf("Error creating weights\n"); fts_free(sys); return 1; }

    double deg = -0.5;
    hw->degree = deg;
    hw->is_negative = true;
    fts_hom_print(hw);

    /* --- FTS Detection --- */
    printf("\n--- FTS Detection ---\n");
    bool is_fts = fts_hom_is_finite_time(hw, deg);
    printf("Negative degree (%.1f < 0): %s\n", deg, is_fts ? "YES -> FTS" : "NO");

    /* --- Dilation --- */
    printf("\n--- Dilation Operator ---\n");
    FTSState x; x.dim = 1; x.x[0] = 3.0;
    printf("x = (%.2f)\n", x.x[0]);
    for (int lam = 1; lam <= 4; lam++) {
        FTSState x_dil = fts_hom_dilate(x, hw, (double)lam);
        printf("Delta_%d(x) = (%.2f)\n", lam, x_dil.x[0]);
    }

    /* Dilations form a group:
     * Delta_2(Delta_3(x)) = Delta_6(x) */
    FTSState x_23 = fts_hom_dilate(fts_hom_dilate(x, hw, 3.0), hw, 2.0);
    FTSState x_6 = fts_hom_dilate(x, hw, 6.0);
    printf("\nGroup property: Delta_2(Delta_3(x)) = %.2f, Delta_6(x) = %.2f\n",
           x_23.x[0], x_6.x[0]);

    /* --- Homogeneous Norm --- */
    printf("\n--- Homogeneous Norm ---\n");
    for (double v = 0.5; v <= 3.0; v += 0.5) {
        FTSState xv; xv.dim = 1; xv.x[0] = v;
        double nrm = fts_hom_norm(xv, hw);
        printf("||(%.1f)||_r = %.4f\n", v, nrm);
    }

    /* --- Convergence Test --- */
    printf("\n--- Convergence Test ---\n");
    FTSState x0; x0.dim = 1; x0.x[0] = 1.0;
    FTSHomogeneousTest* ht = fts_hom_test_convergence(sys, x0, hw, 10.0, 0.001);
    if (ht) {
        fts_hom_test_print(ht);
        fts_hom_test_free(ht);
    }

    /* Settling time bound via homogeneity */
    double T_hom = fts_hom_settling_time_bound(hw, deg, x0);
    printf("Homogeneity settling time bound: T <= %.4f (K=1)\n", T_hom);

    /* Standard weights for double integrator */
    printf("\n--- Standard Weights (Double Integrator) ---\n");
    FTSHomogeneousWeights* sw = fts_hom_standard_weights(2);
    if (sw) {
        fts_hom_print(sw);
        printf("FTS with deg=-1: %s\n",
               fts_hom_is_finite_time(sw, -1.0) ? "YES" : "NO");
        fts_hom_weights_free(sw);
    }

    fts_hom_weights_free(hw);
    fts_free(sys);
    printf("\nExample 2 completed.\n");
    return 0;
}
