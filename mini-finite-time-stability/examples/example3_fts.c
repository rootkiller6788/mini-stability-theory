#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fts_core.h"
#include "fts_sliding.h"
#include "fts_fixed_time.h"
#include "fts_convergence.h"

/* ================================================================
 * Example 3: Sliding Mode Control & Fixed-Time Stability
 *
 * Demonstrates:
 * 1. Terminal sliding surface design
 * 2. SMC simulation with double integrator
 * 3. Reaching time and sliding time computation
 * 4. Chattering analysis
 * 5. Fixed-time stability parameters
 * 6. Convergence classification
 * ================================================================ */

int main(void) {
    printf("=== Example 3: SMC & Fixed-Time Stability ===\n\n");

    /* --- Part 1: Terminal Sliding Mode Control --- */
    printf("--- Part 1: Terminal Sliding Mode Control ---\n");

    /* Create terminal sliding surface:
     * s(x) = x1 + x2 + beta * |x1|^alpha * sign(x1)
     * beta = 1.0, alpha = 0.5 */
    FTSSlidingSurface* sigma = fts_smc_surface_create(2, 1.0, 0.5, 2);
    if (!sigma) { printf("Error creating surface\n"); return 1; }

    printf("Terminal sliding surface:\n");
    printf("  s = x1 + x2 + %.1f * |x1|^%.1f * sign(x1)\n",
           sigma->gain, sigma->exponent);
    printf("  Terminal: %s\n\n", sigma->is_terminal ? "yes (0<alpha<1)" : "no");

    /* Evaluate surface at initial condition */
    FTSState x0; x0.dim = 2; x0.x[0] = 1.0; x0.x[1] = 0.0;
    double s0 = fts_smc_eval_surface(sigma, x0);
    printf("s(x0) = %.4f at x0 = (1.0, 0.0)\n", s0);

    /* Reaching time bound: T_r <= |s0| / K */
    double T_reach = fts_smc_reaching_time_bound(sigma, x0, 2.0);
    printf("Reaching time bound (K=2): T_r <= %.4f\n", T_reach);

    /* Sliding time bound: T_s <= |x1|^(1-alpha) / (beta*(1-alpha)) */
    double T_slide = fts_smc_sliding_time_bound(sigma, x0);
    printf("Sliding time bound: T_s <= %.4f\n\n", T_slide);

    /* --- SMC Simulation --- */
    printf("--- SMC Simulation (Double Integrator) ---\n");

    FTSTerminalSMC* smc = fts_smc_create(sigma, 2);
    if (!smc) { printf("Error creating SMC\n"); fts_smc_surface_free(sigma); return 1; }

    FTSSystem* di = fts_create_double_integrator(0.001);
    if (!di) { printf("Error creating system\n"); fts_smc_free(smc); fts_smc_surface_free(sigma); return 1; }

    /* Simulate */
    int steps = fts_smc_simulate(smc, di, x0, 10.0, 0.001);
    printf("Simulation steps: %d\n", steps);

    fts_smc_print(smc);

    /* Control effort and chattering */
    double effort = fts_smc_control_effort(smc);
    double chatter = fts_smc_chattering_measure(smc);
    printf("Control effort: %.4f\n", effort);
    printf("Chattering measure: %.4f (%.1f%% sign changes)\n",
           chatter, chatter * 100.0);

    fts_smc_surface_free(sigma);
    fts_smc_free(smc);
    fts_free(di);

    /* --- Part 2: Fixed-Time Stability --- */
    printf("\n--- Part 2: Fixed-Time Stability (Polyakov) ---\n");

    /* Fixed-time parameters:
     * alpha=0.5, beta=2.0, p=2.0, q=2.0
     * T_max = 1/(2*(1-0.5)) + 1/(2*(2-1)) = 1/1 + 1/2 = 1.5 */
    double alpha = 0.5, beta = 2.0, p = 2.0, q = 2.0;
    FTSFixedTimeParams* ftp = fts_ft_params_create(alpha, beta, p, q);
    if (!ftp) { printf("Error creating params\n"); return 1; }

    fts_ft_params_print(ftp);

    double T_max = fts_ft_settling_bound(ftp);
    printf("\nFixed-time property:\n");
    printf("  T_max = 1/(%.1f*(1-%.1f)) + 1/(%.1f*(%.1f-1)) = %.4f\n",
           p, alpha, q, beta, T_max);
    printf("  This bound is UNIFORM - independent of initial condition!\n");

    /* Compare with finite-time */
    printf("\nComparison (for x0 = 10.0):\n");
    FTSState x_big; x_big.dim = 1; x_big.x[0] = 10.0;
    double T_finite = pow(x_big.x[0], 1.0-alpha) / (p * (1.0 - alpha));
    printf("  Finite-time T(x0=10):  %.4f (grows with ||x0||)\n", T_finite);
    printf("  Fixed-time T_max:      %.4f (constant!)\n", T_max);
    printf("  Fixed-time is better:  %s\n",
           (T_max < T_finite) ? "YES" : "NO");

    /* --- Part 3: Convergence Classification --- */
    printf("\n--- Part 3: Convergence Classification ---\n");

    FTSSystem* sfts = fts_create_scalar_fts(0.5, 0.001);
    if (!sfts) { fts_ft_params_free(ftp); return 1; }
    sfts->state.x[0] = 1.0;

    FTSSettlingTime* st = fts_st_create(10);
    if (st) {
        fts_st_compute(sfts, x_big, 10.0, 0.001, 0.01, st);
        fts_st_print(st);

        FTSConvergenceType ct = fts_st_classify(st);
        printf("Classification: %s\n",
               ct == FTS_FINITE_TIME ? "Finite-Time" :
               ct == FTS_FIXED_TIME  ? "Fixed-Time"  :
               ct == FTS_EXPONENTIAL ? "Exponential" :
               ct == FTS_ASYMPTOTIC  ? "Asymptotic"  : "Practical");
        fts_st_free(st);
    }

    fts_ft_params_free(ftp);
    fts_free(sfts);
    printf("\nExample 3 completed.\n");
    return 0;
}
