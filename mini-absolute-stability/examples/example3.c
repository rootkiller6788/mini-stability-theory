/* example3.c — Aizerman Conjecture and KYP Lemma demonstration
 *
 * Demonstrates the Aizerman conjecture test and KYP lemma
 * on a Lur'e system. Shows the relationship between:
 *   - Hurwitz sector (linear stability)
 *   - Circle criterion (frequency-domain)
 *   - Popov criterion (frequency-domain, time-invariant)
 *   - KYP LMI (state-space)
 *
 * System: 2nd-order with transfer function:
 *   G(s) = (s+2) / (s^2 + 3s + 2)
 *
 * Aizerman holds for n=2, so the Hurwitz sector should
 * match the absolute stability sector.
 */

#include "abs_core.h"
#include "abs_lyapunov.h"
#include "abs_circle.h"
#include "abs_popov.h"
#include "abs_kyp.h"
#include "abs_aizerman.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

int main(void)
{
    printf("=== Example 3: Aizerman + KYP ===\n\n");

    /* 2nd-order system: G(s) = (s+2)/(s^2 + 3s + 2) = (s+2)/((s+1)(s+2)) = 1/(s+1)
     * Wait, that simplifies. Let's use a different system.
     *
     * G(s) = 2 / (s^2 + 2s + 2) with poles at -1±j
     * Controllable canonical form:
     *   A = [0, 1; -2, -2]
     *   b = [0; 1]
     *   c = [2, 0]
     */
    double A[] = { 0.0,  1.0,
                  -2.0, -2.0};
    double b[] = {0.0, 1.0};
    double c[] = {2.0, 0.0};

    double alpha = 0.0, beta = 3.0;

    AbsLureSystem *sys = abs_lure_create(A, b, c, 2, alpha, beta);
    if (!sys) {
        printf("Failed to create system.\n");
        return 1;
    }

    printf("System: G(s) = 2/(s^2 + 2s + 2)  (poles at -1±j)\n");
    printf("Sector: phi in [%g, %g]\n\n", alpha, beta);

    /* ── 1. Aizerman Conjecture Test ── */
    printf("=== 1. Aizerman Conjecture ===\n");
    AbsAizermanResult *aiz = abs_aizerman_test(sys);
    if (aiz) {
        abs_aizerman_print_result(aiz);
        printf("\n");
        abs_aizerman_result_free(aiz);
    }

    /* ── 2. Hurwitz Sector ── */
    printf("=== 2. Hurwitz Sector ===\n");
    double k_min, k_max;
    if (abs_aizerman_hurwitz_sector(A, b, c, 2, &k_min, &k_max)) {
        printf("Hurwitz range: [%.4f, %.4f]\n", k_min, k_max);
        printf("Sector [%.2f, %.2f] inside Hurwitz: %s\n\n",
               alpha, beta,
               (alpha >= k_min && beta <= k_max) ? "YES" : "NO");
    }

    /* ── 3. Lyapunov Equation ── */
    printf("=== 3. Lyapunov Equation ===\n");
    AbsLyapResult *lyap = abs_lyap_for_lure(A, 2);
    if (lyap && lyap->P) {
        printf("Lyapunov matrix P:\n");
        printf("  [ %.4f  %.4f ]\n", lyap->P[0], lyap->P[1]);
        printf("  [ %.4f  %.4f ]\n", lyap->P[2], lyap->P[3]);
        printf("Min eig P: %.4f, Max eig P: %.4f\n",
               lyap->min_eig_P, lyap->max_eig_P);
        printf("Separation: %.4f\n", abs_lyap_separation(A, 2));

        /* Verify stability */
        double margin = abs_lyap_verify_stability(A, b, c, 2,
                                                   lyap->P, alpha, beta, 50);
        printf("Stability margin: %.6f (negative => stable)\n\n", margin);
        abs_lyap_result_free(lyap);
    }

    /* ── 4. KYP Lemma ── */
    printf("=== 4. KYP Lemma ===\n");
    AbsKYPResult *kyp = abs_kyp_circle_to_lmi(sys);
    if (kyp && kyp->P) {
        printf("KYP solution found: P %s\n",
               kyp->is_feasible ? "FEASIBLE" : "INFEASIBLE");
        printf("Min eig P: %.4f\n", kyp->min_eig_P);

        bool verified = abs_kyp_verify_solution(sys, kyp->P, 1e-4);
        printf("LMI verified: %s\n", verified ? "YES" : "NO");
        abs_kyp_result_free(kyp);
    }

    /* Check SPR */
    bool is_spr = abs_kyp_is_spr(A, b, c, NULL, 2, 1, 1);
    printf("System is SPR: %s\n\n", is_spr ? "YES" : "NO");

    /* ── 5. Criteria Comparison ── */
    printf("=== 5. Criteria Comparison ===\n");

    /* Circle criterion */
    printf("Circle criterion:\n");
    AbsCircleResult *circle = abs_circle_test(sys, 50, 0.1, 50.0);
    if (circle) {
        printf("  Stable: %s, margin: %.4f at w=%.3f\n",
               circle->is_stable ? "YES" : "NO",
               circle->min_margin, circle->critical_freq);
        abs_circle_result_free(circle);
    }

    /* Popov criterion */
    printf("Popov criterion:\n");
    AbsPopovResult *popov = abs_popov_test(sys, 50, 0.1, 50.0);
    if (popov) {
        printf("  Status: %d, optimal q=%.4f, margin=%.4f\n",
               popov->status, popov->optimal_q, popov->min_margin);
        abs_popov_result_free(popov);
    }

    /* Conservatism gap */
    double circle_gap, popov_gap;
    abs_aizerman_conservatism_gap(sys, &circle_gap, &popov_gap);
    printf("\nConservatism gaps:\n");
    printf("  Circle gap: %.4f\n", circle_gap);
    printf("  Popov gap:  %.4f\n", popov_gap);

    abs_lure_free(sys);

    printf("\nDone.\n");
    return 0;
}
