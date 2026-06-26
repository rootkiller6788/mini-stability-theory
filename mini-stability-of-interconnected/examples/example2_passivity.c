#include "intercon_core.h"
#include "intercon_passivity.h"
#include "intercon_lyapunov.h"
#include <stdio.h>

int main(void) {
    printf("=======================================\n");
    printf("  Interconnected Stability - Demo 2\n");
    printf("  Passivity & Lyapunov Methods\n");
    printf("=======================================\n\n");

    /* 1. Passivity test */
    printf("--- 1. Passivity Test ---\n");
    ICSubsystem* P1 = ic_subsys_create(1, 1, 1, false);
    double A[] = {-1.0}; ic_subsys_set_A(P1, A);
    double B[] = {1.0}; ic_subsys_set_B(P1, B);
    double C[] = {1.0}; ic_subsys_set_C(P1, C);
    double D[] = {0.5}; ic_subsys_set_D(P1, D);
    PassivityResult* pr = ic_passivity_test(P1);
    ic_passivity_result_print(pr);
    ic_passivity_result_free(pr);

    /* 2. Build interconnected system */
    printf("\n--- 2. Composite Lyapunov ---\n");
    ICSystem* sys = ic_system_create();
    ICSubsystem* s1 = ic_subsys_create(1,1,1,false);
    double A1[] = {-2.0}; ic_subsys_set_A(s1, A1);
    double B1[] = {1.0}; ic_subsys_set_B(s1, B1);
    double C1[] = {1.0}; ic_subsys_set_C(s1, C1);
    double D1[] = {0.0}; ic_subsys_set_D(s1, D1);
    ic_system_add_subsys(sys, s1);

    ICSubsystem* s2 = ic_subsys_create(1,1,1,false);
    double A2[] = {-3.0}; ic_subsys_set_A(s2, A2);
    double B2[] = {1.0}; ic_subsys_set_B(s2, B2);
    double C2[] = {1.0}; ic_subsys_set_C(s2, C2);
    double D2[] = {0.0}; ic_subsys_set_D(s2, D2);
    ic_system_add_subsys(sys, s2);

    ic_system_connect(sys, IC_FEEDBACK_NEG, 0, 1, 0.3);
    double w[] = {1.0, 1.0};
    CompositeLyapunov* cl = ic_composite_lyapunov(sys, w);
    ic_composite_lyapunov_print(cl);
    ic_composite_lyapunov_free(cl);

    /* 3. M-matrix */
    printf("\n--- 3. M-Matrix Analysis ---\n");
    MMatrixResult* mr = ic_build_M_matrix(sys);
    ic_M_matrix_print(mr);
    printf("M-matrix stability: %s\n",
           ic_M_matrix_stability(mr) ? "STABLE" : "UNSTABLE");
    ic_M_matrix_free(mr);

    /* 4. Gershgorin */
    printf("\n--- 4. Gershgorin Test ---\n");
    double A_comp[IC_MAX_STATES * IC_MAX_STATES];
    ic_build_composite_A(sys, A_comp, IC_MAX_STATES);
    GershgorinResult* gr = ic_gershgorin_test(A_comp, 2);
    ic_gershgorin_print(gr);
    printf("Gershgorin stable: %s\n",
           ic_gershgorin_stable(gr) ? "YES" : "NO");
    ic_gershgorin_free(gr);

    /* 5. Stability margins */
    printf("\n--- 5. Stability Margins ---\n");
    double gm = ic_gain_margin(s1, s2);
    double pm = ic_phase_margin(s1, s2);
    printf("Gain margin: %.2f  Phase margin: %.1f deg\n", gm, pm);

    ic_system_free(sys);
    ic_subsys_free(P1);

    printf("\n=======================================\n");
    return 0;
}