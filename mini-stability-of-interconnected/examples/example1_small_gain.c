#include "intercon_core.h"
#include "intercon_small_gain.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=======================================\n");
    printf("  Interconnected Stability - Demo 1\n");
    printf("  Small-Gain & Circle Criterion\n");
    printf("=======================================\n\n");

    /* 1. Create two stable subsystems */
    printf("--- 1. Small-Gain Theorem ---\n");
    ICSubsystem* G1 = ic_subsys_create(1, 1, 1, false);
    double A1[] = {-2.0}; ic_subsys_set_A(G1, A1);
    double B1[] = {1.0}; ic_subsys_set_B(G1, B1);
    double C1[] = {0.5}; ic_subsys_set_C(G1, C1);
    double D1[] = {0.0}; ic_subsys_set_D(G1, D1);
    printf("G1: 1st-order, pole=-2, DCgain=0.25\n");

    ICSubsystem* G2 = ic_subsys_create(1, 1, 1, false);
    double A2[] = {-1.0}; ic_subsys_set_A(G2, A2);
    double B2[] = {1.0}; ic_subsys_set_B(G2, B2);
    double C2[] = {1.0}; ic_subsys_set_C(G2, C2);
    double D2[] = {0.0}; ic_subsys_set_D(G2, D2);
    printf("G2: 1st-order, pole=-1, DCgain=1.0\n\n");

    SmallGainResult* sg = ic_small_gain_test(G1, G2, IC_FEEDBACK_NEG);
    ic_small_gain_result_print(sg);
    ic_small_gain_result_free(sg);

    /* 2. Vary feedback gain */
    printf("\n--- 2. Feedback Gain Variation ---\n");
    for (double k = 0.2; k <= 5.0; k *= 2.0) {
        ICSystem* sys = ic_system_create();
        ICSubsystem* s1 = ic_subsys_create(1,1,1,false);
        ic_subsys_set_A(s1, A1); ic_subsys_set_B(s1, B1);
        ic_subsys_set_C(s1, C1); ic_subsys_set_D(s1, D1);
        ic_system_add_subsys(sys, s1);
        ICSubsystem* s2 = ic_subsys_create(1,1,1,false);
        ic_subsys_set_A(s2, A2); ic_subsys_set_B(s2, B2);
        ic_subsys_set_C(s2, C2); ic_subsys_set_D(s2, D2);
        ic_system_add_subsys(sys, s2);
        ic_system_connect(sys, IC_FEEDBACK_NEG, 0, 1, k);
        SmallGainResult* sgk = ic_small_gain_test(s1, s2, IC_FEEDBACK_NEG);
        printf("  k=%.1f  product=%.3f  %s\n",
               k, sgk->product, sgk->is_stable ? "STABLE" : "MAY BE UNSTABLE");
        ic_small_gain_result_free(sgk);
        ic_system_free(sys);
    }

    /* 3. Circle criterion */
    printf("\n--- 3. Circle Criterion ---\n");
    CircleCriterion* cc = ic_circle_criterion(G1, 0.1, 1.0);
    ic_circle_criterion_print(cc);
    ic_circle_criterion_free(cc);

    ic_subsys_free(G1);
    ic_subsys_free(G2);

    printf("\n=======================================\n");
    return 0;
}