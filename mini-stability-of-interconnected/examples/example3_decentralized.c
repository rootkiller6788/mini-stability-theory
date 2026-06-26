#include "intercon_core.h"
#include "intercon_decentralized.h"
#include "intercon_lyapunov.h"
#include <stdio.h>

int main(void) {
    printf("=======================================\n");
    printf("  Interconnected Stability - Demo 3\n");
    printf("  Decentralized Control & Interaction\n");
    printf("=======================================\n\n");

    /* Build multi-subsystem plant */
    ICSystem* sys = ic_system_create();

    /* Subsystem 1 */
    ICSubsystem* s1 = ic_subsys_create(2, 1, 1, false);
    double A1[] = {-1.0, 0.2, 0.1, -2.0};
    ic_subsys_set_A(s1, A1);
    double B1[] = {1.0, 0.0};
    ic_subsys_set_B(s1, B1);
    double C1[] = {1.0, 0.0};
    ic_subsys_set_C(s1, C1);
    ic_system_add_subsys(sys, s1);

    /* Subsystem 2 */
    ICSubsystem* s2 = ic_subsys_create(1, 1, 1, false);
    double A2[] = {-3.0};
    ic_subsys_set_A(s2, A2);
    double B2[] = {1.0};
    ic_subsys_set_B(s2, B2);
    double C2[] = {1.0};
    ic_subsys_set_C(s2, C2);
    ic_system_add_subsys(sys, s2);

    /* Interconnect */
    ic_system_connect(sys, IC_FEEDBACK_NEG, 0, 1, 0.5);
    ic_system_connect(sys, IC_FEEDBACK_NEG, 1, 0, 0.2);

    printf("--- 1. Interaction Analysis ---\n");
    double ib = ic_interaction_norm_bound(sys);
    printf("Interaction bound: %.4f\n", ib);

    DecentralizedStability* ds = ic_decentralized_stability_test(sys);
    ic_decentralized_stability_print(ds);
    ic_decentralized_stability_free(ds);

    printf("\n--- 2. Controllability Check ---\n");
    printf("Decentralized controllable: %s\n",
           ic_is_decentralized_controllable(sys) ? "YES" : "NO");
    printf("Decentralized observable: %s\n",
           ic_is_decentralized_observable(sys) ? "YES" : "NO");
    printf("Fixed modes: %d\n", ic_fixed_modes_count(sys));

    printf("\n--- 3. Sparse H-infinity Design ---\n");
    SparseHinfController* sc = ic_sparse_hinf_design(sys, 2.0);
    ic_sparse_hinf_print(sc);
    ic_sparse_hinf_free(sc);

    printf("\n--- 4. Overlapping Decomposition ---\n");
    double A_comp[IC_MAX_STATES*IC_MAX_STATES];
    ic_build_composite_A(sys, A_comp, IC_MAX_STATES);
    OverlappingDecomp* od = ic_overlapping_decompose(A_comp, sys->total_states);
    ic_overlapping_decomp_print(od);
    ic_overlapping_decomp_free(od);

    printf("\n--- 5. Connective Stability ---\n");
    printf("Connectively stable: %s\n",
           ic_connective_stability_test(sys) ? "YES" : "NO");

    ic_system_free(sys);
    printf("\n=======================================\n");
    return 0;
}