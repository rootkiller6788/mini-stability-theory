#include "ios_core.h"
#include <stdio.h>

/* Example 3: Passivity-Based Stability.
 * Demonstrates passivity verification and the passivity theorem:
 * Negative feedback interconnection of passive systems is stable.
 * Shows checking passivity, strict passivity, and feedback stability. */
int main(void) {
    printf("=== Example 3: Passivity-Based Stability ===\n\n");

    /* Passive system: G(s) = 1/(s+1) */
    printf("Passive System: G(s) = 1/(s+1)\n");
    IOSTransferFcn* tf1 = ios_tf_create(1, 2);
    if (tf1) {
        tf1->num[0] = 1.0;
        tf1->den[0] = 1.0; tf1->den[1] = 1.0;
        IOSPassivityResult r = ios_check_passivity_tf(tf1);
        printf("  Passive:          %s\n", r.passive ? "YES" : "NO");
        printf("  Strictly Input:   %s\n", r.strictly_input ? "YES" : "NO");
        printf("  Strictly Output:  %s\n", r.strictly_output ? "YES" : "NO");
        printf("  Dissipation:      %.6f\n", r.dissipation);
        ios_tf_free(tf1);
    }

    /* Non-passive system: G(s) = -1/(s+1) */
    printf("\nNon-Passive System: G(s) = -1/(s+1)\n");
    IOSTransferFcn* tf2 = ios_tf_create(1, 2);
    if (tf2) {
        tf2->num[0] = -1.0;
        tf2->den[0] = 1.0; tf2->den[1] = 1.0;
        IOSPassivityResult r = ios_check_passivity_tf(tf2);
        printf("  Passive:          %s\n", r.passive ? "YES" : "NO");
        ios_tf_free(tf2);
    }

    /* Feedback of two passive state-space systems */
    printf("\nFeedback Passivity Test:\n");
    IOSStateSpace* G = ios_ss_create(1, 1, 1);
    IOSStateSpace* H = ios_ss_create(1, 1, 1);
    if (G && H) {
        G->A[0] = -2.0; G->B[0] = 1.0; G->C[0] = 1.0;
        H->A[0] = -5.0; H->B[0] = 1.0; H->C[0] = 1.0;
        printf("  G passive: %s\n",
               ios_check_passivity_ss(G).passive ? "YES" : "NO");
        printf("  H passive: %s\n",
               ios_check_passivity_ss(H).passive ? "YES" : "NO");
        printf("  Feedback stable: %s\n",
               ios_feedback_passivity_stable(G, H) ? "YES" : "INCONCLUSIVE");
    }
    ios_ss_free(G); ios_ss_free(H);

    return 0;
}
