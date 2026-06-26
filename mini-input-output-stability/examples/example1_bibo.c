#include "ios_core.h"
#include <stdio.h>

/* Example 1: BIBO stability analysis.
 * Demonstrates checking BIBO stability for a first-order transfer function
 * and computing its H-infinity norm. */
int main(void) {
    printf("=== Example 1: BIBO Stability ===\n\n");
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return 1;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;

    printf("Transfer function: G(s) = 1 / (s + 1)\n");
    printf("  BIBO status:    %s\n",
           ios_check_bibo_tf(tf) == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE");
    printf("  Hinf norm:      %.4f\n", ios_Hinf_norm_tf(tf));
    printf("  H2 norm:        %.4f\n", ios_H2_norm_tf(tf));
    printf("  Bandwidth:      %.4f rad/s\n", ios_bandwidth(tf));

    /* Demonstrate unstable case */
    IOSTransferFcn* tf_unstable = ios_tf_create(1, 2);
    if (tf_unstable) {
        tf_unstable->num[0] = 1.0;
        tf_unstable->den[0] = 1.0; tf_unstable->den[1] = -1.0;
        printf("\nUnstable: G(s) = 1 / (s - 1)\n");
        printf("  BIBO status:    %s\n",
               ios_check_bibo_tf(tf_unstable) == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE");
        ios_tf_free(tf_unstable);
    }
    ios_tf_free(tf);
    return 0;
}
