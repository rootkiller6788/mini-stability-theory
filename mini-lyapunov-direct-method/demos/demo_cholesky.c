#include <stdio.h>
#include <assert.h>
#include "ldm_core.h"
#include "ldm_quadratic.h"
int main(void) {
    printf("=== Demo: Cholesky Definiteness ===\n");
    double P1[4]={2,0,0,3};
    printf("P=[2 0;0 3] PD=%s\n", ldm_is_positive_definite(P1,2)?"YES":"NO");
    double P2[4]={1,3,3,1};
    printf("P=[1 3;3 1] PD=%s (det=-8)\n", ldm_is_positive_definite(P2,2)?"YES":"NO");
    double P3[9]={4,1,0,1,4,1,0,1,4};
    printf("P=diag-3x3 PD=%s\n", ldm_is_positive_definite(P3,3)?"YES":"NO");
    return 0;
}
