#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ldm_core.h"
#include "ldm_linear.h"
#include "ldm_nonlinear.h"
int main(void) {
    printf("=== Demo: ROA via Sampling ===\n");
    double A[4]={-1,0,0,-2};
    LDM_System* s=ldm_system_create_linear(2,A);
    LDM_RegionEstimate re=ldm_roa_sampling(s,2,2000,5.0,10.0);
    printf("Linear system A=[-1 0;0 -2]\n");
    printf("Sampled RoA volume: %.4f (range=5)\n",re.volume_estimate);
    printf("Point in RoA (0.5,0.5): %s\n",
        ldm_point_in_roa(s,(double[]){0.5,0.5},2,10.0,0.01)?"YES":"NO");
    printf("Point in RoA (10,10): %s\n",
        ldm_point_in_roa(s,(double[]){10,10},2,0.5,0.01)?"YES":"NO");
    ldm_region_free(&re);
    ldm_system_free(s);
    return 0;
}
