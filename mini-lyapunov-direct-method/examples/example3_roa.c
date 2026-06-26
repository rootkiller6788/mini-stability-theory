#include <stdio.h>
#include "ldm_core.h"
#include "ldm_linear.h"
#include "ldm_region.h"
#include <stdlib.h>
int main(void){printf("=== Region of Attraction ===\n\n");
double A[4]={-1,0,0,-2};LDM_System* s=ldm_system_create_linear(2,A);
double P[4]={1,0,0,1};LDM_RegionEstimate re=ldm_roa_sublevel(s,P,2,5.0);
ldm_region_print(&re);
LDM_RegionEstimate re2=ldm_roa_sampling(s,2,500,3.0,5.0);
printf("Sampling volume: %.4f\n",re2.volume_estimate);
printf("x0=(0.5,0.5): %s\n",ldm_point_in_roa(s,(double[]){0.5,0.5},2,5.0,0.01)?"IN":"OUT");
ldm_system_free(s);ldm_region_free(&re);ldm_region_free(&re2);
printf("Key: Sublevel sets give conservative RoA estimates.\n");return 0;}
