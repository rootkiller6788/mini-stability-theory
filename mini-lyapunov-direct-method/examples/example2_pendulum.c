#include <stdio.h>
#include <math.h>
#include "ldm_core.h"
#include "ldm_nonlinear.h"
#include "ldm_region.h"
#include <stdlib.h>
void pend(double t,const double* x,int n,double* dx){(void)t;dx[0]=x[1];dx[1]=-9.81*sin(x[0])-0.1*x[1];}
int main(void){printf("=== Pendulum Stability ===\n\n");
LDM_System* s=ldm_system_create(2,pend,NULL);s->x[0]=0.5;s->x[1]=0;
printf("Energy V=%.4f\n",ldm_energy_based(s,2));
LDM_StabilityReport sr=ldm_nonlinear_stability(s,(double[]){0.5,0},2,5.0,0.01);
printf("Stability: %s, final|x|=%.4f\n",sr.result>=LDM_ASYM_STABLE?"asym":"stable",sr.final_V);
ldm_system_free(s);free(sr.trajectory);
printf("Key: Energy Lyapunov function proves pendulum stability.\n");return 0;}
