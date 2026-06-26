#include <stdio.h>
#include "ldm_core.h"
#include "ldm_linear.h"
#include "ldm_quadratic.h"
int main(void){printf("=== Linear Stability via Lyapunov ===\n\n");
double A[4]={-1,2,-3,-4};LDM_LinearAnalysis* la=ldm_linear_analyze(A,2);
printf("Stable: %s, eigenvalues: %.3f,%.3f\n",la->is_stable?"YES":"NO",la->eigenvalues[0],la->eigenvalues[1]);
double Q[4]={1,0,0,1};LDM_QuadraticForm* qf=ldm_lyap_solve_2x2(A,Q);
printf("P=[%.4f %.4f; %.4f %.4f], PD=%s\n",qf->P[0],qf->P[1],qf->P[2],qf->P[3],qf->is_pd?"YES":"NO");
printf("Convergence: %.4f, Settling(2%%): %.3fs\n",ldm_quadratic_convergence_rate(A,qf->P,Q,2),ldm_linear_settling_time_estimate(A,2,0.02));
ldm_linear_analysis_free(la);ldm_quadratic_free(qf);
printf("Key: P>0 from Lyapunov eq proves asymptotic stability.\n");return 0;}
