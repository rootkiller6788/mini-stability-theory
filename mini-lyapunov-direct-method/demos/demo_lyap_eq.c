#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "ldm_core.h"
#include "ldm_quadratic.h"
#include "ldm_linear.h"
int main(void) {
    printf("=== Demo: Lyapunov Equation Solver ===\n");
    double A[4]={-2,1,0,-3};
    double Q[4]={1,0,0,1};
    LDM_QuadraticForm* qf=ldm_lyap_solve_2x2(A,Q);
    printf("A=[-2 1;0 -3]  Q=I\n");
    if(qf&&qf->is_pd) {
        printf("P=[%.4f %.4f; %.4f %.4f]  PD=YES\n",
            qf->P[0],qf->P[1],qf->P[2],qf->P[3]);
        printf("lambda(P)=[%.4f,%.4f]  cond=%.2f\n",
            qf->min_eigenvalue,qf->max_eigenvalue,qf->condition_number);
        printf("residual=%.2e\n",ldm_lyapunov_residual(A,qf->P,Q,2));
        printf("convergence=%.4f\n",ldm_quadratic_convergence_rate(A,qf->P,Q,2));
    }
    ldm_quadratic_free(qf);
    return 0;
}
