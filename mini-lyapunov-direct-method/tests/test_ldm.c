#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ldm_core.h"
#include "ldm_quadratic.h"
#include "ldm_linear.h"
#include "ldm_nonlinear.h"
#include "ldm_region.h"
#include "ldm_converse.h"

#define EPS 1e-6
#define NEAR(a,b,e) assert(fabs((a)-(b))<=(e))

int main(void) {
    int test_count=0;
    printf("\n=== Lyapunov Direct Method Tests ===\n\n");

    /* 1. System create/free */
    printf("  [%d] system_create... ", ++test_count); fflush(stdout);
    { LDM_System* s=ldm_system_create(2,NULL,NULL);
      assert(s); assert(s->n==2); ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 2. Linear system create */
    printf("  [%d] system_create_linear... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}; LDM_System* s=ldm_system_create_linear(2,A);
      assert(s); assert(s->is_linear); assert(s->A);
      assert(s->A[0]==-1); assert(s->A[3]==-2); ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 3. Quadratic form */
    printf("  [%d] quadratic_form... ", ++test_count); fflush(stdout);
    { double P[4]={1,0,0,1}, x[2]={3,4};
      NEAR(ldm_quadratic_form(P,x,2),25.0,EPS); }
    printf("PASSED\n"); fflush(stdout);

    /* 4. Positive definite (Cholesky for n>2) */
    printf("  [%d] positive_definite... ", ++test_count); fflush(stdout);
    { assert(ldm_is_positive_definite((double[]){2,0,0,2},2));
      assert(ldm_is_positive_definite((double[]){3,1,0,1,3,0,0,0,3},3));
      assert(!ldm_is_positive_definite((double[]){-1,0,0,-1},2));
      assert(!ldm_is_positive_definite((double[]){0,0,0,0},2)); }
    printf("PASSED\n"); fflush(stdout);

    /* 5. Positive semidefinite */
    printf("  [%d] positive_semidefinite... ", ++test_count); fflush(stdout);
    { assert(ldm_is_positive_semidefinite((double[]){1,0,0,0},2)); }
    printf("PASSED\n"); fflush(stdout);

    /* 6. Lyap solve 2x2 */
    printf("  [%d] lyap_solve_2x2... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}, Q[4]={1,0,0,1};
      LDM_QuadraticForm* qf=ldm_lyap_solve(A,Q,2);
      assert(qf); assert(qf->is_pd);
      NEAR(qf->P[0],0.5,EPS); NEAR(qf->P[3],0.25,EPS);
      ldm_quadratic_free(qf); }
    printf("PASSED\n"); fflush(stdout);

    /* 7. Lyap solve 3x3 */
    printf("  [%d] lyap_solve_3x3... ", ++test_count); fflush(stdout);
    { double A[9]={-1,0,0,0,-2,0,0,0,-3}, Q[9]={1,0,0,0,1,0,0,0,1};
      LDM_QuadraticForm* qf=ldm_lyap_solve(A,Q,3);
      assert(qf); assert(qf->is_pd); ldm_quadratic_free(qf); }
    printf("PASSED\n"); fflush(stdout);

    /* 8. Linear system stability (Hurwitz) */
    printf("  [%d] linear_is_hurwitz... ", ++test_count); fflush(stdout);
    { assert(ldm_linear_is_hurwitz((double[]){-1,0,0,-2},2));
      assert(!ldm_linear_is_hurwitz((double[]){1,0,0,1},2)); }
    printf("PASSED\n"); fflush(stdout);

    /* 9. Spectral abscissa */
    printf("  [%d] spectral_abscissa... ", ++test_count); fflush(stdout);
    { double sa=ldm_linear_spectral_abscissa((double[]){-1,0,0,-2},2);
      assert(sa<-0.9); }
    printf("PASSED\n"); fflush(stdout);

    /* 10. Quadratic is Lyapunov */
    printf("  [%d] quadratic_is_lyapunov... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}, P[4]={1,0,0,1};
      assert(ldm_quadratic_is_lyapunov(A,P,2)); }
    printf("PASSED\n"); fflush(stdout);

    /* 11. Quadratic convergence rate */
    printf("  [%d] convergence_rate... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}, Q[4]={1,0,0,1};
      LDM_QuadraticForm* qf=ldm_lyap_solve(A,Q,2);
      assert(qf); double cr=ldm_quadratic_convergence_rate(A,qf->P,Q,2);
      assert(cr>0); ldm_quadratic_free(qf); }
    printf("PASSED\n"); fflush(stdout);

    /* 12. System step (RK4) */
    printf("  [%d] system_step... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-1}; LDM_System* s=ldm_system_create_linear(2,A);
      s->x[0]=1; s->x[1]=1; ldm_system_step(s,0.1);
      assert(s->x[0]<1.0); assert(s->x[1]<1.0); ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 13. Kronecker sum solver */
    printf("  [%d] lyap_solve_kleinman... ", ++test_count); fflush(stdout);
    { double A[9]={-1,0,0,0,-2,0,0,0,-3}, Q[9]={1,0,0,0,1,0,0,0,1};
      LDM_QuadraticForm* qf=ldm_lyap_solve_kleinman(A,Q,3,100,1e-8);
      assert(qf); assert(qf->is_pd); ldm_quadratic_free(qf); }
    printf("PASSED\n"); fflush(stdout);

    /* 14. Nonlinear stability */
    printf("  [%d] nonlinear_stability... ", ++test_count); fflush(stdout);
    { LDM_System* s=ldm_system_create_linear(1,(double[]){-1});
      LDM_StabilityReport sr=ldm_nonlinear_stability(s,(double[]){1},1,10.0,0.01);
      assert(sr.result>=LDM_STABLE); free(sr.trajectory); ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 15. Energy based */
    printf("  [%d] energy_based... ", ++test_count); fflush(stdout);
    { LDM_System* s=ldm_system_create(2,NULL,NULL);
      s->x[0]=3; s->x[1]=4; NEAR(ldm_energy_based(s,2),12.5,EPS);
      ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 16. Krasovskii candidate */
    printf("  [%d] krasovskii_candidate... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}; LDM_System* s=ldm_system_create_linear(2,A);
      s->x[0]=1; s->x[1]=1; double V=ldm_krasovskii_candidate(s,2);
      assert(V>=0); ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 17. Linear analysis */
    printf("  [%d] linear_analyze... ", ++test_count); fflush(stdout);
    { double A[4]={-1,2,-3,-4};
      LDM_LinearAnalysis* la=ldm_linear_analyze(A,2);
      assert(la); assert(la->eigenvalues); ldm_linear_analysis_free(la); }
    printf("PASSED\n"); fflush(stdout);

    /* 18. RoA sublevel */
    printf("  [%d] roa_sublevel... ", ++test_count); fflush(stdout);
    { double P[4]={1,0,0,1}; LDM_RegionEstimate re=ldm_roa_sublevel(NULL,P,2,10.0);
      assert(re.c_max==10.0); assert(re.volume_estimate>0); ldm_region_free(&re); }
    printf("PASSED\n"); fflush(stdout);

    /* 19. Converse sampling */
    printf("  [%d] converse_sampling... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}; LDM_System* s=ldm_system_create_linear(2,A);
      LDM_ConverseResult cr=ldm_converse_sampling(s,2,100,2.0);
      assert(cr.V_grid); assert(cr.is_valid); ldm_converse_free(&cr);
      ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 20. Settling time */
    printf("  [%d] settling_time... ", ++test_count); fflush(stdout);
    { double Ts=ldm_linear_settling_time_estimate((double[]){-1,0,0,-2},2,0.02);
      assert(Ts>0); assert(Ts<100); }
    printf("PASSED\n"); fflush(stdout);

    /* 21. Overshoot bound */
    printf("  [%d] overshoot_bound... ", ++test_count); fflush(stdout);
    { NEAR(ldm_linear_overshoot_bound((double[]){2,0,0,2},2),1.0,0.1); }
    printf("PASSED\n"); fflush(stdout);

    /* 22. Massera condition */
    printf("  [%d] massera_condition... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}; LDM_System* s=ldm_system_create_linear(2,A);
      assert(ldm_massera_condition(s,2)); ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 23. Variable gradient */
    printf("  [%d] variable_gradient... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}; LDM_System* s=ldm_system_create_linear(2,A);
      double P[4]; bool ok=ldm_variable_gradient_construct(s,2,P);
      assert(ok); assert(ldm_is_positive_definite(P,2));
      ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 24. Indirect method (linearization) */
    printf("  [%d] indirect_method... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}; LDM_System* s=ldm_system_create_linear(2,A);
      double x_eq[2]={0,0}; LDM_StabilityResult r=ldm_indirect_method(s,x_eq,2);
      assert(r==LDM_ASYM_STABLE); ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    /* 25. Jacobian computation */
    printf("  [%d] jacobian... ", ++test_count); fflush(stdout);
    { double A[4]={-1,0,0,-2}; LDM_System* s=ldm_system_create_linear(2,A);
      double J[4], x0[2]={1,1}; ldm_jacobian(s,x0,J,2);
      NEAR(J[0],-1,EPS); NEAR(J[3],-2,EPS); ldm_system_free(s); }
    printf("PASSED\n"); fflush(stdout);

    printf("\n=== %d tests PASSED ===\n", test_count);
    return 0;
}
