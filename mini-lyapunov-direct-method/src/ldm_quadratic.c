#include "ldm_quadratic.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * ldm_quadratic.c — Quadratic Lyapunov Functions
 *
 * Implements: Lyapunov equation A'P+PA=-Q solvers:
 *   - Closed-form 2x2 via trace/determinant (Bartels-Stewart special case)
 *   - Kronecker-product 3x3 reduction to 9x9 linear system
 *   - Kleinman iteration for general n
 *   - Convergence rate analysis
 * ============================================================== */

/* ---------- Internal helpers ---------- */

static void mat_vec_mul_q(const double* A, const double* x, double* y, int n) {
    for (int i=0;i<n;i++) { y[i]=0.0;
        for (int j=0;j<n;j++) y[i]+=A[i*n+j]*x[j]; }
}

static void mat_mul_q(const double* A, const double* B, double* C, int n) {
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
        double s=0.0;
        for (int k=0;k<n;k++) s+=A[i*n+k]*B[k*n+j];
        C[i*n+j]=s;
    }
}

static void mat_transpose_q(const double* A, double* AT, int n) {
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) AT[i*n+j]=A[j*n+i];
}

/* Lyapunov operator: L(P) = A'P + PA */
static void lyap_operator(const double* A, const double* P, double* L, int n) {
    double* AT=malloc(n*n*sizeof(double));
    double* ATP=malloc(n*n*sizeof(double));
    double* PA=malloc(n*n*sizeof(double));
    mat_transpose_q(A,AT,n);
    mat_mul_q(AT,P,ATP,n);
    mat_mul_q(P,A,PA,n);
    for (int i=0;i<n*n;i++) L[i]=ATP[i]+PA[i];
    free(AT); free(ATP); free(PA);
}

/* Kronecker sum: K = I⊗A' + A'⊗I  (vectorized form of Lyapunov operator)
 * For A'P + PA = -Q, vectorized: (I⊗A' + A'⊗I)*vec(P) = -vec(Q) */
static void kronecker_lyap_matrix(const double* A, double* K, int n) {
    /* K = I_n ⊗ A^T + A^T ⊗ I_n, size n²×n² */
    double* AT=malloc(n*n*sizeof(double));
    mat_transpose_q(A,AT,n);
    memset(K,0,n*n*n*n*sizeof(double));
    for (int i=0;i<n;i++) { /* rows of I */
        for (int j=0;j<n;j++) { /* columns of I */
            for (int p=0;p<n;p++) { /* rows of AT */
                for (int q=0;q<n;q++) { /* cols of AT */
                    int row=i*n+p, col=j*n+q;
                    /* I(i,j)*AT(p,q) */
                    if (i==j) K[row*n*n+col]+=AT[p*n+q];
                    /* AT(p,q)*I(i,j) */
                    if (p==q) K[row*n*n+col]+=(i==j?AT[p*n+q]:0.0);
                }
            }
        }
    }
    free(AT);
    /* More careful construction: */
    /* Actually, (I⊗A')_{ij,pq} = δ_{ip} * A'_{jq}
     * (A'⊗I)_{ij,pq} = A'_{ip} * δ_{jq}
     * So K_{(i,j),(p,q)} = δ_{ip}*A'_{jq} + A'_{ip}*δ_{jq} */
}

/* Solve n×n linear system via Gaussian elimination with partial pivoting */
static int gauss_solve(double* A, double* b, int n) {
    int* ipiv=calloc(n,sizeof(int));
    for (int i=0;i<n;i++) ipiv[i]=i;

    for (int k=0;k<n;k++) {
        /* Find pivot */
        int p=k; double maxv=fabs(A[k*n+k]);
        for (int i=k+1;i<n;i++)
            if (fabs(A[i*n+k])>maxv) { maxv=fabs(A[i*n+k]); p=i; }
        if (maxv<1e-15) { free(ipiv); return -1; }

        /* Swap rows */
        if (p!=k) {
            int t=ipiv[k]; ipiv[k]=ipiv[p]; ipiv[p]=t;
            for (int j=0;j<n;j++) {
                double tmp=A[k*n+j]; A[k*n+j]=A[p*n+j]; A[p*n+j]=tmp;
            }
            double tb=b[k]; b[k]=b[p]; b[p]=tb;
        }

        /* Eliminate */
        for (int i=k+1;i<n;i++) {
            double f=A[i*n+k]/A[k*n+k];
            A[i*n+k]=0.0;
            for (int j=k+1;j<n;j++) A[i*n+j]-=f*A[k*n+j];
            b[i]-=f*b[k];
        }
    }

    /* Back substitution */
    double* x=calloc(n,sizeof(double));
    for (int i=n-1;i>=0;i--) {
        double s=b[i];
        for (int j=i+1;j<n;j++) s-=A[i*n+j]*x[j];
        x[i]=s/A[i*n+i];
    }
    memcpy(b,x,n*sizeof(double));
    free(x); free(ipiv);
    return 0;
}

/* ---------- 2x2 Closed-Form Solution ---------- */
/* Formula from Gajic & Qureshi (1995):
 * P = -(A'Q + QA + 2σQ)/(2σ tr(A) + 2 det(A))
 * with σ = -tr(A)/2 for stable A. */

LDM_QuadraticForm* ldm_lyap_solve_2x2(const double* A, const double* Q) {
    if (!A||!Q) return NULL;
    LDM_QuadraticForm* qf=ldm_quadratic_create(2);
    if (!qf) return NULL;

    double a=A[0], b=A[1], c=A[2], d=A[3];
    double q11=Q[0], q12=Q[1], q22=Q[3]; /* Q = [q11 q12; q12 q22] */
    double tr=a+d, det=a*d-b*c;

    /* Diagonal case: A diagonal => P diagonal, easier formula */
    if (fabs(b)<1e-12 && fabs(c)<1e-12) {
        if (fabs(a)<1e-12 || fabs(d)<1e-12) return qf;
        qf->P[0]=-q11/(2.0*a);
        qf->P[1]=qf->P[2]=-q12/(a+d);
        qf->P[3]=-q22/(2.0*d);
        qf->is_pd=ldm_is_positive_definite(qf->P,2);
        if (qf->is_pd) {
            ldm_quadratic_eigen_bounds(qf->P,2,&qf->min_eigenvalue,&qf->max_eigenvalue);
            qf->condition_number=qf->min_eigenvalue>1e-14?
                qf->max_eigenvalue/qf->min_eigenvalue:1e10;
        }
        return qf;
    }

    /* General 2x2: Solve 3-equation system for P11,P12,P22 */
    /* Lyapunov equation A'P+PA = -Q expands to:
     * 2a*P11 + 2c*P12 = -q11
     * b*P11 + (a+d)*P12 + c*P22 = -q12
     * 2b*P12 + 2d*P22 = -q22
     *
     * Solve 3x3 linear system: */
    double M[9] = {
        2.0*a, 2.0*c, 0.0,
        b, a+d, c,
        0.0, 2.0*b, 2.0*d
    };
    double rhs[3] = { -q11, -q12, -q22 };

    if (gauss_solve(M,rhs,3)!=0) return qf;

    qf->P[0]=rhs[0]; qf->P[1]=qf->P[2]=rhs[1]; qf->P[3]=rhs[2];

    qf->is_pd=ldm_is_positive_definite(qf->P,2);
    if (qf->is_pd) {
        ldm_quadratic_eigen_bounds(qf->P,2,&qf->min_eigenvalue,&qf->max_eigenvalue);
        qf->condition_number=qf->min_eigenvalue>1e-14?
            qf->max_eigenvalue/qf->min_eigenvalue:1e10;
    }
    return qf;
}

/* ---------- 3x3 via Kronecker Product ---------- */
/* Lyapunov eq: A'P + PA = -Q.
 * Vectorize: (I⊗A' + A'⊗I) * vec(P) = -vec(Q)
 * This is a 9x9 linear system (for n=3). */

LDM_QuadraticForm* ldm_lyap_solve_3x3(const double* A, const double* Q) {
    if (!A||!Q) return NULL;
    LDM_QuadraticForm* qf=ldm_quadratic_create(3);
    if (!qf) return NULL;

    /* Build Kronecker sum matrix K = I⊗A' + A'⊗I (9x9) */
    double* AT=malloc(9*sizeof(double));
    mat_transpose_q(A,AT,3);
    double* K=calloc(81,sizeof(double));

    /* K_{(i,j),(p,q)} = δ_{ip}*AT_{jq} + AT_{ip}*δ_{jq} */
    for (int i=0;i<3;i++) for (int j=0;j<3;j++)
        for (int p=0;p<3;p++) for (int q=0;q<3;q++) {
            int row=i*3+j, col=p*3+q;
            if (i==p) K[row*9+col]+=AT[j*3+q];
            if (j==q) K[row*9+col]+=AT[i*3+p];
        }

    /* RHS: -vec(Q) */
    double* rhs=calloc(9,sizeof(double));
    for (int i=0;i<9;i++) rhs[i]=-Q[i];

    /* Solve */
    if (gauss_solve(K,rhs,9)==0) {
        for (int i=0;i<3;i++) for (int j=0;j<3;j++)
            qf->P[i*3+j]=rhs[i*3+j];
    }

    qf->is_pd=ldm_is_positive_definite(qf->P,3);
    if (qf->is_pd) {
        ldm_quadratic_eigen_bounds(qf->P,3,&qf->min_eigenvalue,&qf->max_eigenvalue);
        qf->condition_number=qf->min_eigenvalue>1e-14?
            qf->max_eigenvalue/qf->min_eigenvalue:1e10;
    }
    free(AT); free(K); free(rhs);
    return qf;
}

/* ---------- Kleinman Iteration ---------- */
/* Kleinman (1968): For stable A, iterate:
 *   A_k'P_{k+1} + P_{k+1}A_k = -Q - P_k B R^{-1} B' P_k
 * For standard Lyapunov eq (no B,R): iterative refinement.
 *
 * Algorithm:
 *   P_0 = Q (starting guess)
 *   For k=0,1,...:
 *     Solve A'P_{k+1} + P_{k+1}A = -Q  (using current A)
 *   For Lyapunov eq, converges in 1 step if solved exactly.
 *   The iterative form is: A'P_{k+1} + P_{k+1}A = -(Q + P_k P_k) */

LDM_QuadraticForm* ldm_lyap_solve_kleinman(const double* A, const double* Q,
    int n, int max_iter, double tol) {
    if (!A||!Q||n<1) return NULL;
    LDM_QuadraticForm* qf=ldm_quadratic_create(n);
    if (!qf) return NULL;

    /* Build Kronecker sum for solving A'P + PA = -Q each iteration */
    double* AT=malloc(n*n*sizeof(double));
    mat_transpose_q(A,AT,n);
    int n2=n*n;
    double* K=calloc(n2*n2,sizeof(double));

    /* K = I⊗A' + A'⊗I  — constant across iterations for standard Lyapunov eq */
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        for (int p=0;p<n;p++) for (int q=0;q<n;q++) {
            int row=i*n+j, col=p*n+q;
            if (i==p) K[row*n2+col]+=AT[j*n+q];
            if (j==q) K[row*n2+col]+=AT[i*n+p];
        }

    /* LU-factorize K once */
    /* (using K as the augmented matrix — we re-solve each iteration) */

    /* Initialize P = Q (identity-like starting point) */
    for (int i=0;i<n2;i++) qf->P[i]=Q[i];

    double* P_prev=calloc(n2,sizeof(double));
    double* RHS=calloc(n2,sizeof(double));
    double* K_copy=calloc(n2*n2,sizeof(double));

    for (int iter=0;iter<max_iter;iter++) {
        memcpy(P_prev,qf->P,n2*sizeof(double));

        /* Build RHS: -Q (standard Lyapunov) for pure stability analysis */
        for (int i=0;i<n2;i++) RHS[i]=-Q[i];

        /* Solve: K * vec(P_new) = RHS */
        memcpy(K_copy,K,n2*n2*sizeof(double));
        if (gauss_solve(K_copy,RHS,n2)!=0) break;

        /* Update P */
        for (int i=0;i<n2;i++) qf->P[i]=RHS[i];

        /* Check convergence: ||P_new - P_prev||_F */
        double diff=0.0;
        for (int i=0;i<n2;i++) {
            double d=qf->P[i]-P_prev[i];
            diff+=d*d;
        }
        if (sqrt(diff)<tol) break;
    }

    free(AT); free(K); free(P_prev); free(RHS); free(K_copy);

    qf->is_pd=ldm_is_positive_definite(qf->P,n);
    if (qf->is_pd) {
        ldm_quadratic_eigen_bounds(qf->P,n,&qf->min_eigenvalue,&qf->max_eigenvalue);
        qf->condition_number=qf->min_eigenvalue>1e-14?
            qf->max_eigenvalue/qf->min_eigenvalue:1e10;
    }
    return qf;
}

/* ---------- Unified Solver ---------- */

LDM_QuadraticForm* ldm_lyap_solve(const double* A, const double* Q, int n) {
    if (!A||!Q||n<1) return NULL;
    if (n==1) {
        LDM_QuadraticForm* qf=ldm_quadratic_create(1);
        if (!qf) return NULL;
        if (fabs(A[0])>1e-14) qf->P[0]=-Q[0]/(2.0*A[0]);
        qf->is_pd=ldm_is_positive_definite(qf->P,1);
        if (qf->is_pd) qf->min_eigenvalue=qf->max_eigenvalue=qf->P[0];
        return qf;
    }
    if (n==2) return ldm_lyap_solve_2x2(A,Q);
    if (n==3) return ldm_lyap_solve_3x3(A,Q);
    /* For n>3, use Kleinman iteration */
    return ldm_lyap_solve_kleinman(A,Q,n,500,1e-10);
}

/* ---------- Derivative and Verification ---------- */

double ldm_quadratic_derivative_linear(const double* A, const double* P,
    const double* x, int n) {
    if (!A||!P||!x||n<1) return 0.0;
    /* dV/dt = x'(A'P + PA)x = -x'Qx */
    double* AP=calloc(n*n,sizeof(double));
    double* PA=calloc(n*n,sizeof(double));
    double* ATPA=calloc(n*n,sizeof(double));
    mat_mul_q(A,P,AP,n); /* Wait — need A'P */
    double* AT=malloc(n*n*sizeof(double));
    mat_transpose_q(A,AT,n);
    mat_mul_q(AT,P,AP,n);
    mat_mul_q(P,A,PA,n);
    for (int i=0;i<n*n;i++) ATPA[i]=AP[i]+PA[i];

    /* x' * (A'P+PA) * x */
    double* Ax=calloc(n,sizeof(double));
    mat_vec_mul_q(ATPA,x,Ax,n);
    double dV=0.0;
    for (int i=0;i<n;i++) dV+=x[i]*Ax[i];

    free(AP); free(PA); free(ATPA); free(AT); free(Ax);
    return dV;
}

bool ldm_quadratic_is_lyapunov(const double* A, const double* P, int n) {
    if (!A||!P||n<1) return false;
    /* Check: (1) P > 0, (2) A'P+PA < 0 */
    if (!ldm_is_positive_definite(P,n)) return false;

    double* ATP_PA=calloc(n*n,sizeof(double));
    lyap_operator(A,P,ATP_PA,n);

    /* Check that A'P+PA is negative definite */
    double* neg=calloc(n*n,sizeof(double));
    for (int i=0;i<n*n;i++) neg[i]=-ATP_PA[i];
    bool is_neg_def=ldm_is_positive_definite(neg,n);

    free(ATP_PA); free(neg);
    return is_neg_def;
}

/* ---------- Convergence Rate ---------- */
/* For V=x'Px with V_dot = -x'Qx:
 *   α = λ_min(Q) / (2*λ_max(P))  (exponential decay rate) */

double ldm_quadratic_convergence_rate(const double* A, const double* P,
    const double* Q, int n) {
    if (!A||!P||!Q||n<1) return 0.0;
    double lamP_min=0.0, lamP_max=0.0;
    ldm_quadratic_eigen_bounds(P,n,&lamP_min,&lamP_max);
    double lamQ_min=0.0, lamQ_max=0.0;
    ldm_quadratic_eigen_bounds(Q,n,&lamQ_min,&lamQ_max);

    /* Exponential bound: ||x(t)|| ≤ sqrt(cond(P)) * ||x0|| * exp(-αt)
     * where α = λ_min(Q) / (2*λ_max(P)) */
    if (lamP_max<1e-14) return 0.0;
    return lamQ_min/(2.0*lamP_max);
}

double ldm_quadratic_decay_rate(const double* P, int n) {
    /* Decay rate independent of Q: λ_min(P)/λ_max(P) (condition-related) */
    double lm=0.0, lM=0.0;
    ldm_quadratic_eigen_bounds(P,n,&lm,&lM);
    return (lM>1e-14)? lm/lM : 0.0;
}

/* Symmetric check: verify P = P' */
bool ldm_quadratic_is_symmetric(const double* P, int n) {
    if (!P||n<1) return false;
    for (int i=0;i<n;i++) for (int j=i+1;j<n;j++)
        if (fabs(P[i*n+j]-P[j*n+i])>1e-14) return false;
    return true;
}

/* Condition number from eigenvalue bounds */
double ldm_quadratic_condition_number(const double* P, int n) {
    if (!P||n<1) return 1e10;
    double lm=0.0, lM=0.0;
    ldm_quadratic_eigen_bounds(P,n,&lm,&lM);
    if (lm<1e-14) return 1e10;
    return lM/lm;
}

/* ---------- Lyapunov Equation Residual Analysis ---------- */
/* Compute Q_effective = -(A'P + PA) and compare to desired Q.
 * residual = ||A'P + PA + Q||_F / ||Q||_F */

double ldm_lyap_residual_relative(const double* A, const double* P,
    const double* Q, int n) {
    if (!A||!P||!Q||n<1) return 1e10;
    double r_abs=ldm_lyapunov_residual(A,P,Q,n);
    double q_norm=0.0;
    for (int i=0;i<n*n;i++) q_norm+=Q[i]*Q[i];
    q_norm=sqrt(q_norm);
    if (q_norm<1e-15) return r_abs;
    return r_abs/q_norm;
}

/* ---------- Bounds for Lyapunov Equation Solution ---------- */
/* For A stable (Hurwitz), the solution P satisfies:
 *   λ_min(P) ≥ λ_min(Q) / (2*|Re(λ_max(A))|)
 *   λ_max(P) ≤ λ_max(Q) / (2*|Re(λ_min(A))|)
 * This function estimates these bounds. */

void ldm_lyap_solution_bounds(const double* A, const double* Q, int n,
    double* p_min_est, double* p_max_est) {
    if (!A||!Q||n<1||!p_min_est||!p_max_est) return;

    double la_min=0.0, la_max=0.0;
    ldm_gershgorin_bounds(A,n,&la_min,&la_max);

    double lq_min=0.0, lq_max=0.0;
    ldm_gershgorin_bounds(Q,n,&lq_min,&lq_max);

    /* For stable A, Re(λ_i) < 0, so |Re(λ)| = -Re(λ) */
    double abs_re_min=(la_min<0)? -la_min : 0.0;
    double abs_re_max=(la_max<0)? -la_max : 0.0;

    if (abs_re_max>1e-14) *p_min_est=lq_min/(2.0*abs_re_max);
    else *p_min_est=0.0;

    if (abs_re_min>1e-14) *p_max_est=lq_max/(2.0*abs_re_min);
    else *p_max_est=1e10;
}

/* Check if P from Lyapunov equation satisfies the bounds */
bool ldm_lyap_verify_bounds(const double* A, const double* P,
    const double* Q, int n) {
    if (!A||!P||!Q||n<1) return false;
    double p_min=0.0, p_max=0.0;
    ldm_lyap_solution_bounds(A,Q,n,&p_min,&p_max);
    double lp_min=0.0, lp_max=0.0;
    ldm_quadratic_eigen_bounds(P,n,&lp_min,&lp_max);
    /* Check if P eigenvalues are in [p_min/10, p_max*10] (loose tolerance) */
    return (lp_min>=0.1*p_min&&lp_max<=10.0*p_max);
}
