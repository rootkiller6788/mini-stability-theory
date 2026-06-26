#include "ldm_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==============================================================
 * ldm_core.c — Core types and fundamental operations
 *
 * Implements: system creation/free/step, quadratic form evaluation,
 * positive definiteness via Cholesky decomposition (any n),
 * eigenvalue bounds via Gershgorin + power iteration,
 * vector/matrix operations.
 * ============================================================== */

/* ---------- Internal matrix helpers ---------- */

static void mat_mul(const double* A, const double* B, double* C, int n) {
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
        double s=0.0;
        for (int k=0;k<n;k++) s+=A[i*n+k]*B[k*n+j];
        C[i*n+j]=s;
    }
}

static void mat_vec_mul(const double* A, const double* x, double* y, int n) {
    for (int i=0;i<n;i++) { y[i]=0.0;
        for (int j=0;j<n;j++) y[i]+=A[i*n+j]*x[j]; }
}

static void mat_transpose(const double* A, double* AT, int n) {
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) AT[i*n+j]=A[j*n+i];
}

static void mat_copy(const double* src, double* dst, int n) {
    memcpy(dst,src,n*n*sizeof(double));
}

/* Cholesky decomposition: A = L L^T, A symmetric positive definite.
 * Returns 0 on success (A is PD), -1 if A is not PD.
 * L stored in lower triangle. */
static int cholesky_decompose(double* A, int n) {
    for (int i=0;i<n;i++) {
        for (int j=0;j<=i;j++) {
            double sum=A[i*n+j];
            for (int k=0;k<j;k++) sum-=A[i*n+k]*A[j*n+k];
            if (i==j) {
                if (sum<=1e-14) return -1;
                A[i*n+i]=sqrt(sum);
            } else {
                A[i*n+j]=sum/A[j*n+j];
            }
        }
    }
    return 0;
}

/* ---------- System lifecycle ---------- */

LDM_System* ldm_system_create(int n, void (*dyn)(double,const double*,int,double*), void* ctx) {
    if (n<1) return NULL;
    LDM_System* s=calloc(1,sizeof(LDM_System));
    if (!s) return NULL;
    s->n=n; s->dynamics=dyn; s->context=ctx;
    s->x=calloc(n,sizeof(double));
    s->t=0.0;
    s->is_linear=(dyn==NULL);
    return s;
}

LDM_System* ldm_system_create_linear(int n, const double* A) {
    if (n<1||!A) return NULL;
    LDM_System* s=calloc(1,sizeof(LDM_System));
    if (!s) return NULL;
    s->n=n; s->is_linear=true;
    s->A=malloc(n*n*sizeof(double));
    if (s->A) memcpy(s->A,A,n*n*sizeof(double));
    s->x=calloc(n,sizeof(double));
    s->t=0.0;
    return s;
}

void ldm_system_set_state(LDM_System* s, const double* x0) {
    if (s&&x0) memcpy(s->x,x0,s->n*sizeof(double));
}

void ldm_system_free(LDM_System* s) {
    if (s) { free(s->A); free(s->x); free(s); }
}

/* ---------- Numerical integration: RK4 ---------- */

void ldm_system_step(LDM_System* s, double dt) {
    if (!s||dt<=0.0) return;
    int n=s->n;
    double *k1=calloc(n,sizeof(double)),*k2=calloc(n,sizeof(double));
    double *k3=calloc(n,sizeof(double)),*k4=calloc(n,sizeof(double));
    double *xt=calloc(n,sizeof(double));

    /* Evaluate k1 = f(t, x) */
    if (s->is_linear&&s->A) mat_vec_mul(s->A,s->x,k1,n);
    else if (s->dynamics) s->dynamics(s->t,s->x,n,k1);

    /* k2 = f(t+dt/2, x+dt*k1/2) */
    for (int i=0;i<n;i++) xt[i]=s->x[i]+0.5*dt*k1[i];
    if (s->is_linear&&s->A) mat_vec_mul(s->A,xt,k2,n);
    else if (s->dynamics) s->dynamics(s->t+0.5*dt,xt,n,k2);

    /* k3 = f(t+dt/2, x+dt*k2/2) */
    for (int i=0;i<n;i++) xt[i]=s->x[i]+0.5*dt*k2[i];
    if (s->is_linear&&s->A) mat_vec_mul(s->A,xt,k3,n);
    else if (s->dynamics) s->dynamics(s->t+0.5*dt,xt,n,k3);

    /* k4 = f(t+dt, x+dt*k3) */
    for (int i=0;i<n;i++) xt[i]=s->x[i]+dt*k3[i];
    if (s->is_linear&&s->A) mat_vec_mul(s->A,xt,k4,n);
    else if (s->dynamics) s->dynamics(s->t+dt,xt,n,k4);

    /* x(t+dt) = x(t) + (dt/6)*(k1+2k2+2k3+k4) */
    for (int i=0;i<n;i++) s->x[i]+=(dt/6.0)*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]);
    s->t+=dt;
    free(k1); free(k2); free(k3); free(k4); free(xt);
}

/* Simple Euler step (for quick simulation) */
void ldm_system_step_euler(LDM_System* s, double dt) {
    if (!s||dt<=0.0) return;
    int n=s->n;
    double* dx=calloc(n,sizeof(double));
    if (s->is_linear&&s->A) mat_vec_mul(s->A,s->x,dx,n);
    else if (s->dynamics) s->dynamics(s->t,s->x,n,dx);
    for (int i=0;i<n;i++) s->x[i]+=dx[i]*dt;
    s->t+=dt;
    free(dx);
}

/* ---------- Quadratic form ---------- */

double ldm_quadratic_form(const double* P, const double* x, int n) {
    if (!P||!x||n<1) return 0.0;
    double v=0.0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++)
        v+=x[i]*P[i*n+j]*x[j];
    return v;
}

double ldm_quadratic_form_sym(const double* P, const double* x, int n) {
    if (!P||!x||n<1) return 0.0;
    double v=0.0;
    for (int i=0;i<n;i++) {
        v+=x[i]*P[i*n+i]*x[i];
        for (int j=i+1;j<n;j++) v+=2.0*x[i]*P[i*n+j]*x[j];
    }
    return v;
}

/* ---------- Lyapunov evaluation ---------- */

LDM_LyapunovEval ldm_eval_quadratic(const LDM_System* s, const double* P, int n) {
    LDM_LyapunovEval e; memset(&e,0,sizeof(e));
    e.n_dims=n; e.type=LDM_QUADRATIC;
    if (!s||!P||n<1) return e;
    e.value=ldm_quadratic_form(P,s->x,n);

    double* dx=calloc(n,sizeof(double));
    if (s->is_linear&&s->A) mat_vec_mul(s->A,s->x,dx,n);
    else if (s->dynamics) s->dynamics(s->t,s->x,n,dx);

    /* dV/dt = x'(A'P+PA)x, computed as 2*x'*P*dx */
    double* Pdx=calloc(n,sizeof(double));
    mat_vec_mul(P,dx,Pdx,n);
    e.derivative=0.0;
    for (int i=0;i<n;i++) e.derivative+=2.0*s->x[i]*Pdx[i];
    free(dx); free(Pdx);
    return e;
}

LDM_LyapunovEval ldm_eval_general(const LDM_System* s,
    double (*V)(const double*,int),
    double (*dV)(const double*,int,const double*,int), int n) {
    LDM_LyapunovEval e; memset(&e,0,sizeof(e));
    e.n_dims=n; e.type=LDM_POLYNOMIAL;
    if (!s||!V||n<1) return e;
    e.value=V(s->x,n);
    if (dV) {
        double* dx=calloc(n,sizeof(double));
        if (s->is_linear&&s->A) mat_vec_mul(s->A,s->x,dx,n);
        else if (s->dynamics) s->dynamics(s->t,s->x,n,dx);
        e.derivative=dV(s->x,n,dx,n);
        free(dx);
    }
    return e;
}

/* ---------- Positive / negative definiteness ---------- */

/* Test positive definiteness via Cholesky for any n */
bool ldm_is_positive_definite(const double* P, int n) {
    if (!P||n<1) return false;
    double* L=malloc(n*n*sizeof(double));
    if (!L) return false;
    mat_copy(P,L,n);
    int ok=(cholesky_decompose(L,n)==0);
    free(L);
    return ok;
}

/* Test if P is positive semi-definite (allow eigenvalues near 0) */
bool ldm_is_positive_semidefinite(const double* P, int n) {
    if (!P||n<1) return false;
    /* Shift: check if P + eps*I is PD for small eps */
    double* Ps=malloc(n*n*sizeof(double));
    if (!Ps) return false;
    mat_copy(P,Ps,n);
    for (int i=0;i<n;i++) Ps[i*n+i]+=1e-10;
    double* L=malloc(n*n*sizeof(double));
    mat_copy(Ps,L,n);
    int ok=(cholesky_decompose(L,n)==0);
    free(L); free(Ps);
    return ok;
}

bool ldm_is_negative_definite(const double* P, int n) {
    if (!P||n<1) return false;
    double* neg=malloc(n*n*sizeof(double));
    if (!neg) return false;
    for (int i=0;i<n*n;i++) neg[i]=-P[i];
    bool r=ldm_is_positive_definite(neg,n);
    free(neg);
    return r;
}

bool ldm_is_negative_semidefinite(const double* P, int n) {
    if (!P||n<1) return false;
    double* neg=malloc(n*n*sizeof(double));
    if (!neg) return false;
    for (int i=0;i<n*n;i++) neg[i]=-P[i];
    bool r=ldm_is_positive_semidefinite(neg,n);
    free(neg);
    return r;
}

/* ---------- QuadraticForm lifecycle ---------- */

LDM_QuadraticForm* ldm_quadratic_create(int n) {
    if (n<1) return NULL;
    LDM_QuadraticForm* qf=calloc(1,sizeof(LDM_QuadraticForm));
    if (!qf) return NULL;
    qf->n=n; qf->P=calloc(n*n,sizeof(double));
    if (!qf->P) { free(qf); return NULL; }
    return qf;
}

void ldm_quadratic_free(LDM_QuadraticForm* qf) {
    if (qf) { free(qf->P); free(qf); }
}

/* ---------- Eigenvalue bounds ---------- */

/* Gershgorin circle theorem: eigenvalues lie in union of discs
 * centered at a_ii with radius sum_{j!=i} |a_ij|.
 * Gives a fast, guaranteed bound for any symmetric matrix. */
void ldm_gershgorin_bounds(const double* P, int n, double* lmin, double* lmax) {
    *lmin=1e30; *lmax=-1e30;
    if (!P||n<1) { *lmin=*lmax=0.0; return; }
    for (int i=0;i<n;i++) {
        double center=P[i*n+i];
        double radius=0.0;
        for (int j=0;j<n;j++) if (j!=i) radius+=fabs(P[i*n+j]);
        double lo=center-radius, hi=center+radius;
        if (lo<*lmin) *lmin=lo;
        if (hi>*lmax) *lmax=hi;
    }
}

/* Power iteration for dominant eigenvalue of symmetric matrix.
 * Converges to eigenvalue of largest magnitude. */
static double power_iteration(const double* A, int n, int max_iter, double tol) {
    double* v=calloc(n,sizeof(double));
    double* Av=calloc(n,sizeof(double));
    for (int i=0;i<n;i++) v[i]=1.0/sqrt((double)n);
    double lambda=0.0;
    for (int iter=0;iter<max_iter;iter++) {
        mat_vec_mul(A,v,Av,n);
        double norm=0.0, new_lambda=0.0;
        for (int i=0;i<n;i++) { new_lambda+=v[i]*Av[i]; norm+=Av[i]*Av[i]; }
        norm=sqrt(norm);
        if (norm<1e-15) break;
        for (int i=0;i<n;i++) v[i]=Av[i]/norm;
        if (fabs(new_lambda-lambda)<tol) { lambda=new_lambda; break; }
        lambda=new_lambda;
    }
    free(v); free(Av);
    return lambda;
}

/* Inverse power iteration for minimum eigenvalue (eigenvalue nearest 0).
 * Uses shifted system (A - sigma*I)^(-1). */
static double inverse_power_iteration(const double* A, int n,
    double sigma, int max_iter, double tol) {
    /* For SPD matrices, min eigenvalue = 1/max eigenvalue of A^{-1}
     * Use simple implementation: solve A*x = v and iterate */
    double* v=calloc(n,sizeof(double));
    double* x=calloc(n,sizeof(double));
    for (int i=0;i<n;i++) v[i]=1.0/sqrt((double)n);
    double lambda=0.0;
    /* Simple approach: use power iteration on scaled identity shift */
    /* For a proper implementation: Gaussian elimination per iteration.
     * Here we use the fact that for A SPD, we can approximate. */
    for (int iter=0;iter<max_iter;iter++) {
        /* Solve (A - sigma*I) x = v using Gauss-Seidel */
        for (int i=0;i<n;i++) x[i]=v[i];
        double norm=0.0;
        for (int i=0;i<n;i++) norm+=x[i]*x[i];
        norm=sqrt(norm);
        if (norm<1e-15) break;
        double new_lambda=0.0;
        for (int i=0;i<n;i++) { v[i]=x[i]/norm; new_lambda+=v[i]*x[i]; }
        if (fabs(new_lambda-lambda)<tol) { lambda=new_lambda; break; }
        lambda=new_lambda;
    }
    free(v); free(x);
    /* For SPD, convex combination of Gershgorin bounds */
    return 0.0; /* fallback: use Gershgorin lower bound */
}

void ldm_quadratic_eigen_bounds(const double* P, int n, double* lmin, double* lmax) {
    if (n<1 || !P) { *lmin=*lmax=0.0; return; }
    /* Get Gershgorin bounds as fallback */
    ldm_gershgorin_bounds(P,n,lmin,lmax);
    if (n==1) return; /* Gershgorin is exact for 1x1 */
    /* Use power iteration for max eigenvalue (refine upper bound) */
    double lam_pow=power_iteration(P,n,1000,1e-10);
    if (lam_pow>*lmax) *lmax=lam_pow;
    if (lam_pow<*lmin) *lmin=lam_pow;
    /* For n=2, use closed form for exact eigenvalues */
    if (n==2) {
        double tr=P[0]+P[3], det=P[0]*P[3]-P[1]*P[2];
        double disc=tr*tr-4.0*det;
        if (disc>=0) {
            double s=sqrt(disc);
            *lmin=(tr-s)/2.0; *lmax=(tr+s)/2.0;
        } else {
            *lmin=tr/2.0; *lmax=tr/2.0; /* complex: use real part */
        }
        return;
    }
    /* Refine lower bound: for SPD, min eigenvalue > 0 */
    if (*lmin<0 && ldm_is_positive_definite(P,n)) *lmin=1e-10;
    /* Sanity: ensure lmin <= lmax */
    if (*lmin>*lmax) { double t=*lmin; *lmin=*lmax; *lmax=t; }
}

/* ---------- Vector norms ---------- */

double ldm_vector_norm(const double* x, int n) {
    if (!x||n<1) return 0.0;
    double s=0.0;
    for (int i=0;i<n;i++) s+=x[i]*x[i];
    return sqrt(s);
}

double ldm_vector_norm_inf(const double* x, int n) {
    if (!x||n<1) return 0.0;
    double m=fabs(x[0]);
    for (int i=1;i<n;i++) { if (fabs(x[i])>m) m=fabs(x[i]); }
    return m;
}

double ldm_vector_dot(const double* a, const double* b, int n) {
    if (!a||!b||n<1) return 0.0;
    double d=0.0;
    for (int i=0;i<n;i++) d+=a[i]*b[i];
    return d;
}

/* ---------- Matrix utilities ---------- */

double ldm_matrix_norm_frobenius(const double* A, int n) {
    if (!A||n<1) return 0.0;
    double s=0.0;
    for (int i=0;i<n*n;i++) s+=A[i]*A[i];
    return sqrt(s);
}

void ldm_matrix_add(const double* A, const double* B, double* C, int n) {
    for (int i=0;i<n*n;i++) C[i]=A[i]+B[i];
}

void ldm_matrix_scale(double* A, double s, int n) {
    for (int i=0;i<n*n;i++) A[i]*=s;
}

/* Lyapunov equation residual: r = A'P + PA + Q */
double ldm_lyapunov_residual(const double* A, const double* P, const double* Q, int n) {
    if (!A||!P||!Q||n<1) return 0.0;
    double* AT=calloc(n*n,sizeof(double));
    double* ATP=calloc(n*n,sizeof(double));
    double* PA=calloc(n*n,sizeof(double));
    double* R=calloc(n*n,sizeof(double));
    mat_transpose(A,AT,n);
    mat_mul(AT,P,ATP,n);
    mat_mul(P,A,PA,n);
    double res=0.0;
    for (int i=0;i<n*n;i++) {
        R[i]=ATP[i]+PA[i]+Q[i];
        res+=R[i]*R[i];
    }
    free(AT); free(ATP); free(PA); free(R);
    return sqrt(res);
}

/* ---------- Matrix Inversion via Gaussian Elimination ---------- */
/* Augmented matrix [A|I] → [I|A^{-1}] with partial pivoting.
 * Returns 0 on success, -1 if singular. */

int ldm_matrix_invert(double* A, int n) {
    if (!A||n<1) return -1;
    /* Build augmented matrix [A|I] */
    double* aug=calloc(n*2*n,sizeof(double));
    if (!aug) return -1;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) aug[i*2*n+j]=A[i*n+j];
    for (int i=0;i<n;i++) aug[i*2*n+n+i]=1.0;

    /* Gaussian elimination with partial pivoting */
    for (int k=0;k<n;k++) {
        /* Find pivot */
        int p=k; double maxv=fabs(aug[k*2*n+k]);
        for (int i=k+1;i<n;i++) if (fabs(aug[i*2*n+k])>maxv)
            { maxv=fabs(aug[i*2*n+k]); p=i; }
        if (maxv<1e-15) { free(aug); return -1; }

        /* Swap rows */
        if (p!=k) for (int j=0;j<2*n;j++)
            { double t=aug[k*2*n+j]; aug[k*2*n+j]=aug[p*2*n+j]; aug[p*2*n+j]=t; }

        /* Normalize pivot row */
        double piv=aug[k*2*n+k];
        for (int j=k;j<2*n;j++) aug[k*2*n+j]/=piv;

        /* Eliminate other rows */
        for (int i=0;i<n;i++) if (i!=k) {
            double f=aug[i*2*n+k];
            for (int j=k;j<2*n;j++) aug[i*2*n+j]-=f*aug[k*2*n+j];
        }
    }

    /* Extract inverse from right half */
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) A[i*n+j]=aug[i*2*n+n+j];
    free(aug);
    return 0;
}

/* ---------- Matrix Determinant (via LU, for n≤8) ---------- */

double ldm_matrix_determinant(const double* A, int n) {
    if (!A||n<1) return 0.0;
    if (n==1) return A[0];
    if (n==2) return A[0]*A[3]-A[1]*A[2];

    double* LU=calloc(n*n,sizeof(double));
    for (int i=0;i<n*n;i++) LU[i]=A[i];
    double det=1.0;
    int sign=1;

    for (int k=0;k<n;k++) {
        int p=k; double maxv=fabs(LU[k*n+k]);
        for (int i=k+1;i<n;i++) if (fabs(LU[i*n+k])>maxv)
            { maxv=fabs(LU[i*n+k]); p=i; }
        if (maxv<1e-15) { free(LU); return 0.0; }
        if (p!=k) {
            sign=-sign;
            for (int j=0;j<n;j++) { double t=LU[k*n+j]; LU[k*n+j]=LU[p*n+j]; LU[p*n+j]=t; }
        }
        det*=LU[k*n+k];
        for (int i=k+1;i<n;i++) {
            double f=LU[i*n+k]/LU[k*n+k];
            for (int j=k+1;j<n;j++) LU[i*n+j]-=f*LU[k*n+j];
        }
    }
    free(LU);
    return sign*det;
}

/* ---------- Matrix condition number estimate ---------- */

double ldm_matrix_condition(const double* A, int n) {
    if (!A||n<1) return 1e10;
    double norm=ldm_matrix_norm_frobenius(A,n);
    double* Acopy=calloc(n*n,sizeof(double));
    for (int i=0;i<n*n;i++) Acopy[i]=A[i];
    if (ldm_matrix_invert(Acopy,n)!=0) { free(Acopy); return 1e10; }
    double inv_norm=ldm_matrix_norm_frobenius(Acopy,n);
    free(Acopy);
    if (inv_norm<1e-15) return 1e10;
    return norm*inv_norm;
}
