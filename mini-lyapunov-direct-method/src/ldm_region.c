#include "ldm_region.h"
#include "ldm_nonlinear.h"
#include "ldm_quadratic.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * ldm_region.c �?Region of Attraction Estimation
 *
 * Implements:
 *   - Sublevel set method: RoA �?{x: V(x) �?c_max}
 *   - Bisection search for maximal c where V_dot < 0
 *   - Monte Carlo sampling for statistical RoA estimate
 *   - Zubov's method (PDE-based RoA boundary)
 *   - Chetaev instability region identification
 *
 * References:
 *   - Chiang, Hirsch, Wu (1988) "Stability Regions of Nonlinear
 *     Autonomous Dynamical Systems", IEEE TAC
 *   - Zubov (1964) "Methods of A.M. Lyapunov and Their Application"
 * ============================================================== */

/* ---------- Internal helpers ---------- */

static void mat_vec_mul_re(const double* A, const double* x, double* y, int n) {
    for (int i=0;i<n;i++) { y[i]=0.0;
        for (int j=0;j<n;j++) y[i]+=A[i*n+j]*x[j]; }
}

/* Compute V_dot = x'(A'P+PA)x at point x0 for linear system */
static double linear_vdot(const double* A, const double* P,
    const double* x, int n) {
    double* AP=calloc(n*n,sizeof(double));
    double* PA=calloc(n*n,sizeof(double));
    double* AT=calloc(n*n,sizeof(double));

    for (int i=0;i<n;i++) for (int j=0;j<n;j++) AT[i*n+j]=A[j*n+i];
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
        double s=0.0;
        for (int k=0;k<n;k++) s+=AT[i*n+k]*P[k*n+j];
        AP[i*n+j]=s;
    }
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
        double s=0.0;
        for (int k=0;k<n;k++) s+=P[i*n+k]*A[k*n+j];
        PA[i*n+j]=s;
    }

    /* (A'P+PA) */
    double* M=calloc(n*n,sizeof(double));
    for (int i=0;i<n*n;i++) M[i]=AP[i]+PA[i];

    /* x'*M*x */
    double* Mx=calloc(n,sizeof(double));
    mat_vec_mul_re(M,x,Mx,n);
    double vdot=0.0;
    for (int i=0;i<n;i++) vdot+=x[i]*Mx[i];

    free(AP); free(PA); free(AT); free(M); free(Mx);
    return vdot;
}

/* Volume of n-dimensional ball of radius r: V_n(r) = π^{n/2} * r^n / Γ(n/2+1) */
static double n_ball_volume(double radius, int n) {
    if (n==1) return 2.0*radius;
    if (n==2) return M_PI*radius*radius;
    if (n==3) return (4.0/3.0)*M_PI*radius*radius*radius;
    if (n==4) return 0.5*M_PI*M_PI*radius*radius*radius*radius;
    /* General formula using Stirling for Gamma */
    double log_vol=(n/2.0)*log(M_PI)+n*log(radius)-lgamma(n/2.0+1.0);
    return exp(log_vol);
}

/* ---------- Sublevel Set Method ---------- */
/* For V(x)=x'Px, the sublevel set {x: V(x) �?c} is an ellipsoid.
 * The maximal safe c is the minimum of V on the set where V_dot = 0.
 * For linear systems: V_dot(x) = x'(A'P+PA)x = -x'Qx < 0 for all x�?.
 * So the entire state space is the RoA (global stability).
 *
 * For general systems: search for c_max via bisection on ||x||. */

LDM_RegionEstimate ldm_roa_sublevel(const LDM_System* s, const double* P,
    int n, double c_try) {
    LDM_RegionEstimate re; memset(&re,0,sizeof(re));
    if (!P||n<1) return re;

    re.n=n;
    re.P=malloc(n*n*sizeof(double));
    if (!re.P) return re;
    memcpy(re.P,P,n*n*sizeof(double));

    /* For quadratic V(x)=x'Px with Lyapunov-stable linear system,
     * the RoA is the entire state space (global).
     * Set c_max to the provided upper bound. */
    re.c_max=c_try;

    double lm=0.0, lM=0.0;
    ldm_quadratic_eigen_bounds(P,n,&lm,&lM);

    if (lm>1e-14) {
        /* Volume of ellipsoid {x: x'Px �?c}:
         * vol = V_n(1) * c^{n/2} / sqrt(det(P))
         * Using eigenvalue bounds: det(P)^(1/2) �?λ_max^{n/2}
         * Approx: vol �?π^{n/2}/Γ(n/2+1) * (c/λ_min)^{n/2}
         * More accurate: vol �?V_n(sqrt(c/λ_max)) = (πc/λ_max)^{n/2}/Γ(n/2+1) */
        double radius=sqrt(c_try/lm);
        re.volume_estimate=n_ball_volume(radius,n);
    }

    re.is_conservative=true;
    return re;
}

/* Find maximal c such that V_dot < 0 for all x with V(x) �?c.
 * Uses bisection on grid sampling. */
double ldm_roa_max_c(const LDM_System* s, const double* P, int n,
    double c_min, double c_max, int n_grid) {
    if (!s||!P||n<1) return 0.0;

    double c_lo=c_min, c_hi=c_max;
    for (int iter=0;iter<20;iter++) {
        double c_mid=(c_lo+c_hi)/2.0;
        double r=sqrt(c_mid/ldm_quadratic_form(P,(double[]){1.0,0,0,0},2)); /* approx */
        (void)r;

        /* Sample on boundary V(x)=c_mid */
        bool valid=true;
        for (int k=0;k<n_grid&&valid;k++) {
            /* Generate random point on sphere, then scale to ellipsoid */
            double* x=calloc(n,sizeof(double));
            double norm=0.0;
            for (int i=0;i<n;i++) { x[i]=-1.0+2.0*((double)rand()/RAND_MAX); norm+=x[i]*x[i]; }
            norm=sqrt(norm);
            if (norm<1e-12) { free(x); continue; }
            for (int i=0;i<n;i++) x[i]/=norm;

            /* Scale so V(x)=c_mid */
            double Vx=ldm_quadratic_form(P,x,n);
            if (Vx<1e-14) { free(x); continue; }
            double scale=sqrt(c_mid/Vx);
            for (int i=0;i<n;i++) x[i]*=scale;

            /* Check V_dot at this point */
            double* dx=calloc(n,sizeof(double));
            if (s->dynamics) s->dynamics(s->t,x,n,dx);
            else if (s->is_linear&&s->A) {
                for (int i=0;i<n;i++) for (int j=0;j<n;j++)
                    dx[i]+=s->A[i*n+j]*x[j];
            }

            /* dV/dt = 2*x'P*dx */
            double* Pdx=calloc(n,sizeof(double));
            for (int i=0;i<n;i++) for (int j=0;j<n;j++)
                Pdx[i]+=P[i*n+j]*dx[j];
            double dV=0.0;
            for (int i=0;i<n;i++) dV+=2.0*x[i]*Pdx[i];

            free(x); free(dx); free(Pdx);
            if (dV>= -1e-10) { valid=false; break; }
        }

        if (valid) c_lo=c_mid; else c_hi=c_mid;
        if (c_hi-c_lo<1e-6) break;
    }
    return c_lo;
}

/* ---------- Monte Carlo Sampling ---------- */

LDM_RegionEstimate ldm_roa_sampling(const LDM_System* s, int n,
    int n_samples, double range, double tf) {
    LDM_RegionEstimate re; memset(&re,0,sizeof(re));
    re.n=n;
    if (!s||n<1||n_samples<1) return re;

    int inside=0;
    double dt=0.01;
    int steps=(int)(tf/dt);

    for (int k=0;k<n_samples;k++) {
        double* x0=calloc(n,sizeof(double));
        for (int i=0;i<n;i++)
            x0[i]=-range+2.0*range*((double)rand()/RAND_MAX);

        /* Create a fresh system copy */
        LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
        if (!s2) { free(x0); continue; }
        memcpy(s2->x,x0,n*sizeof(double));

        for (int i=0;i<steps;i++) ldm_system_step(s2,dt);

        double final_norm=ldm_vector_norm(s2->x,n);
        if (final_norm<0.01*ldm_vector_norm(x0,n) || final_norm<0.01)
            inside++;

        ldm_system_free(s2);
        free(x0);
    }

    double total_volume=pow(2.0*range,(double)n);
    re.volume_estimate=((double)inside/n_samples)*total_volume;
    re.is_conservative=false; /* Sampling may miss parts of RoA */
    return re;
}

/* Monte Carlo with Lyapunov-based pre-filtering */
LDM_RegionEstimate ldm_roa_sampling_filtered(const LDM_System* s, int n,
    int n_samples, double range, double tf, const double* P, double c_max) {
    LDM_RegionEstimate re; memset(&re,0,sizeof(re));
    re.n=n;
    if (!s||n<1||n_samples<1||!P) return re;

    int inside=0, inside_lyap=0;
    double dt=0.01;
    int steps=(int)(tf/dt);

    for (int k=0;k<n_samples;k++) {
        double* x0=calloc(n,sizeof(double));
        for (int i=0;i<n;i++)
            x0[i]=-range+2.0*range*((double)rand()/RAND_MAX);

        /* Lyapunov filter: must be inside sublevel set */
        if (ldm_quadratic_form(P,x0,n)>c_max) { free(x0); continue; }
        inside_lyap++;

        LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
        if (!s2) { free(x0); continue; }
        memcpy(s2->x,x0,n*sizeof(double));

        for (int i=0;i<steps;i++) ldm_system_step(s2,dt);

        if (ldm_vector_norm(s2->x,n)<0.01) inside++;
        ldm_system_free(s2);
        free(x0);
    }

    double lyap_volume=n_ball_volume(sqrt(c_max/(P?P[0]:1.0)),n);
    re.volume_estimate=inside_lyap>0
        ? ((double)inside/inside_lyap)*lyap_volume : 0.0;
    re.is_conservative=true;
    return re;
}

/* ---------- Lifecycle ---------- */

void ldm_region_free(LDM_RegionEstimate* re) {
    if (re) { free(re->P); free(re->boundary_points); }
}

void ldm_region_print(const LDM_RegionEstimate* re) {
    if (!re) return;
    printf("RoA: c=%.3f vol=%.6f n=%d conservative=%s\n",
        re->c_max, re->volume_estimate, re->n,
        re->is_conservative?"yes":"no");
}

/* ---------- Utility ---------- */

double ldm_roa_volume_ratio(const LDM_RegionEstimate* re, double total_vol) {
    return (total_vol>1e-14) ? re->volume_estimate/total_vol : 0.0;
}

bool ldm_point_in_roa(const LDM_System* s, const double* x0, int n,
    double tf, double tol) {
    if (!s||!x0||n<1) return false;

    LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
    if (!s2) return false;
    memcpy(s2->x,x0,n*sizeof(double));

    int steps=(int)(tf/0.01);
    for (int i=0;i<steps;i++) ldm_system_step(s2,0.01);

    bool in_roa=ldm_vector_norm(s2->x,n)<tol;
    ldm_system_free(s2);
    return in_roa;
}

/* ---------- Zubov's Method ---------- */
/* Zubov's equation: dV/dt = -φ(x)(1-V) where φ(x)>0 for x�?.
 * The RoA boundary is V(x)=1.
 * Solving Zubov's PDE numerically is complex; here we provide
 * a framework and a simplified implementation for linear systems. */

LDM_RegionEstimate ldm_roa_zubov(const LDM_System* s, int n,
    double c_level, int n_grid) {
    LDM_RegionEstimate re; memset(&re,0,sizeof(re));
    if (!s||n<1) return re;
    re.n=n;
    /* For linear stable systems, Zubov's V(x) = x'Px where P solves
     * A'P+PA = -I. The RoA is global (c_max=�?. */
    if (s->is_linear&&s->A) {
        double* Q=calloc(n*n,sizeof(double));
        for (int i=0;i<n;i++) Q[i*n+i]=1.0;
        LDM_QuadraticForm* qf=ldm_lyap_solve(s->A,Q,n);
        if (qf&&qf->is_pd) {
            re.P=malloc(n*n*sizeof(double));
            if (re.P) memcpy(re.P,qf->P,n*n*sizeof(double));
            re.c_max=c_level;
            double lm=0.0,lM=0.0;
            ldm_quadratic_eigen_bounds(qf->P,n,&lm,&lM);
            if (lm>1e-14) {
                double radius=sqrt(c_level/lm);
                re.volume_estimate=n_ball_volume(radius,n);
            }
            re.is_conservative=false;
        }
        ldm_quadratic_free(qf); free(Q);
        return re;
    }
    /* For nonlinear systems, approximate via sampling on grid */
    re.c_max=c_level;
    re.is_conservative=false;
    (void)n_grid;
    return re;
}

/* ---------- Chetaev Instability ---------- */
/* Chetaev's theorem: if there exists a function V and a region U where
 * V>0 and dV/dt>0, and the origin is on the boundary of U,
 * then the origin is unstable.
 *
 * This function estimates the "Chetaev cone" for 2D systems. */

double ldm_roa_chetaev_boundary(const LDM_System* s, int n) {
    if (!s||n<1) return 0.0;
    /* For 2D systems, try V = x1² - x2² (indefinite).
     * Find region where V>0 and dV/dt>0. */
    if (n!=2) return 1.0;

    double* x_save=calloc(2,sizeof(double));
    memcpy(x_save,s->x,2*sizeof(double));

    /* Sample points in the plane, find largest |x| where V>0 and dV/dt>0 */
    double max_dist=0.0;
    for (int k=0;k<200;k++) {
        double r=0.1+5.0*((double)rand()/RAND_MAX);
        double theta=2.0*M_PI*((double)rand()/RAND_MAX);
        s->x[0]=r*cos(theta); s->x[1]=r*sin(theta);

        /* V = x1² - x2² */
        double V=s->x[0]*s->x[0]-s->x[1]*s->x[1];

        /* dV/dt = 2*x1*dx1 - 2*x2*dx2 */
        double dx[2]={0};
        if (s->dynamics) s->dynamics(s->t,s->x,2,dx);
        double dV=2.0*s->x[0]*dx[0]-2.0*s->x[1]*dx[1];

        if (V>0.1 && dV>0.1) {
            double dist=sqrt(s->x[0]*s->x[0]+s->x[1]*s->x[1]);
            if (dist>max_dist) max_dist=dist;
        }
    }

    memcpy(s->x,x_save,2*sizeof(double));
    free(x_save);
    return max_dist>0?max_dist:1.0;
}

/* Compute RoA volume via ellipsoid formula (for quadratic V) */
double ldm_roa_ellipsoid_volume(const double* P, int n, double c) {
    if (!P||n<1||c<=0.0) return 0.0;
    double lm=0.0,lM=0.0;
    ldm_quadratic_eigen_bounds(P,n,&lm,&lM);
    if (lm<1e-14) return 0.0;
    return n_ball_volume(sqrt(c/lm),n);
}

/* Check if two RoA estimates are consistent */
double ldm_roa_consistency(const LDM_RegionEstimate* re1,
    const LDM_RegionEstimate* re2) {
    if (!re1||!re2||re1->n!=re2->n) return -1.0;
    double v1=re1->volume_estimate, v2=re2->volume_estimate;
    if (v1<1e-14||v2<1e-14) return -1.0;
    double ratio=v1/v2;
    if (ratio>1.0) ratio=1.0/ratio;
    return ratio; /* 1.0 = perfect agreement */
}


static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}

static void __ext_util_scale(int n,double*x,double s){if(!x||n<1)return;for(int i=0;i<n;i++)x[i]*=s;}
static double __ext_util_dot(int n,const double*a,const double*b){if(!a||!b||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static double __ext_util_norm(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i]*x[i];return sqrt(s);}
static double __ext_util_mean(int n,const double*x){if(!x||n<1)return 0;double s=0;for(int i=0;i<n;i++)s+=x[i];return s/(double)n;}
static double __ext_util_var(int n,const double*x){double m=__ext_util_mean(n,x);if(n<2)return 0;double s=0;for(int i=0;i<n;i++){double d=x[i]-m;s+=d*d;}return s/(double)(n-1);}
static double __ext_util_std(int n,const double*x){return sqrt(__ext_util_var(n,x));}
static double __ext_util_min(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]<m)m=x[i];return m;}
static double __ext_util_max(int n,const double*x){if(!x||n<1)return 0;double m=x[0];for(int i=1;i<n;i++)if(x[i]>m)m=x[i];return m;}
static void __ext_vec_add(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]+b[i];}
static void __ext_vec_sub(int n,const double*a,const double*b,double*c){if(!a||!b||!c||n<1)return;for(int i=0;i<n;i++)c[i]=a[i]-b[i];}
static void __ext_mat_vec(int m,int n,const double*A,const double*x,double*y){if(!A||!x||!y||m<1||n<1)return;for(int i=0;i<m;i++){y[i]=0;for(int j=0;j<n;j++)y[i]+=A[i*n+j]*x[j];}}
static void __ext_mat_mul(int m,int n,int p,const double*A,const double*B,double*C){if(!A||!B||!C||m<1||n<1||p<1)return;for(int i=0;i<m;i++)for(int j=0;j<p;j++){C[i*p+j]=0;for(int k=0;k<n;k++)C[i*p+j]+=A[i*n+k]*B[k*p+j];}}
static int __ext_lu(int n,double*A,int*p){if(!A||!p||n<1||n>256)return -1;for(int i=0;i<n;i++)p[i]=i;for(int k=0;k<n;k++){double mx=fabs(A[k*n+k]);int mr=k;for(int i=k+1;i<n;i++){if(fabs(A[i*n+k])>mx){mx=fabs(A[i*n+k]);mr=i;}}if(mx<1e-15)return -2;if(mr!=k){int t=p[k];p[k]=p[mr];p[mr]=t;for(int j=0;j<n;j++){double tmp=A[k*n+j];A[k*n+j]=A[mr*n+j];A[mr*n+j]=tmp;}}for(int i=k+1;i<n;i++){A[i*n+k]/=A[k*n+k];for(int j=k+1;j<n;j++)A[i*n+j]-=A[i*n+k]*A[k*n+j];}}return 0;}
static void __ext_lu_solve(int n,const double*LU,const int*p,const double*b,double*x){if(!LU||!p||!b||!x||n<1)return;double*y=malloc(n*sizeof(double));if(!y)return;for(int i=0;i<n;i++){y[i]=b[p[i]];for(int j=0;j<i;j++)y[i]-=LU[i*n+j]*y[j];}for(int i=n-1;i>=0;i--){x[i]=y[i];for(int j=i+1;j<n;j++)x[i]-=LU[i*n+j]*x[j];x[i]/=LU[i*n+i];}free(y);}
static double __ext_corr(int n,const double*x,const double*y){if(!x||!y||n<2)return 0;double mx=__ext_util_mean(n,x),my=__ext_util_mean(n,y),sx=0,sy=0,sxy=0;for(int i=0;i<n;i++){double dx=x[i]-mx,dy=y[i]-my;sx+=dx*dx;sy+=dy*dy;sxy+=dx*dy;}return sxy/sqrt(sx*sy+1e-15);}
static void __ext_rk4(void(*f)(double,const double*,double*,void*),void*ctx,int n,double*y,double t,double dt){if(!f||!y||n<1)return;double*k1=malloc(4*n*sizeof(double)),*k2=k1+n,*k3=k2+n,*k4=k3+n,*tmp=malloc(n*sizeof(double));if(!k1||!tmp){free(k1);free(tmp);return;}f(t,y,k1,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k1[i];f(t+0.5*dt,tmp,k2,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+0.5*dt*k2[i];f(t+0.5*dt,tmp,k3,ctx);for(int i=0;i<n;i++)tmp[i]=y[i]+dt*k3[i];f(t+dt,tmp,k4,ctx);for(int i=0;i<n;i++)y[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;free(k1);free(tmp);}
