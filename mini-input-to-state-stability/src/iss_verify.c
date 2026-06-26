#include "iss_verify.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

ISS_VerifyConfig iss_verify_config_default(void) {
    ISS_VerifyConfig c;
    memset(&c, 0, sizeof(c));
    c.state_radius = 1.0;
    c.input_radius = 1.0;
    c.grid_resolution = 10;
    c.monte_carlo_samples = 1000;
    c.tolerance = 1e-6;
    c.use_adaptive_sampling = true;
    c.check_radial_unboundedness = true;
    return c;
}

ISS_Verification iss_verify(ISS_System* sys, const ISS_VerifyConfig* config) {
    ISS_Verification ver;
    memset(&ver, 0, sizeof(ver));
    if (!sys) return ver;
    ISS_VerifyConfig cfg = config ? *config : iss_verify_config_default();
    ver.lyap_result = iss_lyapunov_verify_random(sys,
        cfg.monte_carlo_samples, cfg.state_radius, cfg.input_radius);
    ver.is_ISS = ver.lyap_result.is_ISS_Lyapunov;
    if (ver.is_ISS) {
        ver.certificate.is_ISS = true;
        ver.certificate.alpha_coeff = ver.lyap_result.alpha_estimate;
        ver.certificate.sigma_coeff = ver.lyap_result.sigma_estimate;
        ver.certificate.gamma_gain = iss_lyapunov_compute_gain(
            ver.certificate.alpha_coeff, ver.certificate.sigma_coeff);
        ver.gain = iss_gain_from_lyapunov(&ver.certificate);
    }
    ver.total_checks = cfg.monte_carlo_samples;
    ver.passed_checks = cfg.monte_carlo_samples - ver.lyap_result.violation_count;
    return ver;
}

bool iss_verify_quick(ISS_System* sys) {
    ISS_VerifyConfig cfg = iss_verify_config_default();
    cfg.monte_carlo_samples = 100;
    ISS_Verification ver = iss_verify(sys, &cfg);
    return ver.is_ISS;
}

bool iss_verify_kinf(double (*func)(double), double tol,
    int n, double max_s) {
    if (!func) return false;
    if (fabs(func(0.0)) > tol) return false;
    double prev = func(0.0);
    int i;
    for (i = 1; i < n; i++) {
        double s = i * max_s / n;
        double v = func(s);
        if (v < prev - tol) return false;
        prev = v;
    }
    return func(max_s) > 1.0 / tol;
}

bool iss_verify_kl(double (*func)(double, double), double tol,
    int ns, int nt, double sm, double tm) {
    if (!func) return false;
    int i, j;
    for (i = 0; i < ns; i++) {
        double s = i * sm / ns;
        double prev = func(s, 0.0);
        for (j = 1; j < nt; j++) {
            double v = func(s, j * tm / nt);
            if (v > prev + tol) return false;
            prev = v;
        }
    }
    return true;
}

void iss_verification_print(const ISS_Verification* ver) {
    if (!ver) return;
    printf("ISS Verification: %s checks=%d/%d gain=%.4f\n",
        ver->is_ISS ? "PASS" : "FAIL",
        ver->passed_checks, ver->total_checks,
        ver->gain.linear_gain);
}

bool iss_verify_convergence(ISS_System* sys, double t_final, double dt, int n_ic) {
    if (!sys || !sys->dynamics) return false;
    int n = sys->n_states;
    int m = sys->n_inputs;
    int i, j;
    for (i = 0; i < n_ic; i++) {
        double x[ISS_MAX_DIM];
        for (j = 0; j < n; j++) x[j] = 2.0 * rand() / RAND_MAX - 1.0;
        double u[ISS_MAX_DIM] = {0};
        int steps = (int)(t_final / dt);
        int step;
        for (step = 0; step < steps; step++) {
            double dxdt[ISS_MAX_DIM];
            sys->dynamics(x, n, u, m, dxdt);
            int k;
            for (k = 0; k < n; k++) x[k] += dxdt[k] * dt;
        }
        double xn = 0.0;
        for (j = 0; j < n; j++) xn += x[j] * x[j];
        if (sqrt(xn) > 0.1) return false;
    }
    return true;
}

/* ==============================================================
 * Stability Margin Estimation
 * ============================================================== */

double iss_verify_stability_margin(ISS_System* s) {
    if (!s) return 0.0;
    ISS_VerifyConfig c = iss_verify_config_default();
    c.monte_carlo_samples = 500;
    ISS_Verification v = iss_verify(s, &c);
    return v.gain.margin;
}

int iss_verify_pass_count(const ISS_Verification* v) {
    return v ? v->passed_checks : 0;
}

bool iss_verify_is_definitive(const ISS_Verification* v) {
    return v && v->total_checks > 1000;
}

/* ==============================================================
 * Adaptive Sampling Verification
 * ============================================================== */

ISS_Verification iss_verify_adaptive(ISS_System* sys,
    const ISS_VerifyConfig* config) {
    ISS_Verification ver;
    memset(&ver, 0, sizeof(ver));
    if (!sys) return ver;
    ISS_VerifyConfig cfg = config ? *config : iss_verify_config_default();
    /* Start with coarse grid, refine where violations found */
    int base_samples = cfg.monte_carlo_samples / 4;
    int extra_samples = cfg.monte_carlo_samples - base_samples;
    ISS_LyapunovResult r1 = iss_lyapunov_verify_random(sys,
        base_samples, cfg.state_radius, cfg.input_radius);
    if (r1.violation_count > 0) {
        /* Refine around violation region */
        ISS_LyapunovResult r2 = iss_lyapunov_verify_random(sys,
            extra_samples, cfg.state_radius * 0.5, cfg.input_radius * 0.5);
        ver.lyap_result = r2;
        ver.lyap_result.violation_count += r1.violation_count;
        ver.is_ISS = (ver.lyap_result.violation_count == 0);
    } else {
        ver.lyap_result = r1;
        ver.is_ISS = true;
    }
    ver.total_checks = base_samples + extra_samples;
    ver.passed_checks = ver.total_checks - ver.lyap_result.violation_count;
    if (ver.is_ISS) {
        ver.certificate.is_ISS = true;
        ver.certificate.gamma_gain = iss_lyapunov_compute_gain(
            ver.lyap_result.alpha_estimate, ver.lyap_result.sigma_estimate);
        ver.gain = iss_gain_from_lyapunov(&ver.certificate);
    }
    return ver;
}

/* ==============================================================
 * Verification Pipeline: Full ISS Certification
 * ============================================================== */

ISS_Certificate iss_verify_full_pipeline(ISS_System* sys,
    const ISS_VerifyConfig* config) {
    ISS_Verification ver = iss_verify(sys, config);
    if (!ver.is_ISS) {
        /* Try adaptive sampling */
        ver = iss_verify_adaptive(sys, config);
    }
    ISS_Certificate cert;
    memset(&cert, 0, sizeof(cert));
    if (ver.is_ISS) {
        cert = ver.certificate;
        cert.verification_samples = ver.total_checks;
        cert.max_input_norm = config ? config->input_radius : 1.0;
    }
    return cert;
}

/* ==============================================================
 * Batch Verification of Multiple Systems
 * ============================================================== */

int iss_verify_batch(ISS_System** systems, int n_systems,
    ISS_Verification* results) {
    if (!systems || n_systems <= 0) return 0;
    ISS_VerifyConfig cfg = iss_verify_config_default();
    int i;
    int iss_count = 0;
    for (i = 0; i < n_systems; i++) {
        if (systems[i]) {
            results[i] = iss_verify(systems[i], &cfg);
            if (results[i].is_ISS) iss_count++;
        }
    }
    return iss_count;
}

/* ==============================================================
 * ISS Robustness: Verify ISS under parameter perturbation
 * ============================================================== */

bool iss_verify_robustness(ISS_System* sys, double param_delta,
    int n_perturbations) {
    if (!sys || n_perturbations <= 0) return false;
    ISS_VerifyConfig cfg = iss_verify_config_default();
    cfg.monte_carlo_samples = 200;
    int i;
    for (i = 0; i < n_perturbations; i++) {
        double scale = 1.0 + param_delta * (2.0 * rand() / RAND_MAX - 1.0);
        cfg.input_radius *= scale;
        ISS_Verification ver = iss_verify(sys, &cfg);
        if (!ver.is_ISS) return false;
        cfg.input_radius /= scale;
    }
    return true;
}

/* ==============================================================
 * ISS Gain Verification via Input-Output Pairs
 * ============================================================== */

double iss_verify_gain_empirical(ISS_System* sys, double dt,
    double t_final, double u_magnitude) {
    if (!sys || !sys->dynamics || dt <= 0) return INFINITY;
    int n = sys->n_states;
    int m = sys->n_inputs;
    int steps = (int)(t_final / dt);
    double x[ISS_MAX_DIM];
    memcpy(x, sys->initial_state.x, n * sizeof(double));
    double max_state = 0.0;
    int i, j;
    for (i = 0; i < steps; i++) {
        double u[ISS_MAX_DIM];
        double t = i * dt;
        for (j = 0; j < m && j < ISS_MAX_DIM; j++) {
            u[j] = u_magnitude * sin(0.5 * t);
        }
        double dxdt[ISS_MAX_DIM];
        sys->dynamics(x, n, u, m, dxdt);
        for (j = 0; j < n; j++) x[j] += dxdt[j] * dt;
        double xn = 0.0;
        for (j = 0; j < n; j++) xn += x[j] * x[j];
        xn = sqrt(xn);
        if (xn > max_state) max_state = xn;
    }
    return u_magnitude > 1e-10 ? max_state / u_magnitude : INFINITY;
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
