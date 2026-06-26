#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "circle_criterion.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
TransferFunction* tf_create(int n_num, int n_den) {
    TransferFunction* tf = (TransferFunction*)calloc(1, sizeof(TransferFunction));
    if (!tf) return NULL;
    tf->n_num = n_num; tf->n_den = n_den;
    tf->num = (double*)calloc((size_t)n_num, sizeof(double));
    tf->den = (double*)calloc((size_t)n_den, sizeof(double));
    if (!tf->num || !tf->den) { free(tf->num); free(tf->den); free(tf); return NULL; }
    if (n_den > 0) tf->den[0] = 1.0;
    return tf;
}
void tf_free(TransferFunction* tf) {
    if (!tf) return;
    free(tf->num); free(tf->den); free(tf);
}
void tf_set_coeffs(TransferFunction* tf, const double *num, const double *den) {
    if (!tf || !num || !den) return;
    for (int i = 0; i < tf->n_num; i++) tf->num[i] = num[i];
    for (int i = 0; i < tf->n_den; i++) tf->den[i] = den[i];
}

/* Complex evaluation G(j*omega) = N(jw)/D(jw)
 * (jw)^k = w^k * j^k: k%4=0->1, k%4=1->j, k%4=2->-1, k%4=3->-j */
void tf_eval_complex(const TransferFunction* tf, double omega,
                     double *re, double *im) {
    *re = 0.0; *im = 0.0;
    if (!tf || tf->n_den < 1) return;
    double Nr = 0.0, Ni = 0.0, wp = 1.0;
    for (int k = 0; k < tf->n_num; k++) {
        double c = tf->num[k];
        if (c == 0.0) { wp *= omega; continue; }
        int r = k & 3;
        if      (r == 0) Nr += c * wp;
        else if (r == 1) Ni += c * wp;
        else if (r == 2) Nr -= c * wp;
        else             Ni -= c * wp;
        wp *= omega;
    }
    double Dr = 0.0, Di = 0.0; wp = 1.0;
    for (int k = 0; k < tf->n_den; k++) {
        double c = tf->den[k];
        if (c == 0.0) { wp *= omega; continue; }
        int r = k & 3;
        if      (r == 0) Dr += c * wp;
        else if (r == 1) Di += c * wp;
        else if (r == 2) Dr -= c * wp;
        else             Di -= c * wp;
        wp *= omega;
    }
    double denom = Dr * Dr + Di * Di;
    if (denom < 1e-30) {
        *re = (Nr >= 0.0) ? 1e15 : -1e15;
        *im = (Ni >= 0.0) ? 1e15 : -1e15;
        return;
    }
    *re = (Nr * Dr + Ni * Di) / denom;
    *im = (Ni * Dr - Nr * Di) / denom;
}

double tf_eval(const TransferFunction* tf, double omega) {
    double re, im; tf_eval_complex(tf, omega, &re, &im);
    return sqrt(re * re + im * im);
}

double tf_phase(const TransferFunction* tf, double omega) {
    double re, im; tf_eval_complex(tf, omega, &re, &im);
    return atan2(im, re);
}

/* ===== Nyquist Data Management ===== */
NyquistData* nyquist_create(int n) {
    NyquistData* nd = (NyquistData*)calloc(1, sizeof(NyquistData));
    if (!nd) return NULL;
    nd->n = n;
    nd->real = (double*)calloc((size_t)n, sizeof(double));
    nd->imag = (double*)calloc((size_t)n, sizeof(double));
    if (!nd->real || !nd->imag) { free(nd->real); free(nd->imag); free(nd); return NULL; }
    return nd;
}
void nyquist_free(NyquistData* nd) {
    if (!nd) return;
    free(nd->real); free(nd->imag); free(nd);
}
void nyquist_compute(const TransferFunction* tf, NyquistData* nd,
                     double wmin, double wmax) {
    if (!tf || !nd || nd->n < 2) return;
    double dw = (wmax - wmin) / (double)(nd->n - 1);
    for (int i = 0; i < nd->n; i++) {
        double w = wmin + (double)i * dw;
        tf_eval_complex(tf, w, &nd->real[i], &nd->imag[i]);
    }
}
void nyquist_compute_log(const TransferFunction* tf, NyquistData* nd,
                          double wmin, double wmax) {
    if (!tf || !nd || nd->n < 2) return;
    if (wmin < 1e-8) wmin = 1e-8;
    double lmin = log10(wmin), lmax = log10(wmax);
    double dl = (lmax - lmin) / (double)(nd->n - 1);
    for (int i = 0; i < nd->n; i++) {
        double w = pow(10.0, lmin + (double)i * dl);
        tf_eval_complex(tf, w, &nd->real[i], &nd->imag[i]);
    }
}

/* ===== Circle Criterion: Critical Disk Formulae =====
 * Sector [alpha,beta]: center c = -0.5*(1/alpha+1/beta), radius r = 0.5*|1/alpha-1/beta|
 * alpha=0 => half-plane; alpha=beta => point (linear) */
double circle_criterion_center(SectorBounds sb) {
    if (fabs(sb.alpha) < CIRCLE_EPS) return -1.0 / sb.beta;
    if (fabs(sb.beta) < CIRCLE_EPS)  return -1.0 / sb.alpha;
    return -0.5 * (1.0 / sb.alpha + 1.0 / sb.beta);
}
double circle_criterion_radius(SectorBounds sb) {
    if (fabs(sb.alpha) < CIRCLE_EPS || fabs(sb.beta) < CIRCLE_EPS) return 1e12;
    if (fabs(sb.alpha - sb.beta) < CIRCLE_EPS) return 0.0;
    return 0.5 * fabs(1.0 / sb.alpha - 1.0 / sb.beta);
}

/* ===== Winding Number / Encirclement Detection =====
 * W = (1/2pi) * sum Delta_angle_i. Non-zero W => encirclement. */
bool circle_criterion_encircles(const NyquistData* nd, double cx, double cy) {
    if (!nd || nd->n < 3) return false;
    double total = 0.0;
    for (int i = 1; i < nd->n; i++) {
        double dx0 = nd->real[i-1] - cx, dy0 = nd->imag[i-1] - cy;
        double dx1 = nd->real[i]   - cx, dy1 = nd->imag[i]   - cy;
        if ((fabs(dx0) < 1e-14 && fabs(dy0) < 1e-14) ||
            (fabs(dx1) < 1e-14 && fabs(dy1) < 1e-14)) continue;
        double a0 = atan2(dy0, dx0), a1 = atan2(dy1, dx1), d = a1 - a0;
        if (d >  M_PI) d -= 2.0 * M_PI;
        if (d < -M_PI) d += 2.0 * M_PI;
        total += d;
    }
    int windings = (int)round(total / (2.0 * M_PI));
    return windings != 0;
}

/* ===== Distance and Disk Intersection ===== */
double circle_criterion_min_distance(const NyquistData* nd,
                                      double center_x, double center_y) {
    if (!nd || nd->n < 1) return 1e9;
    double best = 1e30;
    for (int i = 0; i < nd->n; i++) {
        double dx = nd->real[i] - center_x, dy = nd->imag[i] - center_y;
        double d = sqrt(dx * dx + dy * dy);
        if (d < best) best = d;
    }
    return best;
}

bool circle_criterion_in_disk(const NyquistData* nd, double center_x, double radius) {
    if (!nd || nd->n < 1) return false;
    double r2 = radius * radius;
    for (int i = 0; i < nd->n; i++) {
        double dx = nd->real[i] - center_x;
        double dy = nd->imag[i];  /* center_y = 0 (disk on real axis) */
        if (dx * dx + dy * dy <= r2 + CIRCLE_EPS) return true;
    }
    return false;
}

/* ===== Main Circle Criterion Check =====
 * Returns 1=ABSOLUTELY STABLE, 0=INCONCLUSIVE.
 * *margin = r - min_dist (>0 => stable with safety buffer). */
int circle_criterion_check(const TransferFunction* tf, SectorBounds sb,
                            double *margin) {
    *margin = 0.0;
    if (!tf) return 0;
    if (fabs(sb.alpha) < CIRCLE_EPS)
        return circle_criterion_halfplane_case(tf, sb.beta, margin);
    if (fabs(sb.alpha - sb.beta) < CIRCLE_EPS)
        return circle_criterion_linear_case(tf, sb.alpha, margin);
    double cx = circle_criterion_center(sb);
    double cr = circle_criterion_radius(sb);
    if (cr > 1e10) return 0;
    NyquistData* nd = nyquist_create(CIRCLE_N_POINTS);
    if (!nd) return 0;
    nyquist_compute_log(tf, nd, 0.01, CIRCLE_MAX_FREQ);
    if (circle_criterion_encircles(nd, cx, 0.0)) {
        nyquist_free(nd); *margin = -1.0; return 0;
    }
    double min_d = circle_criterion_min_distance(nd, cx, 0.0);
    *margin = min_d - cr;
    nyquist_free(nd);
    return (*margin > CIRCLE_EPS) ? 1 : 0;
}

/* Half-plane case: alpha=0 => Re[G(jw)] > -1/beta for all w */
int circle_criterion_halfplane_case(const TransferFunction* tf,
                                     double beta, double *margin) {
    *margin = 0.0;
    if (!tf || beta <= CIRCLE_EPS) return 0;
    double target = -1.0 / beta, worst_re = 1e30;
    for (int i = 0; i < CIRCLE_N_POINTS; i++) {
        double w = 0.01 * pow(10.0, 3.0 * (double)i / (double)(CIRCLE_N_POINTS - 1));
        double re, im; tf_eval_complex(tf, w, &re, &im); (void)im;
        if (re < worst_re) worst_re = re;
    }
    { double re, im; tf_eval_complex(tf, 1e-4, &re, &im); if (re < worst_re) worst_re = re; }
    { double re, im; tf_eval_complex(tf, CIRCLE_MAX_FREQ, &re, &im); if (re < worst_re) worst_re = re; }
    *margin = worst_re - target;
    return (*margin > CIRCLE_EPS) ? 1 : 0;
}

/* Linear case: alpha=beta=gain => classical Nyquist at -1/gain */
int circle_criterion_linear_case(const TransferFunction* tf,
                                  double gain, double *margin) {
    *margin = 0.0;
    if (!tf || fabs(gain) < CIRCLE_EPS) return 0;
    double crit = -1.0 / gain;
    NyquistData* nd = nyquist_create(CIRCLE_N_POINTS);
    if (!nd) return 0;
    nyquist_compute_log(tf, nd, 0.01, CIRCLE_MAX_FREQ);
    if (circle_criterion_encircles(nd, crit, 0.0)) {
        nyquist_free(nd); *margin = -1.0; return 0;
    }
    double min_d = circle_criterion_min_distance(nd, crit, 0.0);
    *margin = min_d; nyquist_free(nd);
    return (min_d > CIRCLE_EPS) ? 1 : 0;
}

/* ===== Plot Data Generation ===== */
void circle_criterion_plot_data(const TransferFunction* tf,
                                 SectorBounds sb, double *x, double *y, int n) {
    (void)sb;
    if (!tf || !x || !y || n < 2) return;
    NyquistData* nd = nyquist_create(n);
    if (!nd) return;
    nyquist_compute(tf, nd, 0.01, CIRCLE_MAX_FREQ);
    for (int i = 0; i < n; i++) { x[i] = nd->real[i]; y[i] = nd->imag[i]; }
    nyquist_free(nd);
}

/* ===== Custom Frequency Range Check ===== */
int circle_criterion_check_freqs(const TransferFunction* tf,
                                  SectorBounds sb, double *margin,
                                  double wmin, double wmax, int npts) {
    *margin = 0.0;
    if (!tf || npts < 3) return 0;
    double cx = circle_criterion_center(sb), cr = circle_criterion_radius(sb);
    if (cr > 1e10) return 0;
    NyquistData* nd = nyquist_create(npts);
    if (!nd) return 0;
    nyquist_compute(tf, nd, wmin, wmax);
    if (circle_criterion_encircles(nd, cx, 0.0)) {
        nyquist_free(nd); *margin = -1.0; return 0;
    }
    double min_d = circle_criterion_min_distance(nd, cx, 0.0);
    *margin = min_d - cr; nyquist_free(nd);
    return (*margin > CIRCLE_EPS) ? 1 : 0;
}

/* ===== Gain Crossover Frequency ===== */
double circle_criterion_crossover_freq(const TransferFunction* tf) {
    if (!tf) return 0.0;
    double best_w = 0.0, best_e = 1e30;
    for (int i = 0; i < CIRCLE_N_POINTS; i++) {
        double w = 0.01 * pow(10.0, 3.0 * (double)i / (double)(CIRCLE_N_POINTS - 1));
        double re, im; tf_eval_complex(tf, w, &re, &im);
        double mag = sqrt(re * re + im * im);
        double err = fabs(mag - 1.0);
        if (err < best_e) { best_e = err; best_w = w; }
    }
    return best_w;
}

/* ===== Self-Test Routine =====
 * Validates internal consistency. Uses assert() for all checks.
 * Returns 1 on success, aborts on failure. */
int circle_criterion_self_test(void) {
    /* G(s) = 1/(s+1) */
    TransferFunction* tf = tf_create(1, 2);
    assert(tf != NULL);
    tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;

    /* G(j*1) = 1/(1+j) = (1-j)/2 = 0.5 - j0.5 */
    double re, im;
    tf_eval_complex(tf, 1.0, &re, &im);
    assert(fabs(re - 0.5) < 0.01);
    assert(fabs(im + 0.5) < 0.01);

    /* DC gain = 1.0, HF gain -> 0 */
    assert(fabs(tf_eval(tf, 0.0) - 1.0) < 0.01);
    assert(tf_eval(tf, 1e6) < 0.01);

    /* Phase at low frequency: ~0 (slightly negative) */
    assert(tf_phase(tf, 0.01) < 0.0);
    assert(tf_phase(tf, 0.01) > -0.1);

    /* Critical disk for sector [0.5, 2.0]:
     *   center = -0.5*(2.0 + 0.5) = -1.25
     *   radius = 0.5*|2.0 - 0.5|  = 0.75 */
    SectorBounds sb = {0.5, 2.0, 0, 0, 0.0, 0.0};
    assert(fabs(circle_criterion_center(sb) + 1.25) < 0.01);
    assert(fabs(circle_criterion_radius(sb) - 0.75) < 0.01);

    /* Half-plane [0, 2]: Re[G(jw)]=1/(1+w^2) >= 0 > -0.5 => PASS */
    SectorBounds sb_hp = {0.0, 2.0, 0, 0, 0.0, 0.0};
    double margin;
    int result = circle_criterion_check(tf, sb_hp, &margin);
    assert(result == 1);
    assert(margin > 0.0);

    /* Linear [2, 2]: Nyquist avoids -0.5 point => PASS */
    SectorBounds sb_lin = {2.0, 2.0, 0, 0, 0.0, 0.0};
    result = circle_criterion_check(tf, sb_lin, &margin);
    assert(result == 1);
    assert(margin > 0.0);

    /* General [1, 10]: check runs without crash */
    SectorBounds sb_gen = {1.0, 10.0, 0, 0, 0.0, 0.0};
    result = circle_criterion_check(tf, sb_gen, &margin);
    (void)result; (void)margin;

    /* Nyquist data tests */
    NyquistData* nd = nyquist_create(100);
    assert(nd != NULL);
    nyquist_compute(tf, nd, 0.01, 10.0);
    assert(nd->n == 100);
    assert(!circle_criterion_encircles(nd, -10.0, 0.0));
    assert(circle_criterion_min_distance(nd, -5.0, 0.0) > 1.0);
    assert(!circle_criterion_in_disk(nd, -100.0, 0.5));

    /* Log-swept Nyquist */
    NyquistData* nd2 = nyquist_create(200);
    assert(nd2 != NULL);
    nyquist_compute_log(tf, nd2, 0.01, 100.0);
    assert(nd2->n == 200);
    assert(fabs(nd2->real[0] - 1.0) < 0.1);

    /* Custom frequency range check */
    double m2;
    int r2 = circle_criterion_check_freqs(tf, sb_hp, &m2, 0.1, 50.0, 500);
    assert(r2 == 1);
    assert(m2 > 0.0);

    /* Crossover frequency should be near 1.0 rad/s for 1/(s+1) */
    double wc = circle_criterion_crossover_freq(tf);
    assert(wc > 0.5 && wc < 2.0);

    /* Set coefficients via convenience function */
    TransferFunction* tf2 = tf_create(2, 3);
    assert(tf2 != NULL);
    double n2[] = {0.0, 1.0};
    double d3[] = {1.0, 0.5, 1.0};
    tf_set_coeffs(tf2, n2, d3);
    assert(tf2->num[0] == 0.0);
    assert(tf2->num[1] == 1.0);
    assert(tf2->den[1] == 0.5);
    assert(tf2->den[2] == 1.0);
    tf_free(tf2);

    nyquist_free(nd);
    nyquist_free(nd2);
    tf_free(tf);
    return 1;
}

/* Additional helper: compute the Nyquist plot with both positive and
 * negative frequency sweeps (for complete Nyquist contour). */
void nyquist_compute_full(const TransferFunction* tf, NyquistData* nd,
                           double wmin, double wmax) {
    if (!tf || !nd || nd->n < 4) return;
    int half = nd->n / 2;
    double dw = (wmax - wmin) / (double)(half - 1);
    /* Negative frequencies (mirror image: G(-jw) = conj(G(jw))) */
    for (int i = 0; i < half; i++) {
        double w = wmax - (double)i * dw;
        tf_eval_complex(tf, w, &nd->real[i], &nd->imag[i]);
        nd->imag[i] = -nd->imag[i];  /* conjugate for -w */
    }
    /* Positive frequencies */
    for (int i = 0; i < half; i++) {
        double w = wmin + (double)i * dw;
        tf_eval_complex(tf, w, &nd->real[half + i], &nd->imag[half + i]);
    }
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
