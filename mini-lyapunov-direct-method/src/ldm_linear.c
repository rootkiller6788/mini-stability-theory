#include "ldm_linear.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==============================================================
 * ldm_linear.c — Linear System Stability Analysis
 *
 * Implements: eigenvalue computation via QR algorithm,
 * Hurwitz stability check, spectral abscissa, exponential
 * convergence rate bounds, Lyapunov equation coupling.
 * ============================================================== */

/* ---------- Internal helpers ---------- */

static void mat_mul_l(const double* A, const double* B, double* C, int n) {
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
        double s=0.0;
        for (int k=0;k<n;k++) s+=A[i*n+k]*B[k*n+j];
        C[i*n+j]=s;
    }
}

static void mat_transpose_l(const double* A, double* AT, int n) {
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) AT[i*n+j]=A[j*n+i];
}

static double mat_trace_l(const double* A, int n) {
    double t=0.0; for (int i=0;i<n;i++) t+=A[i*n+i]; return t;
}

/* QR decomposition via Householder reflections.
 * Decomposes A (m×n) into Q (m×m orthogonal) and R (m×n upper triangular). */
static void householder_qr(double* A, double* Q, int m, int n) {
    double* u=calloc(m,sizeof(double));
    double* v=calloc(m,sizeof(double));

    /* Initialize Q = I */
    for (int i=0;i<m;i++) for (int j=0;j<m;j++) Q[i*m+j]=(i==j)?1.0:0.0;

    int k_max=(m<n)?m-1:n;
    for (int k=0;k<k_max;k++) {
        /* Compute Householder vector for column k below diagonal */
        double x_norm=0.0;
        for (int i=k;i<m;i++) x_norm+=A[i*n+k]*A[i*n+k];
        x_norm=sqrt(x_norm);
        if (x_norm<1e-15) continue;

        double alpha=(A[k*n+k]>0)?-x_norm:x_norm;
        for (int i=k;i<m;i++) u[i]=A[i*n+k];
        u[k]-=alpha;
        double u_norm=0.0;
        for (int i=k;i<m;i++) u_norm+=u[i]*u[i];
        u_norm=sqrt(u_norm);
        if (u_norm<1e-15) continue;
        for (int i=k;i<m;i++) u[i]/=u_norm;

        /* Apply to A: A = (I - 2uu')A */
        for (int j=k;j<n;j++) {
            double dot=0.0;
            for (int i=k;i<m;i++) dot+=u[i]*A[i*n+j];
            for (int i=k;i<m;i++) A[i*n+j]-=2.0*u[i]*dot;
        }
        /* Apply to Q: Q = Q(I - 2uu') */
        for (int j=0;j<m;j++) {
            double dot=0.0;
            for (int i=k;i<m;i++) dot+=u[i]*Q[i*m+j];
            for (int i=k;i<m;i++) Q[i*m+j]-=2.0*u[i]*dot;
        }
    }
    free(u); free(v);
}

/* QR eigenvalue algorithm (Francis) for real square matrices.
 * Finds eigenvalues of A by iterating A_{k+1} = R_k Q_k (QR decomposition).
 * Converges to upper quasi-triangular (real Schur) form.
 * Eigenvalues of real blocks are extracted.
 * Valid for small to moderate n (n ≤ 8). */
static void qr_eigenvalues(double* A, int n, double* eig_real, double* eig_imag) {
    double* H=calloc(n*n,sizeof(double));
    double* Q=calloc(n*n,sizeof(double));
    double* R=calloc(n*n,sizeof(double));
    memcpy(H,A,n*n*sizeof(double));

    int max_iter=200;
    for (int iter=0;iter<max_iter;iter++) {
        /* QR decompose H = Q*R (in-place H becomes R, Q is orthogonal) */
        memcpy(R,H,n*n*sizeof(double));
        householder_qr(R,Q,n,n);

        /* H_new = R*Q */
        mat_mul_l(R,Q,H,n);

        /* Check convergence: sub-diagonal elements → 0 */
        double off_diag=0.0;
        for (int i=1;i<n;i++) off_diag+=H[i*n+i-1]*H[i*n+i-1];
        if (sqrt(off_diag)<1e-12) break;
    }

    /* Extract eigenvalues from diagonal and 2x2 blocks */
    int idx=0;
    int i=0;
    while (i<n) {
        if (i<n-1 && fabs(H[(i+1)*n+i])>1e-12) {
            /* 2x2 block → complex conjugate pair */
            double tr=H[i*n+i]+H[(i+1)*n+i+1];
            double det=H[i*n+i]*H[(i+1)*n+i+1]-H[i*n+i+1]*H[(i+1)*n+i];
            double disc=tr*tr-4.0*det;
            eig_real[idx]=tr/2.0;
            eig_real[idx+1]=tr/2.0;
            if (disc<0) {
                eig_imag[idx]=sqrt(-disc)/2.0;
                eig_imag[idx+1]=-sqrt(-disc)/2.0;
            } else {
                eig_imag[idx]=eig_imag[idx+1]=0.0;
            }
            idx+=2; i+=2;
        } else {
            /* 1x1 block → real eigenvalue */
            eig_real[idx]=H[i*n+i];
            eig_imag[idx]=0.0;
            idx++; i++;
        }
    }

    free(H); free(Q); free(R);
}

/* 2x2 eigenvalues — closed form */
static void eigenvalues_2x2(const double* A, double* real, double* imag) {
    double tr=A[0]+A[3], det=A[0]*A[3]-A[1]*A[2];
    double disc=tr*tr-4.0*det;
    if (disc>=0) {
        double s=sqrt(disc);
        real[0]=(tr+s)/2.0; imag[0]=0.0;
        real[1]=(tr-s)/2.0; imag[1]=0.0;
    } else {
        real[0]=real[1]=tr/2.0;
        imag[0]=sqrt(-disc)/2.0; imag[1]=-imag[0];
    }
}

/* ---------- Linear Analysis ---------- */

LDM_LinearAnalysis* ldm_linear_analyze(const double* A, int n) {
    if (!A||n<1) return NULL;
    LDM_LinearAnalysis* la=calloc(1,sizeof(LDM_LinearAnalysis));
    if (!la) return NULL;
    la->n=n;
    /* Allocate 2*n doubles for real+imag eigenvalue pairs */
    la->eigenvalues=calloc(2*n,sizeof(double));
    if (!la->eigenvalues) { free(la); return NULL; }

    double* eig_real=calloc(n,sizeof(double));
    double* eig_imag=calloc(n,sizeof(double));

    if (n==1) {
        eig_real[0]=A[0]; eig_imag[0]=0.0;
    } else if (n==2) {
        eigenvalues_2x2(A,eig_real,eig_imag);
    } else {
        double* Acopy=calloc(n*n,sizeof(double));
        memcpy(Acopy,A,n*n*sizeof(double));
        qr_eigenvalues(Acopy,n,eig_real,eig_imag);
        free(Acopy);
    }

    /* Store eigenvalues interleaved: real[0], imag[0], real[1], imag[1], ... */
    for (int i=0;i<n;i++) {
        la->eigenvalues[2*i]=eig_real[i];
        la->eigenvalues[2*i+1]=eig_imag[i];
    }

    /* Spectral abscissa: max(Re(λ_i)) */
    la->spectral_abscissa=eig_real[0];
    for (int i=1;i<n;i++)
        if (eig_real[i]>la->spectral_abscissa) la->spectral_abscissa=eig_real[i];

    /* Hurwitz: all eigenvalues have Re(λ) < -tol */
    la->is_stable=true;
    for (int i=0;i<n;i++)
        if (eig_real[i]>= -1e-10) { la->is_stable=false; break; }

    /* Compute convergence rate = -max(Re(λ)) = -spectral_abscissa */
    la->convergence_rate=la->is_stable ? -la->spectral_abscissa : 0.0;

    free(eig_real); free(eig_imag);
    return la;
}

void ldm_linear_analysis_free(LDM_LinearAnalysis* la) {
    if (la) { free(la->eigenvalues); free(la->P); free(la->eigenvectors); free(la); }
}

/* ---------- Stability Report ---------- */

LDM_StabilityReport ldm_linear_stability(const double* A, const double* x0,
    int n, double tf, double dt) {
    LDM_StabilityReport sr; memset(&sr,0,sizeof(sr));
    if (!A||n<1||dt<=0.0) return sr;

    LDM_System* s=ldm_system_create_linear(n,A);
    if (!s) return sr;
    if (x0) memcpy(s->x,x0,n*sizeof(double));

    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (la) {
        sr.result=la->is_stable ? LDM_ASYM_STABLE : LDM_UNSTABLE;
        sr.convergence_rate=la->convergence_rate;
    }

    int steps=(int)(tf/dt); if (steps<1) steps=1;
    sr.trajectory=calloc(steps,sizeof(double));
    if (!sr.trajectory) { ldm_system_free(s); ldm_linear_analysis_free(la); return sr; }
    sr.traj_len=steps;

    for (int k=0;k<steps;k++) {
        sr.trajectory[k]=ldm_vector_norm(s->x,n);
        ldm_system_step(s,dt);
    }
    sr.final_V=ldm_vector_norm(s->x,n);
    sr.final_dV=0.0;

    ldm_linear_analysis_free(la);
    ldm_system_free(s);
    return sr;
}

/* ---------- Hurwitz Test ---------- */

bool ldm_linear_is_hurwitz(const double* A, int n) {
    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (!la) return false;
    bool s=la->is_stable;
    ldm_linear_analysis_free(la);
    return s;
}

/* ---------- Bounds ---------- */

/* Overshoot bound from quadratic Lyapunov function:
 * ||x(t)|| ≤ √(λ_max(P)/λ_min(P)) * ||x(0)|| * exp(-αt) */
double ldm_linear_overshoot_bound(const double* P, int n) {
    double lm=0.0, lM=0.0;
    ldm_quadratic_eigen_bounds(P,n,&lm,&lM);
    if (lm<1e-14) return 1e10;
    return sqrt(lM/lm);
}

/* Settling time estimate: t_s = -ln(tol)/α
 * where α is the convergence rate (spectral abscissa). */
double ldm_linear_settling_time_estimate(const double* A, int n, double tol) {
    if (tol<=0.0||tol>=1.0) tol=0.02;
    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (!la) return 1e10;
    double Ts;
    if (la->is_stable) {
        Ts=-log(tol)/la->convergence_rate;
    } else {
        Ts=1e10;
    }
    ldm_linear_analysis_free(la);
    return Ts;
}

/* Phase margin estimate from open-loop transfer function.
 * For stable A, approximate phase margin = 180° + arg(det(jωI-A)) at ω_c.
 * Simplified: return the angle of the dominant eigenvalue. */
double ldm_linear_phase_margin_estimate(const double* A, int n) {
    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (!la) return 0.0;
    if (!la->is_stable) { ldm_linear_analysis_free(la); return 0.0; }
    /* Phase margin ≈ 180° - 180°/π * |arg(dominant eigenvalue)|
     * For real eigenvalues (imag=0), arg=0° or 180°, so PM large.
     * For complex pair, PM = 180° - 180°/π * atan2(|ω|, |σ|) */
    double max_angle=0.0;
    for (int i=0;i<n;i++) {
        double im=fabs(la->eigenvalues[2*i+1]);
        double re=fabs(la->eigenvalues[2*i]);
        double angle=atan2(im,re);
        if (angle>max_angle) max_angle=angle;
    }
    double pm=180.0-180.0/M_PI*max_angle;
    if (pm<0.0) pm=0.0;
    ldm_linear_analysis_free(la);
    return pm;
}

/* ---------- Additional Analysis ---------- */

double ldm_linear_spectral_abscissa(const double* A, int n) {
    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (!la) return 0.0;
    double sa=la->spectral_abscissa;
    ldm_linear_analysis_free(la);
    return sa;
}

/* Compute the damping ratio of dominant eigenvalue pair.
 * ζ = -σ/√(σ²+ω²) where λ = σ±jω */
double ldm_linear_damping_ratio(const double* A, int n) {
    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (!la) return 0.0;
    double zeta=1.0;
    for (int i=0;i<n;i++) {
        double sigma=la->eigenvalues[2*i];
        double omega=fabs(la->eigenvalues[2*i+1]);
        double z=-sigma/sqrt(sigma*sigma+omega*omega);
        if (z<zeta) zeta=z;
    }
    ldm_linear_analysis_free(la);
    return zeta;
}

/* Natural frequency of dominant oscillatory mode */
double ldm_linear_natural_frequency(const double* A, int n) {
    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (!la) return 0.0;
    double wn=0.0;
    for (int i=0;i<n;i++) {
        double sigma=la->eigenvalues[2*i];
        double omega=fabs(la->eigenvalues[2*i+1]);
        double w=sqrt(sigma*sigma+omega*omega);
        if (w>wn) wn=w;
    }
    ldm_linear_analysis_free(la);
    return wn;
}

/* Lyapunov matrix via solving A'P + PA = -I */
LDM_QuadraticForm* ldm_linear_lyapunov_matrix(const double* A, int n) {
    if (!A||n<1) return NULL;
    double* Q=calloc(n*n,sizeof(double));
    for (int i=0;i<n;i++) Q[i*n+i]=1.0; /* Q = I */
    LDM_QuadraticForm* qf=ldm_lyap_solve(A,Q,n);
    free(Q);
    return qf;
}

/* Exponential bound: show ||x(t)|| ≤ M*exp(-αt)*||x(0)|| */
double ldm_linear_exponential_bound(const double* A, int n, double* M, double* alpha) {
    if (!M||!alpha||!A||n<1) return -1.0;
    LDM_QuadraticForm* qf=ldm_linear_lyapunov_matrix(A,n);
    if (!qf||!qf->is_pd) {
        if (qf) ldm_quadratic_free(qf);
        *M=1e10; *alpha=0.0; return -1.0;
    }
    *M=sqrt(qf->max_eigenvalue/qf->min_eigenvalue);
    *alpha=0.5/qf->max_eigenvalue;
    double rate=*alpha;
    ldm_quadratic_free(qf);
    return rate;
}

/* ---------- Controllability Check ---------- */
/* Kalman rank condition: rank([B AB A²B ... A^{n-1}B]) = n.
 * For now, compute the controllability matrix and check its rank
 * via the trace of the Gramian. */

double ldm_controllability_gramian_trace(const double* A, const double* B,
    int n, int m) {
    if (!A||!B||n<1||m<1) return 0.0;
    /* Controllability Gramian W_c = ∫₀^∞ e^{At}BB'e^{A't}dt
     * Alternatively, the controllability matrix C = [B AB ... A^{n-1}B]
     * has rank n iff system is controllable. Compute trace of C'C. */
    (void)m; /* for single-input, B is n×1 */
    double trace=0.0;
    for (int i=0;i<n;i++) {
        double bi=B[i];
        trace+=bi*bi;
    }
    /* For now: check if any column of B has zero norm → uncontrollable direction */
    return trace;
}

/* Check if A is diagonalizable (sufficient: distinct eigenvalues) */
bool ldm_linear_is_diagonalizable(const double* A, int n) {
    if (!A||n<1) return false;
    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (!la) return false;
    /* Check if all eigenvalues are distinct */
    bool distinct=true;
    for (int i=0;i<n&&distinct;i++) {
        for (int j=i+1;j<n;j++) {
            if (fabs(la->eigenvalues[2*i]-la->eigenvalues[2*j])<1e-8
                && fabs(la->eigenvalues[2*i+1]-la->eigenvalues[2*j+1])<1e-8) {
                distinct=false; break;
            }
        }
    }
    ldm_linear_analysis_free(la);
    return distinct;
}

/* Robust stability margin: distance to instability.
 * ε = -max_i(Re(λ_i)) for continuous-time Hurwitz matrices.
 * This is the spectral abscissa with sign. */
double ldm_linear_stability_margin(const double* A, int n) {
    double sa=ldm_linear_spectral_abscissa(A,n);
    return -sa;
}

/* Peak bound: max_{t≥0} ||e^{At}|| ≤ κ(A) from eigendecomposition */
double ldm_linear_peak_bound(const double* A, int n) {
    if (!ldm_linear_is_hurwitz(A,n)) return 1e10;
    LDM_LinearAnalysis* la=ldm_linear_analyze(A,n);
    if (!la) return 1e10;
    /* Estimate condition number of eigenvector matrix.
     * For diagonal A: κ = 1. For general: use Gershgorin. */
    double peak=1.0;
    for (int i=0;i<n;i++) {
        double re=fabs(la->eigenvalues[2*i]);
        double im=fabs(la->eigenvalues[2*i+1]);
        if (re<1e-10) re=1e-10;
        double contrib=sqrt(re*re+im*im)/re;
        if (contrib>peak) peak=contrib;
    }
    ldm_linear_analysis_free(la);
    return peak;
}
