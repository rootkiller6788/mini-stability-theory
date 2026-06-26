#include "ldm_nonlinear.h"
#include "ldm_linear.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ==============================================================
 * ldm_nonlinear.c �?Nonlinear System Stability
 *
 * Implements:
 *   - Linearization-based local stability (Lyapunov's indirect method)
 *   - Krasovskii's method: V = f'Pf
 *   - Variable gradient method: construct V via curl-free gradient
 *   - Energy-based Lyapunov functions
 *   - LaSalle's invariance principle check
 *   - Estimation of decay rate via simulation
 * ============================================================== */

/* ---------- Internal helpers ---------- */

static void mat_mul_nl(const double* A, const double* B, double* C, int n) {
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
        double s=0.0;
        for (int k=0;k<n;k++) s+=A[i*n+k]*B[k*n+j];
        C[i*n+j]=s;
    }
}

static void mat_vec_mul_nl(const double* A, const double* x, double* y, int n) {
    for (int i=0;i<n;i++) { y[i]=0.0;
        for (int j=0;j<n;j++) y[i]+=A[i*n+j]*x[j]; }
}

static void mat_transpose_nl(const double* A, double* AT, int n) {
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) AT[i*n+j]=A[j*n+i];
}

/* ---------- Jacobian via finite differences ---------- */
/* ∂f_i/∂x_j �?(f_i(x + h*e_j) - f_i(x - h*e_j)) / (2h) */

void ldm_jacobian(const LDM_System* s, const double* x0, double* J, int n) {
    if (!s||!x0||!J||n<1) return;
    double h=1e-6;
    double* xp=calloc(n,sizeof(double));
    double* xm=calloc(n,sizeof(double));
    double* fp=calloc(n,sizeof(double));
    double* fm=calloc(n,sizeof(double));

    for (int j=0;j<n;j++) {
        memcpy(xp,x0,n*sizeof(double));
        memcpy(xm,x0,n*sizeof(double));
        xp[j]+=h; xm[j]-=h;

        /* Evaluate at x+h*e_j */
        if (s->dynamics) s->dynamics(s->t,xp,n,fp);
        else if (s->is_linear&&s->A) mat_vec_mul_nl(s->A,xp,fp,n);

        /* Evaluate at x-h*e_j */
        if (s->dynamics) s->dynamics(s->t,xm,n,fm);
        else if (s->is_linear&&s->A) mat_vec_mul_nl(s->A,xm,fm,n);

        for (int i=0;i<n;i++) J[i*n+j]=(fp[i]-fm[i])/(2.0*h);
    }
    free(xp); free(xm); free(fp); free(fm);
}

/* Simplified Jacobian for linear systems (returns A directly) */
void ldm_jacobian_linear(const double* A, double* J, int n) {
    if (A&&J) memcpy(J,A,n*n*sizeof(double));
}

/* ---------- Stability Analysis ---------- */

LDM_StabilityReport ldm_nonlinear_stability(const LDM_System* s,
    const double* x0, int n, double tf, double dt) {
    LDM_StabilityReport sr; memset(&sr,0,sizeof(sr));
    if (!s||n<1||tf<=0.0||dt<=0.0) return sr;

    int steps=(int)(tf/dt); if (steps<1) steps=1;
    sr.trajectory=calloc(steps,sizeof(double));
    if (!sr.trajectory) return sr;
    sr.traj_len=steps;

    /* Save current state */
    double* x_save=calloc(n,sizeof(double));
    if (x_save) memcpy(x_save,s->x,n*sizeof(double));
    if (x0) memcpy(s->x,x0,n*sizeof(double));

    double max_norm=0.0;
    for (int k=0;k<steps;k++) {
        double norm=ldm_vector_norm(s->x,n);
        sr.trajectory[k]=norm;
        if (norm>max_norm) max_norm=norm;
        ldm_system_step((LDM_System*)s,dt);
    }

    sr.final_V=ldm_vector_norm(s->x,n);
    sr.final_dV=steps>1 ? (sr.trajectory[steps-1]-sr.trajectory[steps-2])/dt : 0.0;

    /* Heuristic stability classification */
    double final_norm=sr.final_V;
    double initial_norm=x0?ldm_vector_norm(x0,n):0.0;

    if (initial_norm<1e-12) {
        sr.result=LDM_STABLE; /* At equilibrium */
    } else if (final_norm<1e-6) {
        sr.result=LDM_ASYM_STABLE;
    } else if (final_norm<0.05*initial_norm) {
        sr.result=LDM_GLOBAL_ASYM;
    } else if (final_norm>10.0*initial_norm) {
        sr.result=LDM_UNSTABLE;
    } else if (final_norm<=initial_norm) {
        sr.result=LDM_STABLE;
    } else {
        sr.result=LDM_UNSTABLE;
    }

    /* Restore state */
    if (x_save) { memcpy(s->x,x_save,n*sizeof(double)); free(x_save); }
    return sr;
}

/* ---------- Krasovskii's Method ---------- */
/* For dx/dt = f(x), try V(x) = f(x)' P f(x).
 * dV/dt = f' (J'P + PJ) f.
 * If J'P + PJ �?-εI < 0 globally, then origin is asymptotically stable.
 * Typically choose P = I, giving V = |f(x)|². */

double ldm_krasovskii_candidate(const LDM_System* s, int n) {
    if (!s||n<1) return 0.0;
    double* dx=calloc(n,sizeof(double));
    if (s->is_linear&&s->A) mat_vec_mul_nl(s->A,s->x,dx,n);
    else if (s->dynamics) s->dynamics(s->t,s->x,n,dx);
    double V=0.0;
    for (int i=0;i<n;i++) V+=dx[i]*dx[i];
    free(dx);
    return V;
}

/* Krasovskii derivative: dV/dt = f' (J'+J) f (with P=I) */
double ldm_krasovskii_derivative(const LDM_System* s, int n) {
    if (!s||n<1) return 0.0;

    /* Compute f(x) */
    double* f=calloc(n,sizeof(double));
    if (s->is_linear&&s->A) mat_vec_mul_nl(s->A,s->x,f,n);
    else if (s->dynamics) s->dynamics(s->t,s->x,n,f);

    /* Compute Jacobian J at x */
    double* J=calloc(n*n,sizeof(double));
    ldm_jacobian(s,s->x,J,n);

    /* Compute JT = J' */
    double* JT=calloc(n*n,sizeof(double));
    mat_transpose_nl(J,JT,n);

    /* Compute J' + J */
    double* JJ=calloc(n*n,sizeof(double));
    for (int i=0;i<n*n;i++) JJ[i]=JT[i]+J[i];

    /* dV/dt = f' (J'+J) f */
    double* Jf=calloc(n,sizeof(double));
    mat_vec_mul_nl(JJ,f,Jf,n);
    double dV=0.0;
    for (int i=0;i<n;i++) dV+=f[i]*Jf[i];

    free(f); free(J); free(JT); free(JJ); free(Jf);
    return dV;
}

/* Full Krasovskii: check if J'P + PJ < 0 at sampled points */
bool ldm_krasovskii_check(const LDM_System* s, int n, int n_samples, double range) {
    if (!s||n<1) return false;

    /* Use P = I for simplicity */
    for (int k=0;k<n_samples;k++) {
        double* x_test=calloc(n,sizeof(double));
        for (int i=0;i<n;i++) x_test[i]=-range+2.0*range*((double)rand()/RAND_MAX);

        /* Save and set state */
        double* x_save=calloc(n,sizeof(double));
        if (x_save) memcpy(x_save,s->x,n*sizeof(double));
        memcpy(s->x,x_test,n*sizeof(double));

        double* J=calloc(n*n,sizeof(double));
        ldm_jacobian(s,s->x,J,n);

        double* JT=calloc(n*n,sizeof(double));
        mat_transpose_nl(J,JT,n);

        /* J' + J (with P=I) */
        double* M=calloc(n*n,sizeof(double));
        for (int i=0;i<n*n;i++) M[i]=JT[i]+J[i];

        /* Check negative definiteness */
        bool stable=ldm_is_negative_definite(M,n);

        /* Restore state */
        if (x_save) memcpy(s->x,x_save,n*sizeof(double));
        free(x_test); free(x_save); free(J); free(JT); free(M);

        if (!stable) return false;
    }
    return true;
}

/* ---------- Variable Gradient Method ---------- */
/* Assume grad(V)(x) = [ Σ a_{ij}(x) x_j ].
 * Enforce curl condition ∂g_i/∂x_j = ∂g_j/∂x_i.
 * Then V(x) = ∫_0^1 g(σx)' x dσ (line integral from origin).
 *
 * For 2D systems, the most common construction:
 *   g1 = a*x1 + b*x2
 *   g2 = c*x1 + d*x2
 * Curl condition: ∂g1/∂x2 = ∂g2/∂x1  �? b = c  (symmetric matrix!)
 *
 * So V(x) = (1/2)x'Px with P = [a b; b d], which is the quadratic form case.
 * General variable gradient: g_i = Σ_j p_{ij}(x) x_j with p_{ij} possibly nonlinear.
 *
 * This function constructs a quadratic Lyapunov function for n�? using
 * the variable gradient approach with constant coefficients.
 */
bool ldm_variable_gradient_construct(const LDM_System* s, int n, double* P_out) {
    if (!s||!P_out||n<1) return false;

    /* Try simple quadratic: P = I */
    for (int i=0;i<n*n;i++) P_out[i]=(i%(n+1)==0)?1.0:0.0;

    if (ldm_is_positive_definite(P_out,n)) {
        /* Check if V_dot < 0 at origin neighborhood */
        double* x_save=calloc(n,sizeof(double));
        if (x_save) memcpy(x_save,s->x,n*sizeof(double));

        bool valid=true;
        double* test_pts[]={
            (double[]){0.1,0.1,0.1},
            (double[]){-0.1,0.2,0.0},
            (double[]){0.2,-0.1,0.0},
            (double[]){-0.1,-0.1,0.1},
        };

        for (int t=0;t<4;t++) {
            memcpy(s->x,test_pts[t],n*sizeof(double));

            double* dx=calloc(n,sizeof(double));
            if (s->is_linear&&s->A) mat_vec_mul_nl(s->A,s->x,dx,n);
            else if (s->dynamics) s->dynamics(s->t,s->x,n,dx);

            /* dV/dt = x'(A'P+PA)x = 2*x'P*dx */
            double* Pdx=calloc(n,sizeof(double));
            mat_vec_mul_nl(P_out,dx,Pdx,n);
            double dV=0.0;
            for (int i=0;i<n;i++) dV+=2.0*s->x[i]*Pdx[i];

            free(dx); free(Pdx);
            if (dV>=0.0) { valid=false; break; }
        }

        if (x_save) memcpy(s->x,x_save,n*sizeof(double));
        free(x_save);
        if (valid) return true;
    }

    /* If I doesn't work, return false (requires more sophisticated construction) */
    return false;
}

/* Variable gradient with specified coefficient form for 2D:
 * g1 = a11*x1 + a12*x2, g2 = a21*x1 + a22*x2
 * Returns V = (1/2)x'Px with P = [a11 a12; a21 a22] symmetric */
bool ldm_variable_gradient_2d(const LDM_System* s, double* P_out) {
    if (!s||!P_out) return false;

    /* Try a family of candidates by varying a parameter */
    double candidates[][4] = {
        {1.0,0.0,0.0,1.0},   /* P = I */
        {2.0,0.0,0.0,1.0},   /* P = diag(2,1) */
        {1.0,0.0,0.0,2.0},   /* P = diag(1,2) */
        {2.0,1.0,1.0,2.0},   /* P = [2 1; 1 2] */
        {3.0,1.0,1.0,1.0},   /* P = [3 1; 1 1] */
    };

    double* x_save=calloc(2,sizeof(double));
    memcpy(x_save,s->x,2*sizeof(double));

    for (int c=0;c<5;c++) {
        memcpy(P_out,candidates[c],4*sizeof(double));

        if (!ldm_is_positive_definite(P_out,2)) continue;

        bool valid=true;
        double test_pts[][2]={{0.1,0.1},{-0.1,0.2},{0.2,-0.1},{-0.1,-0.1},
                              {0.5,0.0},{0.0,0.5},{-0.5,0.3},{0.3,-0.5}};

        for (int t=0;t<8&&valid;t++) {
            s->x[0]=test_pts[t][0]; s->x[1]=test_pts[t][1];

            double dx[2]={0};
            if (s->dynamics) s->dynamics(s->t,s->x,2,dx);

            /* dV/dt = 2*(P11*x1*dx1 + P12*(x1*dx2+x2*dx1) + P22*x2*dx2) */
            double dV=2.0*(P_out[0]*s->x[0]*dx[0]
                         + P_out[1]*(s->x[0]*dx[1]+s->x[1]*dx[0])
                         + P_out[3]*s->x[1]*dx[1]);
            if (dV>=0.0) valid=false;
        }
        if (valid) { free(x_save); return true; }
    }

    memcpy(s->x,x_save,2*sizeof(double));
    free(x_save);
    return false;
}

/* ---------- Energy-Based Lyapunov Function ---------- */
/* For mechanical systems: V = (1/2) Σ v_i² + U(q)
 * This implementation provides the kinetic energy component. */

double ldm_energy_based(const LDM_System* s, int n) {
    if (!s) return 0.0;
    /* Generic: V = (1/2) |x|² (kinetic energy of point mass) */
    double E=0.0;
    for (int i=0;i<n;i++) E+=s->x[i]*s->x[i];
    return 0.5*E;
}

/* Energy derivative: dV/dt = x'f(x) */
double ldm_energy_derivative(const LDM_System* s, int n) {
    if (!s||n<1) return 0.0;
    double* dx=calloc(n,sizeof(double));
    if (s->is_linear&&s->A) mat_vec_mul_nl(s->A,s->x,dx,n);
    else if (s->dynamics) s->dynamics(s->t,s->x,n,dx);
    double dV=0.0;
    for (int i=0;i<n;i++) dV+=s->x[i]*dx[i];
    free(dx);
    return dV;
}

/* ---------- Local / Global Stability Classification ---------- */

bool ldm_is_locally_stable(const LDM_System* s, const double* x_eq, int n) {
    if (!s||n<1) return false;
    /* Linearize at equilibrium and check Hurwitz */
    double* J=calloc(n*n,sizeof(double));
    ldm_jacobian(s,x_eq,J,n);
    bool stable=ldm_linear_is_hurwitz(J,n);
    free(J);
    return stable;
}

bool ldm_is_globally_stable(const LDM_System* s, int n) {
    if (!s||n<1) return false;
    /* Test multiple initial conditions */
    LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
    if (!s2) return false;

    bool global=true;
    double test_pts[][4]={
        {1.0,0.0,0.0,0.0},{-1.0,0.0,0.0,0.0},{0.0,1.0,0.0,0.0},
        {0.0,-1.0,0.0,0.0},{2.0,2.0,0.0,0.0},{-2.0,-2.0,0.0,0.0},
        {5.0,0.0,0.0,0.0},{0.0,5.0,0.0,0.0}};

    int n_test=(n<=4)?8:3;
    for (int t=0;t<n_test&&global;t++) {
        for (int i=0;i<n;i++) s2->x[i]=test_pts[t][i];
        LDM_StabilityReport sr=ldm_nonlinear_stability(s2,s2->x,n,10.0,0.01);
        if (sr.result==LDM_UNSTABLE) global=false;
        free(sr.trajectory);
    }
    ldm_system_free(s2);
    return global;
}

/* ---------- Decay Rate Estimation ---------- */

double ldm_estimate_decay_rate(const LDM_System* s, const double* x0,
    int n, double dt, int steps) {
    if (!s||!x0||n<1||dt<=0.0||steps<2) return 0.0;

    LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
    if (!s2) return 0.0;
    memcpy(s2->x,x0,n*sizeof(double));

    double V0=ldm_vector_norm(s2->x,n);
    if (V0<1e-12) { ldm_system_free(s2); return 0.0; }

    for (int k=0;k<steps;k++) ldm_system_step(s2,dt);
    double Vf=ldm_vector_norm(s2->x,n);

    ldm_system_free(s2);
    if (Vf<1e-12||Vf>=V0) return 0.0;
    return -log(Vf/V0)/(steps*dt);
}

/* ---------- Stability via Linearization ---------- */
/* Lyapunov's indirect method: if the linearization at equilibrium
 * is asymptotically stable, then the nonlinear system is locally
 * asymptotically stable at that equilibrium. */

LDM_StabilityResult ldm_indirect_method(const LDM_System* s,
    const double* x_eq, int n) {
    if (!s||!x_eq||n<1) return LDM_INCONCLUSIVE;

    double* J=calloc(n*n,sizeof(double));
    ldm_jacobian(s,x_eq,J,n);
    bool hurwitz=ldm_linear_is_hurwitz(J,n);
    free(J);

    return hurwitz ? LDM_ASYM_STABLE : LDM_INCONCLUSIVE;
}

/* LaSalle's invariance principle check: for autonomous systems,
 * if V <= 0, trajectories converge to largest invariant set in {x: V_dot=0}. */
bool ldm_lasalle_check(const LDM_System* s, int n, double tol) {
    if (!s||n<1) return false;
    /* Simplified: check if V_dot < 0 at sampled non-equilibrium points.
     * This is a heuristic, not a proof. */
    double* x_save=calloc(n,sizeof(double));
    memcpy(x_save,s->x,n*sizeof(double));

    bool decreasing=true;
    double test_pts[6][3]={
        {0.5,0.0,0.0},{0.0,0.5,0.0},{0.5,0.5,0.0},
        {-0.5,0.3,0.0},{0.3,-0.5,0.0},{-0.5,-0.5,0.0}};

    for (int t=0;t<6&&decreasing;t++) {
        memcpy(s->x,test_pts[t],n*sizeof(double));
        double V=ldm_energy_based(s,n);
        double dV=ldm_energy_derivative(s,n);
        if (V>tol && dV> -1e-10) decreasing=false;
    }

    memcpy(s->x,x_save,n*sizeof(double));
    free(x_save);
    return decreasing;
}

/* Compare two Lyapunov function candidates */
double ldm_compare_candidates(const LDM_System* s, int n) {
    if (!s||n<1) return 0.0;
    double Vq=ldm_quadratic_form((double[]){1,0,0,1},s->x,n>2?2:n);
    double Vk=ldm_krasovskii_candidate(s,n);
    double Ve=ldm_energy_based(s,n);
    /* Return max of three candidates as "best" V estimate */
    double best=Vq;
    if (Vk>best) best=Vk;
    if (Ve>best) best=Ve;
    return best;
}

