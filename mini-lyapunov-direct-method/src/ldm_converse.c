#include "ldm_converse.h"
#include "ldm_nonlinear.h"
#include "ldm_linear.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ==============================================================
 * ldm_converse.c — Converse Lyapunov Theorems
 *
 * Converse theorems: asymptotic stability => existence of a smooth
 * Lyapunov function. Construction methods:
 *   1. Massera integral: V(x) = ∫₀^∞ |φ(t,x)|² dt
 *   2. Kurzweil: V(x) = ∫₀^∞ w(|φ(t,x)|) dt
 *   3. Yoshizawa bounds: sup over time with weight function
 *
 * These are constructive: if we can simulate trajectories, we can
 * build a Lyapunov function numerically.
 * ============================================================== */

/* ---------- Internal helpers ---------- */

static double vec_norm_sq(const double* x, int n) {
    double s=0.0; for(int i=0;i<n;i++) s+=x[i]*x[i]; return s;
}

/* ---------- Massera Integral Construction ---------- */
/* V(x) = ∫₀^∞ |φ(t,x)|² dt
 * Where φ(t,x) is the solution starting from x at t=0.
 * For exponentially stable systems, this integral converges.
 * Implementation: sample trajectories and integrate |x(t)|².
 *
 * Algorithm:
 *   For each grid point (xi, yi):
 *     1. Set initial condition to (xi, yi)
 *     2. Simulate trajectory
 *     3. Integrate |x(t)|² * exp(-γt) dt (windowed for numerical stability)
 *     4. Estimate dV/dt = -|x|² (from the definition) */

LDM_ConverseResult ldm_converse_integral(const LDM_System* s, int n,
    double t_max, double dt, double* grid_params) {
    LDM_ConverseResult cr; memset(&cr,0,sizeof(cr));
    if (!s||n<1||dt<=0.0) return cr;

    /* Grid setup */
    cr.nx=40; cr.ny=40;
    if (grid_params) {
        cr.x_min=grid_params[0]; cr.x_max=grid_params[1];
        cr.y_min=grid_params[2]; cr.y_max=grid_params[3];
    } else {
        cr.x_min=-2.0; cr.x_max=2.0; cr.y_min=-2.0; cr.y_max=2.0;
    }

    int total_cells=cr.nx*cr.ny;
    cr.V_grid=calloc(total_cells,sizeof(double));
    cr.dV_grid=calloc(total_cells,sizeof(double));
    if (!cr.V_grid||!cr.dV_grid) {
        free(cr.V_grid); free(cr.dV_grid);
        cr.V_grid=cr.dV_grid=NULL;
        return cr;
    }

    int steps=(int)(t_max/dt);
    if (steps<1) steps=1;

    double* x_save=calloc(n,sizeof(double));
    if (x_save) memcpy(x_save,s->x,n*sizeof(double));

    LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
    if (!s2) { free(x_save); return cr; }

    /* For each grid cell, compute V */
    for (int iy=0;iy<cr.ny;iy++) {
        for (int ix=0;ix<cr.nx;ix++) {
            int idx=iy*cr.nx+ix;

            /* Grid point */
            double x0[4]={0};
            if (n>=1) x0[0]=cr.x_min+(cr.x_max-cr.x_min)*ix/(cr.nx-1.0);
            if (n>=2) x0[1]=cr.y_min+(cr.y_max-cr.y_min)*iy/(cr.ny-1.0);

            memcpy(s2->x,x0,n*sizeof(double));
            s2->t=0.0;

            /* Integrate |φ(t,x0)|² along trajectory */
            double V=0.0;
            double prev_norm_sq=vec_norm_sq(s2->x,n);

            for (int k=1;k<steps;k++) {
                ldm_system_step(s2,dt);
                double norm_sq=vec_norm_sq(s2->x,n);

                /* Trapezoidal integration: ∫|x(t)|² dt */
                V+=0.5*(prev_norm_sq+norm_sq)*dt;

                prev_norm_sq=norm_sq;

                /* Early termination if converged */
                if (norm_sq<1e-20) break;
            }

            cr.V_grid[idx]=V;

            /* dV/dt = -|x|² (from the definition of Massera's V) */
            cr.dV_grid[idx]=-vec_norm_sq(x0,n);
        }
    }

    cr.is_valid=true;

    /* Restore state */
    if (x_save) { memcpy(s->x,x_save,n*sizeof(double)); free(x_save); }
    ldm_system_free(s2);
    return cr;
}

/* ---------- Sampling-Based Converse ---------- */
/* Randomly sample initial conditions, simulate trajectories,
 * and compute V(x) = ∫₀^∞ |φ(t,x)|² dt at each sample point. */

LDM_ConverseResult ldm_converse_sampling(const LDM_System* s, int n,
    int n_samples, double range) {
    LDM_ConverseResult cr; memset(&cr,0,sizeof(cr));
    if (!s||n<1||n_samples<1) return cr;

    int grid_size=(int)sqrt((double)n_samples)+1;
    cr.nx=grid_size; cr.ny=grid_size;
    cr.x_min=-range; cr.x_max=range;
    cr.y_min=-range; cr.y_max=range;

    int total=cr.nx*cr.ny;
    cr.V_grid=calloc(total,sizeof(double));
    cr.dV_grid=calloc(total,sizeof(double));
    if (!cr.V_grid||!cr.dV_grid) {
        free(cr.V_grid); free(cr.dV_grid);
        cr.V_grid=cr.dV_grid=NULL;
        return cr;
    }

    double* x_save=calloc(n,sizeof(double));
    if (x_save) memcpy(x_save,s->x,n*sizeof(double));

    LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
    if (!s2) { free(x_save); return cr; }

    double dt=0.01;
    int max_steps=2000;

    /* Uniform grid of initial conditions */
    for (int iy=0;iy<cr.ny;iy++) {
        for (int ix=0;ix<cr.nx;ix++) {
            int idx=iy*cr.nx+ix;
            double x0[4]={0};
            if (n>=1) x0[0]=-range+2.0*range*ix/(cr.nx-1.0);
            if (n>=2) x0[1]=-range+2.0*range*iy/(cr.ny-1.0);

            memcpy(s2->x,x0,n*sizeof(double));
            s2->t=0.0;

            double V=0.0, prev_norm_sq=vec_norm_sq(s2->x,n);
            for (int k=1;k<max_steps;k++) {
                ldm_system_step(s2,dt);
                double norm_sq=vec_norm_sq(s2->x,n);
                V+=0.5*(prev_norm_sq+norm_sq)*dt;
                prev_norm_sq=norm_sq;
                if (norm_sq<1e-20) break;
            }
            cr.V_grid[idx]=V;
            cr.dV_grid[idx]=-vec_norm_sq(x0,n);
        }
    }

    cr.is_valid=true;
    if (x_save) { memcpy(s->x,x_save,n*sizeof(double)); free(x_save); }
    ldm_system_free(s2);
    return cr;
}

/* ---------- Lifecycle ---------- */

void ldm_converse_free(LDM_ConverseResult* cr) {
    if (cr) { free(cr->V_grid); free(cr->dV_grid); }
}

void ldm_converse_print(const LDM_ConverseResult* cr) {
    if (!cr) return;
    printf("Converse Lyapunov: %dx%d grid, [%.2f,%.2f]x[%.2f,%.2f] %s\n",
        cr->nx, cr->ny, cr->x_min, cr->x_max, cr->y_min, cr->y_max,
        cr->is_valid?"valid":"invalid");
    if (cr->is_valid&&cr->V_grid) {
        printf("  V(0,0)=%.4f  V(1,0)=%.4f  dV(0,0)=%.4f\n",
            cr->V_grid[0], cr->V_grid[cr->nx/2],
            cr->dV_grid[0]);
    }
}

/* ---------- Accuracy Assessment ---------- */

double ldm_converse_accuracy(const LDM_ConverseResult* cr) {
    if (!cr||!cr->is_valid) return 0.0;

    /* Check: V(x) should be positive away from origin.
     * Accuracy = fraction of grid points with V>0 when |x|>ε */
    int correct=0, away_from_origin=0;

    for (int iy=0;iy<cr->ny;iy++) {
        double y=cr->y_min+(cr->y_max-cr->y_min)*iy/(cr->ny-1.0);
        for (int ix=0;ix<cr->nx;ix++) {
            double x=cr->x_min+(cr->x_max-cr->x_min)*ix/(cr->nx-1.0);
            double dist_sq=x*x+y*y;
            if (dist_sq>0.01) {
                away_from_origin++;
                int idx=iy*cr->nx+ix;
                if (cr->V_grid[idx]>0.0) correct++;
            }
        }
    }

    return away_from_origin>0 ? (double)correct/away_from_origin : 0.0;
}

/* ---------- Massera Condition ---------- */
/* Massera (1949): Asymptotic stability implies the existence of a
 * smooth Lyapunov function. The necessary condition is that the system
 * is asymptotically stable at the origin.
 *
 * This function checks the necessary condition: does the linearization
 * have all eigenvalues with strictly negative real parts? */

bool ldm_massera_condition(const LDM_System* s, int n) {
    if (!s||n<1) return false;

    /* Check linearization at origin */
    double x0[8]={0};
    double* J=calloc(n*n,sizeof(double));
    ldm_jacobian(s,x0,J,n);
    bool satisfied=ldm_linear_is_hurwitz(J,n);
    free(J);
    return satisfied;
}

/* ---------- Yoshizawa Bound ---------- */
/* Yoshizawa (1966): For asymptotically stable systems,
 * ||φ(t,x)|| ≤ β(||x||,0) * σ(t) for class KL function β and
 * class L function σ.
 *
 * This function estimates the maximum norm along a trajectory,
 * giving an upper bound for the Yoshizawa construction. */

double ldm_yoshizawa_bound(const LDM_System* s, int n, double t_max) {
    if (!s||n<1) return 0.0;

    double* x_save=calloc(n,sizeof(double));
    if (x_save) memcpy(x_save,s->x,n*sizeof(double));

    LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
    if (!s2) { free(x_save); return 0.0; }

    /* Start from unit sphere to get normalized bound */
    double max_expansion=0.0;
    int n_test=8;
    for (int t=0;t<n_test;t++) {
        double angle=2.0*M_PI*t/n_test;
        if (n>=2) { s2->x[0]=cos(angle); s2->x[1]=sin(angle); }
        else { s2->x[0]=(t<n_test/2)?1.0:-1.0; }

        double dt=0.01;
        int steps=(int)(t_max/dt);
        double max_norm=0.0;

        for (int k=0;k<steps;k++) {
            ldm_system_step(s2,dt);
            double norm=ldm_vector_norm(s2->x,n);
            if (norm>max_norm) max_norm=norm;
        }
        if (max_norm>max_expansion) max_expansion=max_norm;
    }

    ldm_system_free(s2);
    if (x_save) { memcpy(s->x,x_save,n*sizeof(double)); free(x_save); }

    return max_expansion;
}

/* ---------- Kurzweil's Construction ---------- */
/* V(x) = ∫₀^∞ w(||φ(t,x)||) dt where w is a class-K function.
 * This generalizes Massera's integral (w(s)=s²).
 * Implementation with w(s) = s² (default) or w(s) = s. */

LDM_ConverseResult ldm_converse_kurzweil(const LDM_System* s, int n,
    double t_max, double dt, double p) {
    /* p determines the weight: w(s) = s^p */
    LDM_ConverseResult cr; memset(&cr,0,sizeof(cr));
    if (!s||n<1||dt<=0.0) return cr;

    cr.nx=30; cr.ny=30;
    cr.x_min=-2.0; cr.x_max=2.0;
    cr.y_min=-2.0; cr.y_max=2.0;

    int total=cr.nx*cr.ny;
    cr.V_grid=calloc(total,sizeof(double));
    cr.dV_grid=calloc(total,sizeof(double));
    if (!cr.V_grid||!cr.dV_grid) {
        free(cr.V_grid); free(cr.dV_grid);
        cr.V_grid=cr.dV_grid=NULL;
        return cr;
    }

    double* x_save=calloc(n,sizeof(double));
    if (x_save) memcpy(x_save,s->x,n*sizeof(double));

    LDM_System* s2=ldm_system_create(n,s->dynamics,s->context);
    if (!s2) { free(x_save); return cr; }

    int steps=(int)(t_max/dt);

    for (int iy=0;iy<cr.ny;iy++) {
        for (int ix=0;ix<cr.nx;ix++) {
            int idx=iy*cr.nx+ix;
            double x0[4]={0};
            if (n>=1) x0[0]=cr.x_min+(cr.x_max-cr.x_min)*ix/(cr.nx-1.0);
            if (n>=2) x0[1]=cr.y_min+(cr.y_max-cr.y_min)*iy/(cr.ny-1.0);

            memcpy(s2->x,x0,n*sizeof(double));
            s2->t=0.0;

            double V=0.0;
            double prev_norm=sqrt(vec_norm_sq(s2->x,n));
            double prev_w=pow(prev_norm,p);

            for (int k=1;k<steps;k++) {
                ldm_system_step(s2,dt);
                double norm=sqrt(vec_norm_sq(s2->x,n));
                double w=pow(norm,p);
                V+=0.5*(prev_w+w)*dt;
                prev_w=w;
                if (norm<1e-10) break;
            }
            cr.V_grid[idx]=V;

            /* dV/dt = -w(||x||) = -||x||^p */
            cr.dV_grid[idx]=-pow(sqrt(vec_norm_sq(x0,n)),p);
        }
    }

    cr.is_valid=true;
    if (x_save) { memcpy(s->x,x_save,n*sizeof(double)); free(x_save); }
    ldm_system_free(s2);
    return cr;
}

/* ---------- Convergence Rate from Converse Function ---------- */
/* From V(x) and dV/dt = -|x|², estimate the exponential decay rate:
 * α ≈ -dV/dt / V = |x|² / V(x) */

double ldm_converse_decay_rate(const LDM_ConverseResult* cr, int ix, int iy) {
    if (!cr||!cr->is_valid||!cr->V_grid) return 0.0;
    if (ix<0||ix>=cr->nx||iy<0||iy>=cr->ny) return 0.0;
    int idx=iy*cr->nx+ix;
    if (cr->V_grid[idx]<1e-14) return 0.0;
    return -cr->dV_grid[idx]/cr->V_grid[idx];
}

/* Check if the constructed converse V is positive definite */
bool ldm_converse_is_positive_definite(const LDM_ConverseResult* cr) {
    if (!cr||!cr->is_valid||!cr->V_grid) return false;
    /* Check V(0,0)≈0 and V(x)>0 for x≠0 */
    int cx=cr->nx/2, cy=cr->ny/2;
    double v_origin=cr->V_grid[cy*cr->nx+cx];
    if (v_origin>0.1) return false; /* V(0) should be near 0 */

    /* Sample away from origin */
    for (int i=0;i<4;i++) {
        int ix=cx+(i+1)*3;
        if (ix<cr->nx) {
            double v=cr->V_grid[cy*cr->nx+ix];
            if (v<=0.0) return false;
        }
    }
    return true;
}

/* ---------- Analysis utilities for Converse Results ---------- */

/* Find maximum value in the V grid (useful for normalization) */
double ldm_converse_max_value(const LDM_ConverseResult* cr) {
    if (!cr||!cr->is_valid||!cr->V_grid) return 0.0;
    int N=cr->nx*cr->ny;
    double maxv=cr->V_grid[0];
    for (int i=1;i<N;i++) if (cr->V_grid[i]>maxv) maxv=cr->V_grid[i];
    return maxv;
}

/* Find minimum value in the V grid (should be near 0 at origin) */
double ldm_converse_min_value(const LDM_ConverseResult* cr) {
    if (!cr||!cr->is_valid||!cr->V_grid) return 0.0;
    int N=cr->nx*cr->ny;
    double minv=cr->V_grid[0];
    for (int i=1;i<N;i++) if (cr->V_grid[i]<minv) minv=cr->V_grid[i];
    return minv;
}

/* Average dV/dt over the grid (should be negative for stable systems) */
double ldm_converse_average_dV(const LDM_ConverseResult* cr) {
    if (!cr||!cr->is_valid||!cr->dV_grid) return 0.0;
    int N=cr->nx*cr->ny;
    double sum=0.0;
    for (int i=0;i<N;i++) sum+=cr->dV_grid[i];
    return sum/N;
}

/* Check monotonicity: V should increase with distance from origin */
double ldm_converse_monotonicity(const LDM_ConverseResult* cr) {
    if (!cr||!cr->is_valid||!cr->V_grid||cr->nx<3||cr->ny<3) return -1.0;
    /* Check that V along a ray increases */
    int cx=cr->nx/2, cy=cr->ny/2;
    int violations=0, total=0;
    for (int i=1;i<4;i++) {
        int ix=cx+i;
        if (ix<cr->nx) {
            total++;
            if (cr->V_grid[cy*cr->nx+ix]<=cr->V_grid[cy*cr->nx+cx])
                violations++;
        }
    }
    return total>0 ? 1.0-(double)violations/total : 0.0;
}

/* Smoothness check: estimate Lipschitz constant of V over the grid */
double ldm_converse_lipschitz_estimate(const LDM_ConverseResult* cr) {
    if (!cr||!cr->is_valid||!cr->V_grid||cr->nx<2||cr->ny<2) return 1e10;
    double L=0.0;
    double dx=(cr->x_max-cr->x_min)/(cr->nx-1.0);
    double dy=(cr->y_max-cr->y_min)/(cr->ny-1.0);
    for (int iy=0;iy<cr->ny;iy++) {
        for (int ix=1;ix<cr->nx;ix++) {
            double diff=fabs(cr->V_grid[iy*cr->nx+ix]-cr->V_grid[iy*cr->nx+ix-1])/dx;
            if (diff>L) L=diff;
        }
    }
    for (int ix=0;ix<cr->nx;ix++) {
        for (int iy=1;iy<cr->ny;iy++) {
            double diff=fabs(cr->V_grid[iy*cr->nx+ix]-cr->V_grid[(iy-1)*cr->nx+ix])/dy;
            if (diff>L) L=diff;
        }
    }
    return L;
}

/* Norm of V evaluated at unit circle (discrete approximation) */
double ldm_converse_unit_circle_norm(const LDM_ConverseResult* cr) {
    if (!cr||!cr->is_valid||!cr->V_grid) return 0.0;
    int n_angles=16;
    double sum=0.0;
    for (int k=0;k<n_angles;k++) {
        double angle=2.0*M_PI*k/n_angles;
        double x=cos(angle), y=sin(angle);
        double norm=sqrt(x*x+y*y);
        if (norm<0.01) continue;
        /* Map to grid index */
        int ix=(int)((x-cr->x_min)/(cr->x_max-cr->x_min)*(cr->nx-1));
        int iy=(int)((y-cr->y_min)/(cr->y_max-cr->y_min)*(cr->ny-1));
        if (ix<0) ix=0;
        if (ix>=cr->nx) ix=cr->nx-1;
        if (iy<0) iy=0;
        if (iy>=cr->ny) iy=cr->ny-1;
        sum+=cr->V_grid[iy*cr->nx+ix];
    }
    return sum/n_angles;
}
