/* pendulum_model.c - Canonical nonlinear system models and Lyapunov functions.
 * L1-L6: ODE definitions, state-space models, Lyapunov function construction.
 * Each function implements an independent knowledge point:
 *   damped_pendulum: energy-based Lyapunov, LaSalle GAS proof (L4/L6)
 *   van_der_pol: limit cycle existence via Poincare-Bendixson (L6)
 *   duffing: potential-based Lyapunov, multiple equilibria (L6)
 *   lorenz: volume contraction, global stability for rho<1 (L6/L8)
 *   lotka_volterra: first integral, neutral stability (L7)
 *   bistable: double-well potential, local stability (L6)
 *   quadratic_lyapunov: generic construction for linearized systems (L2/L5)
 *   find_equilibrium: numerical root-finding for ODE equilibria (L5)
 *   simple_pendulum: Hamiltonian system, energy conservation (L6)
 * Ref: Khalil (2002) Nonlinear Systems, Strogatz (2015) Nonlinear Dynamics
 */
#include "pendulum_model.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void damped_pendulum(const double* x, double* dxdt, int n, double t)
{ (void)t; (void)n; dxdt[0] = x[1]; dxdt[1] = -0.5 * x[1] - 9.81 * sin(x[0]); }

void damped_pendulum_ode(const double* x, double* dxdt, int n, double t,
                         const PendulumParams* p)
{
    (void)t; (void)n;
    if (!p) { dxdt[0] = 0; dxdt[1] = 0; return; }
    dxdt[0] = x[1];
    dxdt[1] = -p->b * x[1] - p->g_over_L * sin(x[0]);
}

double pendulum_energy_lyapunov(const double* x, int n)
{ (void)n; return 0.5*x[1]*x[1] + 9.81*(1.0 - cos(x[0])); }

double pendulum_energy_lyapunov_param(const double* x, int n, double gL)
{ (void)n; return 0.5*x[1]*x[1] + gL*(1.0 - cos(x[0])); }

void van_der_pol(const double* x, double* dxdt, int n, double t)
{ (void)t; (void)n; dxdt[0]=x[1]; dxdt[1]=-1.0*(x[0]*x[0]-1.0)*x[1] - x[0]; }

void van_der_pol_ode(const double* x, double* dxdt, int n, double t,
                     const VanDerPolParams* p)
{
    (void)t; (void)n;
    if (!p) { dxdt[0]=0; dxdt[1]=0; return; }
    dxdt[0]=x[1];
    dxdt[1]=-p->mu*(x[0]*x[0]-1.0)*x[1] - x[0];
}

double vdp_candidate_lyapunov(const double* x, int n)
{ (void)n; return 0.5*(x[0]*x[0] + x[1]*x[1]); }

void duffing_oscillator(const double* x, double* dxdt, int n, double t)
{ (void)n; (void)t; dxdt[0]=x[1]; dxdt[1]=-0.3*x[1]+x[0]-x[0]*x[0]*x[0]; }

void duffing_ode(const double* x, double* dxdt, int n, double t,
                 const DuffingParams* p)
{
    (void)n;
    if (!p) { dxdt[0]=0; dxdt[1]=0; return; }
    double frc = p->gamma * cos(p->omega * t);
    dxdt[0] = x[1];
    dxdt[1] = -p->delta*x[1] - p->alpha*x[0] - p->beta*x[0]*x[0]*x[0] + frc;
}

double duffing_energy(const double* x, int n, double alpha, double beta)
{
    (void)n;
    return 0.5*x[1]*x[1] + 0.5*alpha*x[0]*x[0] + 0.25*beta*x[0]*x[0]*x[0]*x[0];
}

double duffing_lyapunov(const double* x, int n, double alpha, double beta, double delta)
{
    (void)delta;
    return duffing_energy(x, n, alpha, beta);
}

void lorenz_system(const double* x, double* dxdt, int n, double t)
{
    (void)t; (void)n;
    dxdt[0] = 10.0 * (x[1] - x[0]);
    dxdt[1] = x[0] * (28.0 - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - (8.0/3.0) * x[2];
}

void lorenz_ode(const double* x, double* dxdt, int n, double t,
                const LorenzParams* p)
{
    (void)t; (void)n;
    if (!p) { dxdt[0]=0; dxdt[1]=0; dxdt[2]=0; return; }
    dxdt[0] = p->sigma * (x[1] - x[0]);
    dxdt[1] = x[0] * (p->rho - x[2]) - x[1];
    dxdt[2] = x[0] * x[1] - p->beta * x[2];
}

double lorenz_lyapunov_origin(const double* x, int n, const LorenzParams* p)
{
    (void)n;
    if (!p) return 0.0;
    return 0.5*(x[0]*x[0] + p->sigma*x[1]*x[1] + p->sigma*x[2]*x[2]);
}

double lorenz_volume_contraction(const LorenzParams* p)
{ if (!p) return 0.0; return -(p->sigma + 1.0 + p->beta); }

void lotka_volterra(const double* x, double* dxdt, int n, double t)
{
    (void)t; (void)n;
    dxdt[0] = 1.0*x[0] - 0.1*x[0]*x[1];
    dxdt[1] = 0.075*x[0]*x[1] - 1.5*x[1];
}

void lotka_volterra_ode(const double* x, double* dxdt, int n, double t,
                        const LotkaVolterraParams* p)
{
    (void)t; (void)n;
    if (!p) { dxdt[0]=0; dxdt[1]=0; return; }
    dxdt[0] = p->alpha*x[0] - p->beta*x[0]*x[1];
    dxdt[1] = p->delta_param*x[0]*x[1] - p->gamma*x[1];
}

double lv_constant_of_motion(const double* x, int n, const LotkaVolterraParams* p)
{
    (void)n;
    if (!p || x[0]<=0 || x[1]<=0) return 0.0;
    return p->delta_param*x[0] - p->gamma*log(x[0])
         + p->beta*x[1] - p->alpha*log(x[1]);
}

double lv_lyapunov(const double* x, int n, const LotkaVolterraParams* p)
{
    (void)n;
    if (!p || x[0]<=0 || x[1]<=0) return 0.0;
    double xs = p->gamma / p->delta_param;
    double ys = p->alpha / p->beta;
    return p->delta_param*(x[0]-xs*log(x[0]/xs))
         + p->beta*(x[1]-ys*log(x[1]/ys));
}

void bistable_system(const double* x, double* dxdt, int n, double t)
{ (void)t; (void)n; dxdt[0]=x[1]; dxdt[1]=-0.3*x[1]+x[0]-x[0]*x[0]*x[0]; }

void bistable_ode(const double* x, double* dxdt, int n, double t, double b)
{ (void)t; (void)n; dxdt[0]=x[1]; dxdt[1]=-b*x[1]+x[0]-x[0]*x[0]*x[0]; }

double bistable_potential(double x)
{ return -0.5*x*x + 0.25*x*x*x*x; }

double bistable_lyapunov_local(const double* x, int n, double eq)
{
    (void)n;
    double Ux = -0.5*x[0]*x[0] + 0.25*x[0]*x[0]*x[0]*x[0];
    double Ueq = -0.5*eq*eq + 0.25*eq*eq*eq*eq;
    return 0.5*x[1]*x[1] + Ux - Ueq;
}

double quadratic_lyapunov(const double* x, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += x[i]*x[i];
    return 0.5*s;
}

double quadratic_lyapunov_weighted(const double* x, int n, const double* P)
{
    if (!x || !P || n<=0) return 0.0;
    double s = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            s += x[i] * P[i*n + j] * x[j];
    return s;
}

void simple_pendulum(const double* x, double* dxdt, int n, double t)
{ (void)t; (void)n; dxdt[0]=x[1]; dxdt[1]=-sin(x[0]); }

double simple_pendulum_energy(const double* x, int n)
{ (void)n; return 0.5*x[1]*x[1] + (1.0-cos(x[0])); }

int find_equilibrium(ODEFunc f, double* x, int n, int max_iter, double tol)
{
    if (!f || !x || n<=0 || max_iter<=0) return -1;
    double* fx = (double*)calloc((size_t)n, sizeof(double));
    double* xnew = (double*)calloc((size_t)n, sizeof(double));
    if (!fx || !xnew) { free(fx); free(xnew); return -1; }

    for (int iter = 0; iter < max_iter; iter++) {
        f(x, fx, n, 0.0);
        double norm = 0.0;
        for (int d = 0; d < n; d++) norm += fx[d]*fx[d];
        norm = sqrt(norm);
        if (norm < tol) { free(fx); free(xnew); return iter; }

        double eps_fd = 1e-6;
        for (int d = 0; d < n; d++) {
            double* xp = (double*)calloc((size_t)n, sizeof(double));
            if (!xp) { free(fx); free(xnew); return -1; }
            memcpy(xp, x, (size_t)n*sizeof(double));
            double h = eps_fd * (1.0 + fabs(x[d]));
            xp[d] += h;
            double* fp = (double*)calloc((size_t)n, sizeof(double));
            if (!fp) { free(xp); free(fx); free(xnew); return -1; }
            f(xp, fp, n, 0.0);
            double Jtf = 0.0;
            for (int i = 0; i < n; i++)
                Jtf += (fp[i] - fx[i]) / h * fx[i];
            xnew[d] = x[d] - 0.01 * Jtf;
            free(xp); free(fp);
        }
        memcpy(x, xnew, (size_t)n*sizeof(double));
    }
    free(fx); free(xnew);
    return -1;
}