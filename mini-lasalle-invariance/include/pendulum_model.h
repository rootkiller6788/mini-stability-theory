/* pendulum_model.h - Canonical nonlinear system models and Lyapunov functions.
 * L1: Canonical ODE definitions for stability theory benchmarks
 * L2: State-space representations of physical systems
 * L3: Parameterized ODE structure, equilibrium classification
 * L5: Lyapunov function construction per system type
 * L6: Damped pendulum, Van der Pol, Duffing, Lorenz, Lotka-Volterra, bistable
 * Ref: Khalil (2002) Nonlinear Systems, Strogatz (2015) Nonlinear Dynamics & Chaos
 */
#ifndef PENDULUM_MODEL_H
#define PENDULUM_MODEL_H
#include "gst_core.h"

typedef struct { double g_over_L; double b; } PendulumParams;
typedef struct { double mu; } VanDerPolParams;
typedef struct { double delta, alpha, beta, gamma, omega; } DuffingParams;
typedef struct { double sigma, rho, beta; } LorenzParams;
typedef struct { double alpha, beta, delta_param, gamma; } LotkaVolterraParams;

void damped_pendulum(const double* x, double* dxdt, int n, double t);
void damped_pendulum_ode(const double* x, double* dxdt, int n, double t, const PendulumParams* p);
double pendulum_energy_lyapunov(const double* x, int n);
double pendulum_energy_lyapunov_param(const double* x, int n, double gL);

void van_der_pol(const double* x, double* dxdt, int n, double t);
void van_der_pol_ode(const double* x, double* dxdt, int n, double t, const VanDerPolParams* p);
double vdp_candidate_lyapunov(const double* x, int n);

void duffing_oscillator(const double* x, double* dxdt, int n, double t);
void duffing_ode(const double* x, double* dxdt, int n, double t, const DuffingParams* p);
double duffing_energy(const double* x, int n, double alpha, double beta);
double duffing_lyapunov(const double* x, int n, double alpha, double beta, double delta);

void lorenz_system(const double* x, double* dxdt, int n, double t);
void lorenz_ode(const double* x, double* dxdt, int n, double t, const LorenzParams* p);
double lorenz_lyapunov_origin(const double* x, int n, const LorenzParams* p);
double lorenz_volume_contraction(const LorenzParams* p);

void lotka_volterra(const double* x, double* dxdt, int n, double t);
void lotka_volterra_ode(const double* x, double* dxdt, int n, double t, const LotkaVolterraParams* p);
double lv_constant_of_motion(const double* x, int n, const LotkaVolterraParams* p);
double lv_lyapunov(const double* x, int n, const LotkaVolterraParams* p);

void bistable_system(const double* x, double* dxdt, int n, double t);
void bistable_ode(const double* x, double* dxdt, int n, double t, double b);
double bistable_potential(double x);
double bistable_lyapunov_local(const double* x, int n, double equilibrium);

double quadratic_lyapunov(const double* x, int n);
double quadratic_lyapunov_weighted(const double* x, int n, const double* P_data);

void simple_pendulum(const double* x, double* dxdt, int n, double t);
double simple_pendulum_energy(const double* x, int n);

int find_equilibrium(ODEFunc f, double* x, int n, int max_iter, double tol);

#endif
