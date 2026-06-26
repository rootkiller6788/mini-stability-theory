#include "lasalle_invariance.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

LasalleResult* lasalle_create(void) {
    LasalleResult* lr = calloc(1, sizeof(LasalleResult));
    if (lr) lr->convergence_tolerance = 1e-6;
    return lr;
}

void lasalle_free(LasalleResult* lr) {
    if (!lr) return;
    free(lr->omega_points);
    free(lr->zero_derivative_set);
    free(lr->largest_invariant);
    free(lr);
}

int lasalle_find_zero_derivative_set(LasalleResult* lr, LyapunovFunc V,
    ODEFunc f, const double* x0, int n, double dt,
    double duration, double eps)
{
    if (!lr || !V || !f || !x0) return -1;
    double* x = malloc((size_t)n * sizeof(double));
    memcpy(x, x0, (size_t)n * sizeof(double));
    int steps = (int)(duration / dt);
    lr->cap_zero = steps;
    lr->zero_derivative_set = malloc((size_t)(steps * n) * sizeof(double));
    if (!lr->zero_derivative_set) { free(x); return -1; }
    lr->n_zero = 0;
    for (int i = 0; i < steps; i++) {
        double vdot = lyapunov_derivative(f, V, x, n, eps);
        if (fabs(vdot) < eps && lr->n_zero < lr->cap_zero) {
            memcpy(lr->zero_derivative_set + lr->n_zero * n,
                   x, (size_t)n * sizeof(double));
            lr->n_zero++;
        }
        rk4_step(f, x, n, 0.0, dt);
    }
    free(x);
    return lr->n_zero;
}

int lasalle_find_largest_invariant(LasalleResult* lr, ODEFunc f,
    const double* x0, int n, int max_iter, double dt, double eps)
{
    if (!lr || lr->n_zero == 0) return -1;
    lr->cap_invariant = lr->n_zero + 1;
    lr->largest_invariant = malloc((size_t)(lr->cap_invariant * n) * sizeof(double));
    if (!lr->largest_invariant) return -1;
    lr->n_invariant = 0;

    /* 1) check each zero-Vdot point for forward invariance */
    for (int i = 0; i < lr->n_zero; i++) {
        double* x = malloc((size_t)n * sizeof(double));
        memcpy(x, lr->zero_derivative_set + i * n, (size_t)n * sizeof(double));
        int is_inv = 1;
        for (int j = 0; j < max_iter && is_inv; j++) {
            rk4_step(f, x, n, 0.0, dt);
            int found = 0;
            for (int k = 0; k < lr->n_zero; k++) {
                double dist = 0.0;
                for (int d = 0; d < n; d++) {
                    double diff = x[d] - lr->zero_derivative_set[k * n + d];
                    dist += diff * diff;
                }
                if (sqrt(dist) < eps) { found = 1; break; }
            }
            if (!found) is_inv = 0;
        }
        if (is_inv) {
            memcpy(lr->largest_invariant + lr->n_invariant * n,
                   lr->zero_derivative_set + i * n, (size_t)n * sizeof(double));
            lr->n_invariant++;
        }
        free(x);
    }

    /* 2) fallback: simulate forward from last zero-Vdot point toward equilibrium */
    if (lr->n_invariant == 0 && lr->n_zero > 0) {
        double* x = malloc((size_t)n * sizeof(double));
        memcpy(x, lr->zero_derivative_set + (lr->n_zero - 1) * n,
               (size_t)n * sizeof(double));

        /* use a practical tolerance for convergence detection */
        double conv_tol = eps > 1e-8 ? eps * 2000.0 : 1e-4;
        int max_sim_steps = max_iter * 20;
        for (int j = 0; j < max_sim_steps; j++) {
            rk4_step(f, x, n, 0.0, dt);
            if (j % 50 == 49) {
                double fx[16];
                f(x, fx, n, 0.0);
                double vnorm = 0.0;
                for (int d = 0; d < n; d++) vnorm += fx[d] * fx[d];
                vnorm = sqrt(vnorm);
                if (vnorm < conv_tol) {
                    memcpy(lr->largest_invariant + lr->n_invariant * n,
                           x, (size_t)n * sizeof(double));
                    lr->n_invariant++;
                    break;
                }
            }
        }
        free(x);
    }

    lr->is_converged = (lr->n_invariant > 0);
    lr->convergence_time = (double)max_iter * dt;
    return lr->n_invariant;
}

bool lasalle_verify_conditions(LyapunovFunc V, ODEFunc f,
    const double* x_sample, int n, int n_samples,
    double radius, double eps)
{
    (void)radius;
    if (!V || !f || !x_sample) return false;
    for (int i = 0; i < n_samples; i++) {
        double vdot = lyapunov_derivative(f, V, x_sample + i * n, n, eps);
        if (vdot > eps) return false;
    }
    return true;
}

double lasalle_convergence_time_estimate(const LasalleResult* lr) {
    return lr ? lr->convergence_time : 0.0;
}

bool lasalle_is_asymptotically_stable(const LasalleResult* lr) {
    return lr && lr->is_converged && lr->n_invariant > 0;
}

void lasalle_print(const LasalleResult* lr) {
    if (!lr) { printf("LasalleResult: NULL\n"); return; }
    printf("=== LaSalle Invariance ===\n");
    printf("  Zero Vdot points: %d\n", lr->n_zero);
    printf("  Invariant set size: %d\n", lr->n_invariant);
    printf("  Asympt. stable: %s\n",
           lasalle_is_asymptotically_stable(lr) ? "YES" : "NO");
}

InvarianceAnalysis* invar_analyze(ODEFunc f, LyapunovFunc V,
    const double* x0, int n, double dt, double duration)
{
    if (!f || !V || !x0) return NULL;
    InvarianceAnalysis* ia = calloc(1, sizeof(InvarianceAnalysis));
    if (!ia) return NULL;
    int steps = (int)(duration / dt);
    ia->n_steps = steps;
    ia->n_dims = n;
    ia->x_trajectory = malloc((size_t)(steps * n) * sizeof(double));
    ia->V_values = malloc((size_t)steps * sizeof(double));
    ia->Vdot_values = malloc((size_t)steps * sizeof(double));
    if (!ia->x_trajectory || !ia->V_values || !ia->Vdot_values) {
        invar_free(ia);
        return NULL;
    }
    double* x = malloc((size_t)n * sizeof(double));
    memcpy(x, x0, (size_t)n * sizeof(double));
    for (int i = 0; i < steps; i++) {
        memcpy(ia->x_trajectory + i * n, x, (size_t)n * sizeof(double));
        ia->V_values[i] = V(x, n);
        ia->Vdot_values[i] = lyapunov_derivative(f, V, x, n, 1e-6);
        rk4_step(f, x, n, 0.0, dt);
    }
    ia->final_V = ia->V_values[steps - 1];
    ia->final_Vdot = ia->Vdot_values[steps - 1];
    ia->approaches_invariant = (ia->final_Vdot > -1e-6);
    free(x);
    return ia;
}

void invar_free(InvarianceAnalysis* ia) {
    if (!ia) return;
    free(ia->x_trajectory);
    free(ia->V_values);
    free(ia->Vdot_values);
    free(ia);
}

double invar_compute_distance(const InvarianceAnalysis* ia,
    const double* inv_set, int n_inv, int n_dims)
{
    if (!ia || !inv_set || n_inv == 0) return 1e10;
    double md = 1e10;
    for (int i = 0; i < ia->n_steps; i++)
        for (int j = 0; j < n_inv; j++) {
            double d = 0.0;
            for (int k = 0; k < n_dims; k++) {
                double df = ia->x_trajectory[i * n_dims + k]
                          - inv_set[j * n_dims + k];
                d += df * df;
            }
            d = sqrt(d);
            if (d < md) md = d;
        }
    ((InvarianceAnalysis*)ia)->distance_to_invariant = md;
    return md;
}

bool invar_within_tolerance(const InvarianceAnalysis* ia, double tol) {
    return ia && ((InvarianceAnalysis*)ia)->distance_to_invariant < tol;
}

void invar_print(const InvarianceAnalysis* ia) {
    if (!ia) return;
    printf("=== Invariance Analysis ===\n");
    printf("  Steps: %d  Final V: %.6f  Vdot: %.6f\n",
           ia->n_steps, ia->final_V, ia->final_Vdot);
    printf("  Approaches invariant: %s\n",
           ia->approaches_invariant ? "YES" : "NO");
}

/* ================================================================
 * L5: Lasalle trajectory energy profile
 * Compute V(t) and Vdot(t) along a trajectory and check if
 * V monotonically decreases (necessary condition for LaSalle).
 * Returns the number of violations of monotone V decrease.
 * ================================================================ */
int lasalle_trajectory_energy_profile(ODEFunc f, LyapunovFunc V,
    const double* x0, int n, double dt, double duration,
    double* V_history, double* Vdot_history, int max_steps)
{
    if (!f || !V || !x0 || n <= 0 || !V_history || !Vdot_history || max_steps <= 0)
        return -1;
    int steps = (int)(duration / dt);
    if (steps > max_steps) steps = max_steps;
    double* x = malloc((size_t)n * sizeof(double));
    if (!x) return -1;
    memcpy(x, x0, (size_t)n * sizeof(double));

    int violations = 0;
    double prev_V = V(x, n);
    V_history[0] = prev_V;
    Vdot_history[0] = lyapunov_derivative(f, V, x, n, 1e-6);

    for (int i = 1; i < steps; i++) {
        rk4_step(f, x, n, 0.0, dt);
        double cur_V = V(x, n);
        V_history[i] = cur_V;
        Vdot_history[i] = lyapunov_derivative(f, V, x, n, 1e-6);
        if (cur_V > prev_V + 1e-6) violations++;
        prev_V = cur_V;
    }
    free(x);
    return violations;
}

/* ================================================================
 * L8: Lasalle-Yoshizawa extension for non-autonomous systems.
 * For x_dot = f(x,t) with V(x,t) s.t. Vdot(x,t) <= -W(x) <= 0,
 * all trajectories are bounded and converge to {x: W(x)=0}.
 * This function checks the condition Vdot(x,t) + W3(x) <= 0
 * where W3 is a user-supplied positive-definite function.
 * Ref: Yoshizawa (1966) Stability Theory by Liapunov's Second Method
 * ================================================================ */
int lasalle_yoshizawa_check(ODEFunc f, LyapunovFunc V,
    LyapunovFunc W3, const double* samples, int n_samples,
    int n_dims, double eps)
{
    if (!f || !V || !W3 || !samples || n_samples <= 0 || n_dims <= 0)
        return -1;
    int violations = 0;
    for (int i = 0; i < n_samples; i++) {
        double vdot = lyapunov_derivative(f, V,
            samples + i * n_dims, n_dims, eps);
        double w_val = W3(samples + i * n_dims, n_dims);
        if (vdot + w_val > eps) violations++;
    }
    return violations;
}

/* ================================================================
 * L8: LaSalle for gradient systems.
 * For x_dot = -grad V(x) (gradient system), V is a strict Lyapunov
 * function with Vdot = -||grad V||^2 <= 0. The invariant set is
 * the set of critical points of V. This function verifies if a
 * given system is a gradient system by checking if the Jacobian is
 * symmetric (f_i = -dV/dx_i => df_i/dx_j = df_j/dx_i).
 * Ref: Hirsch & Smale (1974) Differential Equations, Ch.9
 * ================================================================ */
bool lasalle_check_gradient_system(ODEFunc f, int n,
    const double* x_test, double eps)
{
    if (!f || n <= 1 || !x_test) return false;
    double* fx = malloc((size_t)n * sizeof(double));
    if (!fx) return false;
    f(x_test, fx, n, 0.0);

    /* Check symmetry of Jacobian via finite differences */
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double xi[16]; memcpy(xi, x_test, (size_t)n*sizeof(double));
            double xj[16]; memcpy(xj, x_test, (size_t)n*sizeof(double));
            double h = eps * (1.0 + fabs(x_test[i]));
            xi[i] += h;
            double fi[16]; f(xi, fi, n, 0.0);
            double dfj_dxi = (fi[j] - fx[j]) / h;

            h = eps * (1.0 + fabs(x_test[j]));
            xj[j] += h;
            double fj[16]; f(xj, fj, n, 0.0);
            double dfi_dxj = (fj[i] - fx[i]) / h;

            if (fabs(dfj_dxi - dfi_dxj) > eps * 10.0) {
                free(fx); return false;
            }
        }
    }
    free(fx);
    return true;
}

