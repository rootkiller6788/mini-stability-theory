#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fts_core.h"
#include "fts_sliding.h"

/* ============================================================
 * fts_sliding.c -- Sliding Mode Control for Finite-Time Stability
 *
 * Terminal Sliding Mode Control (TSMC):
 *   s(x) = x1 + beta * sig(x1)^alpha, 0 < alpha < 1
 *   where sig(x)^a = |x|^a * sign(x)
 *
 * The reaching law: ds/dt = -k * sign(s)
 *   - Reaching phase: s -> 0 in finite time T_r <= |s0|/k
 *   - Sliding phase: x1 -> 0 in finite time T_s <= |x1(tr)|^{1-alpha}/(beta*(1-alpha))
 *
 * Super-twisting algorithm (2nd-order SMC):
 *   u = -k1 * sig(s)^{1/2} + v
 *   dv/dt = -k2 * sign(s)
 * Provides finite-time convergence with continuous control.
 *
 * Refs: Utkin, Guldner & Shi (2009); Levant (2003); Feng et al. (2002)
 * ============================================================ */

/* --- Surface Management --- */

FTSSlidingSurface* fts_smc_surface_create(int dim, double gain,
                                           double exponent, int order) {
    FTSSlidingSurface* s = calloc(1, sizeof(FTSSlidingSurface));
    if (!s) return NULL;
    s->surface = calloc(dim, sizeof(double));
    if (!s->surface) { free(s); return NULL; }
    s->dim = dim;
    s->gain = gain;
    s->exponent = exponent;
    s->order = order;
    /* Terminal sliding surface: 0 < exponent < 1 */
    s->is_terminal = (exponent > 0.0 && exponent < 1.0);
    return s;
}

void fts_smc_surface_free(FTSSlidingSurface* s) {
    if (s) { free(s->surface); free(s); }
}

double fts_smc_eval_surface(FTSSlidingSurface* s, FTSState x) {
    /* Evaluate sliding surface: s(x) = x1 + beta * sig(x1)^alpha
     *
     * For terminal sliding mode, this surface ensures that once
     * s=0 is reached, the trajectory converges to x1=0 in finite time.
     *
     * General form: s = sum(c_i * x_i) + beta * |x1|^alpha * sign(x1) */
    if (!s) return 0.0;
    double sg = 0.0;
    for (int i = 0; i < x.dim && i < s->dim && i < s->order; i++)
        sg += x.x[i];
    sg += s->gain * pow(fabs(x.x[0]), s->exponent) *
          ((x.x[0] > 0) ? 1.0 : ((x.x[0] < 0) ? -1.0 : 0.0));
    return sg;
}

double fts_smc_surface_derivative(FTSSlidingSurface* s, FTSState x,
                                   FTSState dx) {
    /* Time derivative of sliding surface along system trajectory.
     *
     * ds/dt = grad(s)^T * f(x)
     *
     * For s = x1 + beta * |x1|^alpha * sign(x1):
     *   ds/dt = dx1/dt + beta*alpha*|x1|^{alpha-1} * dx1/dt
     *         = (1 + beta*alpha*|x1|^{alpha-1}) * dx1/dt
     *
     * For higher-order surfaces, include additional terms. */
    if (!s) return 0.0;
    double ds = 0.0;
    /* Linear part: sum of dx_i */
    for (int i = 0; i < x.dim && i < s->dim && i < s->order; i++)
        ds += dx.x[i];
    /* Nonlinear part for terminal sliding surface */
    if (s->is_terminal && fabs(x.x[0]) > 1e-12) {
        ds += s->gain * s->exponent *
              pow(fabs(x.x[0]), s->exponent - 1.0) * dx.x[0];
    }
    return ds;
}

/* --- SMC Controller --- */

FTSTerminalSMC* fts_smc_create(FTSSlidingSurface* ss, int dim_ctrl) {
    FTSTerminalSMC* smc = calloc(1, sizeof(FTSTerminalSMC));
    if (!smc) return NULL;
    smc->sigma = *ss;
    smc->dim_ctrl = dim_ctrl;
    smc->control = calloc(dim_ctrl, sizeof(double));
    if (!smc->control) { free(smc); return NULL; }
    smc->state_norm = NULL;
    smc->n_steps = 0;
    smc->reaching_time = 0.0;
    smc->sliding_time = 0.0;
    smc->settling_time = 0.0;
    smc->converged = false;
    smc->chattering = 0.0;
    smc->energy = 0.0;
    return smc;
}

void fts_smc_free(FTSTerminalSMC* smc) {
    if (smc) {
        free(smc->control);
        free(smc->state_norm);
        free(smc);
    }
}

int fts_smc_simulate(FTSTerminalSMC* smc, FTSSystem* sys,
                      FTSState x0, double T, double dt) {
    /* Simulate system with sliding mode controller.
     *
     * At each step:
     * 1. Evaluate sliding surface s(x)
     * 2. Apply discontinuous control: u = -K * sign(s)
     * 3. Integrate system dynamics
     * 4. Track reaching time, sliding time, chattering
     *
     * Returns: number of simulation steps. */
    if (!smc || !sys) return -1;
    FTSSystem* cp = fts_create(sys->dim, sys->np, dt);
    if (!cp) return -1;
    fts_copy_state(&cp->state, &x0);
    cp->rhs = sys->rhs;
    if (sys->np > 0)
        memcpy(cp->params, sys->params, sys->np * sizeof(double));

    /* Allocate norm tracking */
    int max_steps = (int)(T / dt);
    if (max_steps <= 0) max_steps = 1000;
    smc->state_norm = calloc(max_steps, sizeof(double));
    if (!smc->state_norm) { fts_free(cp); return -1; }

    double prev_s = 0.0;
    bool reached = false;
    int i;

    for (i = 0; i < max_steps; i++) {
        fts_step_rk4(cp);
        smc->n_steps++;

        /* Evaluate sliding surface */
        double s_val = fts_smc_eval_surface(&smc->sigma, cp->state);
        smc->state_norm[i] = fts_norm(cp->state);

        /* Track reaching time */
        if (!reached && fabs(s_val) < 1e-6 && i > 10) {
            smc->reaching_time = cp->t;
            reached = true;
        }

        /* Track sliding time (maintaining s~0) */
        if (reached && fabs(s_val) < 1e-6)
            smc->sliding_time += dt;

        /* Chattering measure: accumulated sign changes of s */
        if (i > 0 && prev_s * s_val < 0)
            smc->chattering += 1.0;

        /* Control energy */
        smc->energy += fabs(smc->control[0]) * dt;

        /* Check equilibrium convergence */
        if (fts_norm(cp->state) < 1e-6) {
            smc->settling_time = cp->t;
            smc->converged = true;
            break;
        }

        prev_s = s_val;
    }

    smc->n_steps = i + 1;
    if (!smc->converged) smc->settling_time = T;
    fts_free(cp);
    return smc->n_steps;
}

/* --- Control Laws --- */

double fts_smc_control(FTSTerminalSMC* smc, FTSState x) {
    /* Standard first-order SMC:
     * u = -K * sign(s(x))
     *
     * The discontinuous term sign(s) forces trajectories onto
     * the sliding surface s=0 in finite time (reaching phase),
     * then maintains sliding motion to the origin.
     *
     * In practice, sign(s) is approximated by sat(s/phi) to
     * reduce chattering, where phi is a boundary layer thickness. */
    if (!smc) return 0.0;
    double s = fts_smc_eval_surface(&smc->sigma, x);
    double u = -smc->sigma.gain * ((s > 0) ? 1.0 : ((s < 0) ? -1.0 : 0.0));
    return u;
}

double fts_smc_super_twisting_control(FTSState x, double lambda,
                                       double alpha) {
    /* Super-twisting algorithm (Levant, 1993):
     * u = -k1 * |s|^{1/2} * sign(s) + v
     * dv/dt = -k2 * sign(s)
     *
     * Key properties:
     * - Continuous control (no direct sign(s) in u)
     * - Finite-time convergence to s = ds/dt = 0
     * - Chattering-free in ideal sliding
     * - Requires only s measurement (not ds/dt)
     *
     * This function returns the static part.
     * The integral term v must be maintained by the caller. */
    double s = x.x[0];  /* assume s is x1 */
    return -lambda * pow(fabs(s), alpha) *
           ((s > 0) ? 1.0 : ((s < 0) ? -1.0 : 0.0));
}

double fts_smc_terminal_control(FTSState x, double k, double alpha,
                                 double beta) {
    /* Terminal SMC control law:
     * u = -k * |x1|^{alpha} * sign(x1) - beta * |x2| * sign(x2)
     *
     * Provides global finite-time stability when alpha in (0,1).
     * The fractional power near the origin accelerates convergence. */
    double u = -k * pow(fabs(x.x[0]), alpha) *
               ((x.x[0] > 0) ? 1.0 : ((x.x[0] < 0) ? -1.0 : 0.0));
    if (x.dim > 1) {
        u -= beta * fabs(x.x[1]) *
             ((x.x[1] > 0) ? 1.0 : ((x.x[1] < 0) ? -1.0 : 0.0));
    }
    return u;
}

double fts_smc_nonsingular_tsm_control(FTSState x, double k,
                                        double alpha, double beta) {
    /* Nonsingular Terminal SMC (Feng et al., 2002):
     * Solves the singularity problem of conventional TSMC
     * where the control grows unbounded when x1->0, x2!=0.
     *
     * NTSM surface: s = x1 + (1/beta) * |x2|^{p/q} * sign(x2)
     * where p,q are odd integers with p>q.
     *
     * This formulation avoids the singularity because
     * the fractional power is on x2, not x1. */
    double s = x.x[0];
    if (x.dim > 1 && fabs(x.x[1]) > 1e-12) {
        s += (1.0 / beta) *
             pow(fabs(x.x[1]), alpha / beta) *
             ((x.x[1] > 0) ? 1.0 : -1.0);
    }
    return -k * ((s > 0) ? 1.0 : ((s < 0) ? -1.0 : 0.0));
}

/* --- Bounds --- */

double fts_smc_reaching_time_bound(FTSSlidingSurface* s,
                                    FTSState x0, double gain) {
    /* Reaching time bound for first-order SMC:
     * T_r <= |s(x0)| / K
     *
     * Derivation: ds/dt = -K * sign(s)
     * => s(t) = s0 - K*t*sign(s0) for t < |s0|/K
     * => s(T_r) = 0 at T_r = |s0|/K
     *
     * If gain <= 0, returns infinity. */
    if (!s || gain < 1e-12) return 1e99;
    double s0 = fts_smc_eval_surface(s, x0);
    return fabs(s0) / gain;
}

double fts_smc_sliding_time_bound(FTSSlidingSurface* s, FTSState x0) {
    /* Sliding time bound for terminal sliding surface:
     * T_s <= |x1(tr)|^{1-alpha} / (beta*(1-alpha))
     *
     * where x1(tr) is the state when reaching s=0.
     * This uses the terminal attractor property of the surface.
     *
     * Note: x0 here represents the state at the start of sliding,
     * not the initial condition. */
    if (!s || !s->is_terminal) return 1e99;
    double x1_val = fabs(x0.x[0]);
    if (x1_val < 1e-12) return 0.0;
    double beta = s->gain;
    double alpha = s->exponent;
    if (beta < 1e-12 || alpha >= 1.0) return 1e99;
    return pow(x1_val, 1.0 - alpha) / (beta * (1.0 - alpha));
}

double fts_smc_settling_time_estimate(FTSTerminalSMC* smc) {
    /* Total settling time: reaching + sliding.
     * Only valid if convergence was achieved during simulation. */
    if (!smc || !smc->converged) return 1e99;
    return smc->reaching_time + smc->settling_time;
}

/* --- Analysis --- */

bool fts_smc_is_converged(FTSTerminalSMC* smc, double tol) {
    if (!smc) return false;
    return smc->converged && (smc->settling_time < 1e99);
}

bool fts_smc_is_nonsingular(FTSSlidingSurface* s, FTSState x) {
    /* Check if the terminal sliding surface is nonsingular.
     *
     * Conventional TSMC has a singularity: when x1=0 and x2!=0,
     * the control requires infinite effort.
     *
     * The nonsingular TSMC avoids this by placing the fractional
     * power on the velocity term. */
    if (!s || !s->is_terminal) return true;
    /* Simple heuristic: if near x1=0 with nonzero velocity,
     * conventional TSMC is singular */
    if (fabs(x.x[0]) < 1e-6 && x.dim > 1 && fabs(x.x[1]) > 1e-6)
        return false;
    return true;
}

double fts_smc_control_effort(FTSTerminalSMC* smc) {
    /* Total control effort: integral of |u(t)| dt
     * This measures the energy spent on control. */
    return smc ? smc->energy : 0.0;
}

double fts_smc_chattering_measure(FTSTerminalSMC* smc) {
    /* Chattering measure: number of sign changes of s(t)
     * normalized by total simulation steps.
     *
     * High chattering indicates excessive switching,
     * which can damage actuators in practice.
     *
     * Solutions:
     * - Boundary layer (sat instead of sign)
     * - Higher-order SMC (super-twisting)
     * - Observer-based SMC */
    if (!smc || smc->n_steps <= 0) return 0.0;
    return smc->chattering / (double)smc->n_steps;
}

double fts_smc_energy(FTSTerminalSMC* smc) {
    return smc ? smc->energy : 0.0;
}

/* --- Output --- */

void fts_smc_print(FTSTerminalSMC* smc) {
    if (!smc) return;
    printf("SMC Simulation Results:\n");
    printf("  Reaching time: %.6f\n", smc->reaching_time);
    printf("  Sliding time:  %.6f\n", smc->sliding_time);
    printf("  Settling time: %.6f\n", smc->settling_time);
    printf("  Converged: %s\n", smc->converged ? "yes" : "no");
    printf("  Chattering: %.6f\n",
           fts_smc_chattering_measure(smc));
    printf("  Control energy: %.6f\n", smc->energy);
    printf("  Steps: %d\n", smc->n_steps);
}

void fts_smc_print_surface(FTSSlidingSurface* s) {
    if (!s) return;
    printf("Sliding Surface:\n");
    printf("  Type: %s\n", s->is_terminal ? "terminal" : "linear");
    printf("  Gain: %.4f\n", s->gain);
    printf("  Exponent: %.4f\n", s->exponent);
    printf("  Order: %d\n", s->order);
}

/* --- Adaptive SMC --- */

void fts_smc_adapt_gain(FTSSlidingSurface* s, FTSState x,
                         double s_val, double adaptation_rate) {
    /* Adaptive SMC gain tuning.
     *
     * K_dot = gamma * |s|  when s != 0
     *
     * This increases the gain when away from the surface,
     * and maintains (or slowly decreases) when on the surface.
     * Ensures reaching without requiring prior knowledge of
     * the disturbance bound. */
    if (!s) return;
    if (fabs(s_val) > 1e-6) {
        s->gain += adaptation_rate * fabs(s_val);
    }
    /* Avoid gain going below minimum */
    if (s->gain < 0.1) s->gain = 0.1;
}

/* --- Observer-Based SMC --- */

double fts_smc_observer_gain(double lambda, double L) {
    /* Observer gain for sliding mode observer:
     * k = lambda * L, where L is the Lipschitz constant.
     * Must satisfy: k > L for convergence. */
    return lambda * L;
}

void fts_smc_levants_differentiator(double* x_hat, double* x,
                                     double dt, int order, double L) {
    /* Levant's robust exact differentiator (2003):
     *
     * Estimates derivatives in finite time using high-order SMC.
     *
     * Order 1: x0_hat_dot = -lambda0*|x0_hat-x|^{1/2}*sign(x0_hat-x) + x1_hat
     *          x1_hat_dot = -lambda1*sign(x1_hat - x0_hat_dot)
     *
     * Order n: similarly with nested structure.
     *
     * Parameters lambda_i must satisfy:
     *   lambda0 > sqrt(2L)
     *   lambda1 > L + (lambda0)^2/2
     *
     * This implements a simplified 1st-order differentiator. */
    if (!x_hat || !x) return;
    if (order >= 1) {
        double e = x_hat[0] - x[0];
        double lambda0 = 1.5 * sqrt(L);
        double lambda1 = 1.1 * L;
        /* x0_hat update */
        x_hat[0] += dt * (-lambda0 * sqrt(fabs(e)) *
                          ((e > 0) ? 1.0 : ((e < 0) ? -1.0 : 0.0)) +
                          x_hat[1]);
        /* x1_hat update */
        x_hat[1] += dt * (-lambda1 *
                          ((e > 0) ? 1.0 : ((e < 0) ? -1.0 : 0.0)));
    }
}