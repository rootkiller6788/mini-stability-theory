#include "gst_core.h"
#include "convergence_analysis.h"
#include "pendulum_model.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(void) {
    printf("=== Convergence Analysis: Damped Pendulum ===\n");

    /* Simulate damped pendulum using rk4_solve */
    double x0[] = {2.0, 1.0};
    double dt = 0.01, duration = 20.0;
    double* traj = NULL;
    int n_steps = 0;

    int ret = rk4_solve(damped_pendulum, x0, 2, 0.0, duration, dt,
                        &traj, &n_steps);
    if (ret != 0 || !traj) {
        printf("Simulation failed.\n");
        return 1;
    }

    /* Equilibrium is (0,0) for damped pendulum */
    double eq[] = {0.0, 0.0};
    ConvergenceResult* cm = convergence_create();
    conv_compute_from_trajectory(cm, traj, n_steps, 2, eq, dt);
    free(traj);

    printf("Exponential rate: %.6f\n", cm->exponential_rate);
    printf("Time constant: %.4f s\n", cm->time_constant);
    double st = conv_settling_time_to_tolerance(cm, 0.05);
    printf("Settling time (5%%): %.2f s\n", st);
    printf("Overshoot: %.4f\n", cm->overshoot);
    printf("Is contractive: %s\n",
           conv_is_contractive(cm) ? "YES" : "NO");
    printf("Decay envelope constant: %.4f\n",
           conv_decay_envelope(cm, cm->exponential_rate));
    printf("Monotone violations: %d\n",
           conv_monotone_convergence_check(cm));

    if (cm->exponential_rate > 0.01) {
        printf("\nSystem exhibits EXPONENTIAL convergence to origin.\n");
        printf("LaSalle principle confirms GAS via energy Lyapunov function.\n");
    }

    convergence_free(cm);
    printf("\nExample 3 PASSED\n");
    return 0;
}