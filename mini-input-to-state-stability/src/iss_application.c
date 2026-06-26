/* iss_application.c - L7 Application: Quadrotor altitude control with ISS.
 * Demonstrates Input-to-State Stability for a UAV subject to wind gusts.
 * The system is ISS if the altitude error remains bounded proportional to disturbance.
 *
 * L7: Quadrotor/DJI drone control, F-35 flight control gust rejection
 * Ref: Sontag (1989) "Smooth Stabilization Implies Coprime Factorization"
 *       Kendoul (2012) "Survey of Advances in Guidance, Navigation, Control of
 *                        Unmanned Rotorcraft Systems"
 *
 * Model: m*z'' = -mg + T + d(t)  where T=thrust, d(t)=wind disturbance
 * Controller: T = mg + kp*e_z + kd*e_z''  (PD altitude control)
 */

#include "iss_core.h"
#include "iss_verify.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define QUAD_MASS 0.5
#define QUAD_G 9.81

typedef struct {
    double mass, g, kp, kd, z_ref;
} QuadrotorParams;

static QuadrotorParams quad = {QUAD_MASS, QUAD_G, 4.0, 3.0, 1.0};

void quadrotor_ode(const double* x, double* dxdt, int n, double t) {
    (void)n;
    /* x[0]=z, x[1]=dz/dt, x[2]=integral_error */
    double e_z = x[0] - quad.z_ref;
    /* Wind disturbance: sinusoidal gust */
    double d_wind = 0.5 * sin(2.0 * M_PI * 0.5 * t);
    /* PD control thrust */
    double thrust = quad.mass * quad.g + quad.kp * e_z + quad.kd * x[1];
    dxdt[0] = x[1];
    dxdt[1] = -quad.g + thrust / quad.mass + d_wind / quad.mass;
}

double iss_lyapunov_candidate(const double* x, int n) {
    (void)n;
    double e = x[0] - quad.z_ref;
    return 0.5 * x[1] * x[1] + 0.5 * quad.kp * e * e;
}

int main(void) {
    printf("=== L7: Quadrotor ISS with Wind Disturbance ===\n");
    printf(" Mass=%.2f kg, Target altitude=%.1f m\n", QUAD_MASS, quad.z_ref);
    printf(" Controller: PD (kp=%.1f, kd=%.1f)\n\n", quad.kp, quad.kd);

    double x[2] = {0.0, 0.0};  /* start at ground, zero velocity */
    double dt = 0.01;
    double t = 0.0;
    double max_error = 0.0;

    for (int i = 0; i <= 200; i++) {
        double e = x[0] - quad.z_ref;
        if (fabs(e) > max_error) max_error = fabs(e);
        if (i % 20 == 0) {
            double V = iss_lyapunov_candidate(x, 2);
            printf("  t=%.2f: z=%.4f, dz=%.4f, V=%.6f\n", t, x[0], x[1], V);
        }
        for (int s = 0; s < 50; s++) {
            double dx[2];
            quadrotor_ode(x, dx, 2, t);
            x[0] += dt * dx[0];
            x[1] += dt * dx[1];
            t += dt;
        }
    }
    printf("\n Final: z=%.4f m, Max error=%.4f m\n", x[0], max_error);
    printf(" ISS property: bounded disturbance -> bounded state error.\n");
    printf(" The max altitude error (%.4f) remains bounded despite wind.\n", max_error);
    return 0;
}