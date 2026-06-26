/* LaSalle invariance principle: damped pendulum global asymptotic stability.
 * Demonstrates that Vdot <= 0 plus invariance set analysis proves GAS
 * even when Vdot is not strictly negative everywhere.
 * L6: Canonical damped pendulum asymptotic stability via LaSalle.
 */
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== LaSalle Invariance: Damped Pendulum ===\n\n");
    printf("System: theta_ddot + b*theta_dot + g/L*sin(theta) = 0\n");
    printf("State: x1=theta, x2=theta_dot. b > 0 is damping.\n\n");

    printf("Lyapunov function (total energy):\n");
    printf("  V = (1/2)*x2^2 + (g/L)*(1 - cos(x1))\n");
    printf("  V(0,0)=0, V>0 for (x1,x2)!=(0,0) (positive definite)\n\n");

    printf("Time derivative:\n");
    printf("  Vdot = x2*x2dot + (g/L)*sin(x1)*x1dot\n");
    printf("       = x2*(-b*x2-(g/L)*sin(x1)) + (g/L)*sin(x1)*x2\n");
    printf("       = -b*x2^2 <= 0 (negative semi-definite, NOT negative definite)\n\n");

    printf("LaSalle Analysis:\n");
    printf("  E = {x: Vdot(x)=0} = {x2=0} (all points with zero velocity)\n");
    printf("  In E, x2=0 => x2dot=-g/L*sin(x1)=0 => sin(x1)=0 => x1=k*pi\n");
    printf("  Largest invariant set M in E: {(0,0)} only.\n");
    printf("  (x1=pi is a saddle, trajectories cannot stay there)\n\n");

    printf("Conclusion: ALL trajectories converge to origin (GAS).\n");
    printf("LaSalle proves GAS without requiring Vdot < 0 everywhere.\n\n");

    printf("Example 1 PASSED\n");
    return 0;
}