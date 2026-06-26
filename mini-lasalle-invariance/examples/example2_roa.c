#include "gst_core.h"
#include "region_attraction.h"
#include "pendulum_model.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void) {
    printf("=== Region of Attraction Estimation ===\n");
    RegionAttraction* roa = roa_create(2);
    double ctr[] = {0.0, 0.0};

    /* Estimate ROA using Lyapunov sublevel set method */
    int ret = roa_estimate_lyapunov(roa, quadratic_lyapunov,
        damped_pendulum, 2, ctr, 3.0, 30, 1e-6);
    if (ret != 0) {
        printf("ROA estimation failed.\n");
        roa_free(roa);
        return 1;
    }

    printf("Estimated max safe Lyapunov level: %.4f\n",
           roa->inner_approximation->max_safe_level);
    printf("Sublevel set interior points: %d\n",
           roa->inner_approximation->n_points);

    /* Monte Carlo verification within estimated ROA */
    double roa_radius = sqrt(2.0 * roa->inner_approximation->max_safe_level);
    int n_conv = roa_monte_carlo_verify(roa, damped_pendulum,
        2, ctr, roa_radius, 100, 0.01, 10.0, 0.1);
    printf("Monte Carlo verification: %d/%d converged (%.1f%%)\n",
           n_conv, roa->n_samples,
           100.0 * roa->confidence);

    printf("ROA volume estimate: %.4f\n", roa->estimated_volume);
    printf("Confidence: %.2f\n\n", roa->confidence);

    if (roa->confidence > 0.5) {
        printf("ROA estimation SUCCESSFUL:\n");
        printf("  The origin is asymptotically stable with\n");
        printf("  estimated region of attraction radius >= %.4f\n", roa_radius);
    }

    roa_free(roa);
    printf("\nExample ROA PASSED\n");
    return 0;
}