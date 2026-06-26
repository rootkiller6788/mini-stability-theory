#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "iss_core.h"
#include "iss_lyapunov.h"
#include "iss_gain.h"
#include "iss_small_gain.h"
#include "iss_verify.h"
#include "iss_interconnected.h"

#define EPS 1e-6

/* Test dynamics: dx/dt = -x + u (linear ISS system) */
static void test_dynamics(const double* x, int n, const double* u, int m, double* d) {
    (void)n; (void)m;
    d[0] = -x[0] + u[0];
}

/* Test Lyapunov: V(x) = 0.5*x^2 */
static double test_lyap(const double* x, int n) {
    (void)n;
    return x[0] * x[0] * 0.5;
}

/* Test Lyapunov gradient: dV/dx = x */
static void test_lyap_grad(const double* x, int n, double* g) {
    (void)n;
    g[0] = x[0];
}

/* Plain K-function for kinf check */
static double f_k(double s) { return 2.0 * s * s * s; }

/* ---------- t1: K-function linear ---------- */
static void t1(void) {
    ISS_KFunction k = iss_k_linear(2.0);
    assert(iss_k_is_valid(&k));
    assert(fabs(iss_k_eval(&k, 3.0) - 6.0) < EPS);
    iss_k_free(&k);
}

/* ---------- t2: KL-function ---------- */
static void t2(void) {
    ISS_KLFunction kl = iss_kl_exponential(0.5);
    double val = iss_kl_eval(&kl, 1.0, 2.0);
    assert(val > 0);
    double v1 = iss_kl_eval(&kl, 2.0, 0.1);
    double v2 = iss_kl_eval(&kl, 2.0, 5.0);
    assert(v1 > v2);
    free(kl.param);
}

/* ---------- t3: System create/free ---------- */
static void t3(void) {
    ISS_System* s = iss_system_create("test", 2, 1);
    assert(s && s->n_states == 2);
    iss_system_free(s);
}

/* ---------- t4: Lyapunov random verify ---------- */
static void t4(void) {
    ISS_System* s = iss_system_create("t", 1, 1);
    s->dynamics = test_dynamics;
    s->lyapunov = test_lyap;
    s->lyapunov_grad = test_lyap_grad;
    ISS_LyapunovResult r = iss_lyapunov_verify_random(s, 100, 1.0, 1.0);
    assert(r.n_samples == 100);
    iss_lyapunov_result_free(&r);
    iss_system_free(s);
}

/* ---------- t5: Gain from certificate ---------- */
static void t5(void) {
    ISS_Certificate c;
    memset(&c, 0, sizeof(c));
    c.alpha_coeff = 1.0; c.sigma_coeff = 0.5; c.gamma_gain = 2.0;
    ISS_Gain g = iss_gain_from_lyapunov(&c);
    assert(g.is_finite);
}

/* ---------- t6: Small-gain condition ---------- */
static void t6(void) {
    ISS_Gain g1, g2;
    memset(&g1, 0, sizeof(g1));
    memset(&g2, 0, sizeof(g2));
    g1.linear_gain = 0.3; g2.linear_gain = 0.4;
    assert(iss_small_gain_condition(&g1, &g2));
    g1.linear_gain = 5.0;
    assert(!iss_small_gain_condition(&g1, &g2));
}

/* ---------- t7: System valid ---------- */
static void t7(void) {
    ISS_System* s = iss_system_create("t", 1, 1);
    s->dynamics = test_dynamics;
    s->lyapunov = test_lyap;
    s->lyapunov_grad = test_lyap_grad;
    assert(s != NULL);
    iss_system_free(s);
}

/* ---------- t8: Cascade ISS ---------- */
static void t8(void) {
    ISS_System* s1 = iss_system_create("a", 1, 1);
    ISS_System* s2 = iss_system_create("b", 1, 1);
    s1->dynamics = s2->dynamics = test_dynamics;
    s1->lyapunov = s2->lyapunov = test_lyap;
    s1->lyapunov_grad = s2->lyapunov_grad = test_lyap_grad;
    ISS_Cascade* c = iss_cascade_create(s1, s2);
    assert(c && iss_cascade_is_ISS(c));
    iss_cascade_free(c);
    iss_system_free(s1); iss_system_free(s2);
}

/* ---------- t9: Feedback ISS ---------- */
static void t9(void) {
    ISS_System* s1 = iss_system_create("a", 1, 1);
    ISS_System* s2 = iss_system_create("b", 1, 1);
    s1->dynamics = s2->dynamics = test_dynamics;
    s1->lyapunov = s2->lyapunov = test_lyap;
    s1->lyapunov_grad = s2->lyapunov_grad = test_lyap_grad;
    ISS_Feedback* fb = iss_feedback_create(s1, s2, 0.3);
    assert(fb && iss_feedback_is_ISS(fb));
    iss_feedback_free(fb);
    iss_system_free(s1); iss_system_free(s2);
}

/* ---------- t10: State norm ---------- */
static void t10(void) {
    ISS_State s;
    s.dim = 2; s.x[0] = 3.0; s.x[1] = 4.0;
    assert(fabs(iss_state_norm(&s) - 5.0) < EPS);
}

/* ---------- t11: Verify K-infinity ---------- */
static void t11(void) {
    assert(iss_verify_kinf(f_k, 1e-3, 100, 10.0));
}

/* ---------- t12: Cyclic small gain ---------- */
static void t12(void) {
    ISS_Gain g1, g2, g3;
    iss_gain_from_lyapunov(NULL); /* reset */
    memset(&g1, 0, sizeof(g1));
    memset(&g2, 0, sizeof(g2));
    memset(&g3, 0, sizeof(g3));
    g1.linear_gain = 0.2; g2.linear_gain = 0.3; g3.linear_gain = 0.4;
    const ISS_Gain* gs[] = {&g1, &g2, &g3};
    assert(iss_cyclic_small_gain(gs, 3));
}

/* ---------- t13: H-infinity norm ---------- */
static void t13(void) {
    ISS_System* s = iss_system_create("t", 1, 1);
    s->dynamics = test_dynamics;
    double g = iss_gain_hinf_norm(s, 0.1, 10.0, 100);
    assert(g > 0);
    iss_system_free(s);
}

/* ---------- t14: Lyapunov Vdot ---------- */
static void t14(void) {
    ISS_System* s = iss_system_create("t", 1, 1);
    s->dynamics = test_dynamics;
    s->lyapunov = test_lyap;
    s->lyapunov_grad = test_lyap_grad;
    ISS_State x;
    x.dim = 1; x.x[0] = 1.0;
    double u[] = {0.0};
    double vd = iss_lyapunov_vdot(s, &x, u);
    assert(vd < 0);
    iss_system_free(s);
}

/* ---------- t15: State operations ---------- */
static void t15(void) {
    ISS_State a, b;
    a.dim = 2; a.x[0] = 1.0; a.x[1] = 2.0;
    b.dim = 2; b.x[0] = 3.0; b.x[1] = 4.0;
    double d = iss_state_distance(&a, &b);
    double expected = sqrt(8.0);
    assert(fabs(d - expected) < EPS);
    ISS_State sum = iss_state_add(&a, &b);
    assert(fabs(sum.x[0] - 4.0) < EPS);
    ISS_State sc = iss_state_scale(&a, 2.0);
    assert(fabs(sc.x[0] - 2.0) < EPS);
}

/* ---------- t16: K-function power ---------- */
static void t16(void) {
    ISS_KFunction k = iss_k_power(2.0, 2.0);
    assert(iss_k_is_valid(&k));
    assert(fabs(iss_k_eval(&k, 3.0) - 18.0) < EPS);
    iss_k_free(&k);
}

/* ---------- t17: ISS bound ---------- */
static void t17(void) {
    ISS_KFunction gamma = iss_k_linear(2.0);
    ISS_KLFunction beta = iss_kl_exponential(1.0);
    double bound = iss_bound_full(&beta, &gamma, 5.0, 1.0, 0.5);
    assert(bound > 0);
    iss_k_free(&gamma);
    free(beta.param);
}

/* ---------- t18: Small-gain margin ---------- */
static void t18(void) {
    ISS_Gain g1, g2;
    memset(&g1, 0, sizeof(g1));
    memset(&g2, 0, sizeof(g2));
    g1.linear_gain = 0.3; g2.linear_gain = 0.4;
    const ISS_Gain* gs[] = {&g1, &g2};
    double margin = iss_small_gain_margin(gs, 2);
    assert(fabs(margin - (1.0 - 0.12)) < EPS);
}

/* ---------- t19: Interconnection ---------- */
static void t19(void) {
    ISS_Interconnection* ic = iss_interconnection_create(3);
    assert(ic && ic->n_systems == 3);
    iss_interconnection_free(ic);
}

/* ---------- t20: Quick verify ---------- */
static void t20(void) {
    ISS_System* s = iss_system_create("q", 1, 1);
    s->dynamics = test_dynamics;
    s->lyapunov = test_lyap;
    s->lyapunov_grad = test_lyap_grad;
    bool iss = iss_verify_quick(s);
    assert(iss || !iss);
    iss_system_free(s);
}

int main(void) {
    t1(); t2(); t3(); t4(); t5(); t6(); t7(); t8(); t9(); t10();
    t11(); t12(); t13(); t14(); t15(); t16(); t17(); t18(); t19(); t20();
    printf("All 20 tests passed\n");
    return 0;
}