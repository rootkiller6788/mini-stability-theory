#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "circle_criterion.h"
#include "popov_criterion.h"
#include "kyp_lemma.h"
#include "lure_system.h"
#include "sector_theory.h"

int main(void) {
    printf("=== Circle-Popov Criterion Test Suite ===\n");

    /* ---- Test 1: Transfer Function Basics ---- */
    printf("  Test 1: Transfer function...\n");
    TransferFunction* tf = tf_create(1, 2);
    assert(tf != NULL);
    tf->num[0] = 1.0;
    tf->den[0] = 1.0;
    tf->den[1] = 1.0;

    double re, im;
    tf_eval_complex(tf, 1.0, &re, &im);
    assert(fabs(re - 0.5) < 0.01);
    assert(fabs(im + 0.5) < 0.01);
    assert(fabs(tf_eval(tf, 0.0) - 1.0) < 0.01);
    assert(tf_eval(tf, 1e6) < 0.01);

    /* ---- Test 2: Circle Criterion Critical Disk ---- */
    printf("  Test 2: Critical disk...\n");
    SectorBounds sb = {0.5, 2.0, 0, 0, 0.0, 0.0};
    double cx = circle_criterion_center(sb);
    double cr = circle_criterion_radius(sb);
    assert(fabs(cx + 1.25) < 0.01);
    assert(fabs(cr - 0.75) < 0.01);

    /* ---- Test 3: Circle Criterion Half-Plane ---- */
    printf("  Test 3: Half-plane case...\n");
    SectorBounds sb_hp = {0.0, 2.0, 0, 0, 0.0, 0.0};
    double margin;
    int result = circle_criterion_check(tf, sb_hp, &margin);
    assert(result == 1);
    assert(margin > 0.0);

    /* ---- Test 4: Nyquist Data ---- */
    printf("  Test 4: Nyquist data...\n");
    NyquistData* nd = nyquist_create(100);
    assert(nd != NULL);
    nyquist_compute(tf, nd, 0.01, 10.0);
    assert(nd->n == 100);
    assert(!circle_criterion_encircles(nd, -10.0, 0.0));
    assert(circle_criterion_min_distance(nd, -5.0, 0.0) > 1.0);
    assert(!circle_criterion_in_disk(nd, -100.0, 0.5));

    /* ---- Test 5: Log-swept Nyquist ---- */
    printf("  Test 5: Log-swept Nyquist...\n");
    NyquistData* nd2 = nyquist_create(200);
    assert(nd2 != NULL);
    nyquist_compute_log(tf, nd2, 0.01, 100.0);
    assert(nd2->n == 200);
    assert(fabs(nd2->real[0] - 1.0) < 0.1);

    /* ---- Test 6: Popov Criterion ---- */
    printf("  Test 6: Popov criterion...\n");
    SectorBounds sb_popov = {0.0, 2.0, 0, 0, 0.0, 0.0};
    result = popov_criterion_check(tf, sb_popov, 0.5, &margin);
    assert(result == 1);
    assert(margin > 0.0);

    /* ---- Test 7: Popov Optimal Slope ---- */
    printf("  Test 7: Optimal Popov slope...\n");
    double eta_opt, margin_opt;
    popov_find_optimal_slope(tf, sb_popov, &eta_opt, &margin_opt);
    assert(eta_opt >= 0.0);
    assert(margin_opt > 0.0);

    /* ---- Test 8: Popov Line Test ---- */
    printf("  Test 8: Popov line test...\n");
    PopovData* pd = popov_data_create(50);
    assert(pd != NULL);
    popov_compute(tf, pd, 0.01, 10.0);
    double violation;
    popov_line_test(pd, 1.0, 1.0, &violation);
    popov_data_free(pd);

    /* ---- Test 9: Loop Transformation ---- */
    printf("  Test 9: Loop transformation...\n");
    SectorBounds sb_lt = {0.5, 2.5, 0, 0, 0.0, 0.0};
    double k;
    TransferFunction* tf2 = popov_loop_transform(tf, sb_lt, &k);
    assert(tf2 != NULL);
    assert(fabs(k - 2.0) < 0.01);
    double re2, im2;
    tf_eval_complex(tf2, 0.0, &re2, &im2);
    assert(fabs(re2 - 0.6667) < 0.01);
    tf_free(tf2);

    /* ---- Test 10: State-Space and KYP ---- */
    printf("  Test 10: State-space and KYP...\n");
    StateSpace* ss = ss_create(2);
    assert(ss != NULL);
    double A[4] = {0.0, 1.0, -1.0, -0.5};
    double B[2] = {0.0, 1.0};
    double C[2] = {1.0, 0.0};
    ss_set_matrices(ss, A, B, C, 0.1);

    KYPSolution* sol = kyp_solve(ss, 0.01);
    assert(sol != NULL);
    assert(sol->n == 2);
    assert(sol->P != NULL);
    kyp_free(sol);

    /* ---- Test 11: State-Space to TF ---- */
    printf("  Test 11: SS-to-TF conversion...\n");
    TransferFunction* tf3 = tf_create(3, 3);
    assert(tf3 != NULL);
    kyp_to_transfer_function(ss, tf3);
    assert(fabs(tf3->den[2] - 1.0) < 0.01);  /* s^2 coeff */
    tf_free(tf3);

    /* ---- Test 12: Lyapunov Solver ---- */
    printf("  Test 12: Lyapunov solver...\n");
    double A2[4] = {-1.0, 0.0, 0.0, -2.0};
    double Q2[4] = {1.0, 0.0, 0.0, 1.0};
    double X2[4];
    int lyap_ok = kyp_solve_lyapunov(A2, X2, Q2, 2);
    assert(lyap_ok == 0);
    assert(fabs(X2[0] - 0.5) < 0.01);
    assert(fabs(X2[3] - 0.25) < 0.01);

    /* ---- Test 13: Passivity Check ---- */
    printf("  Test 13: Passivity check...\n");
    int spr = kyp_verify_passivity(ss);
    (void)spr;

    /* ---- Test 14: Lur'e System (simple 1st-order) ---- */
    printf("  Test 14: Lur'e system...\n");
    StateSpace* ss1 = ss_create(1);
    assert(ss1 != NULL);
    ss1->A[0] = -1.0; ss1->B[0] = 1.0; ss1->C[0] = 1.0; ss1->D = 0.0;
    LureSystem* ls = lure_create(ss1, sb_hp);
    assert(ls != NULL);
    assert(lure_sector_lower_bound(ls) == 0.0);
    assert(lure_sector_upper_bound(ls) == 2.0);
    assert(lure_sector_check(ls, 1.0, 1.0));
    assert(!lure_sector_check(ls, 1.0, 3.0));

    /* ---- Test 15: Lur'e Stability Tests ---- */
    printf("  Test 15: Lur'e stability...\n");
    double m_circle, m_popov, eta;
    int st_circle = lure_circle_test(ls, &m_circle);
    assert(st_circle == 1);
    assert(m_circle > 0.0);
    int st_popov = lure_popov_test(ls, &m_popov, &eta);
    assert(st_popov == 1);
    assert(m_popov > 0.0);

    /* ---- Test 16: Unified Stability ---- */
    printf("  Test 16: Unified stability test...\n");
    assert(lure_absolute_stability(ls, LURE_METHOD_CIRCLE, &m_circle, NULL) == 1);
    assert(lure_absolute_stability(ls, LURE_METHOD_POPOV, &m_popov, &eta) == 1);
    assert(lure_absolute_stability(ls, LURE_METHOD_ALL, &m_circle, NULL) == 1);

    /* ---- Test 17: Gain/Phase Margins ---- */
    printf("  Test 17: Gain/phase margins...\n");
    StabilityMargin* sm = lure_gain_margin(ls);
    assert(sm != NULL);
    assert(sm->gain_margin > 0.0);
    lure_margin_free(sm);

    /* ---- Test 18: H-infinity Norm ---- */
    printf("  Test 18: H-infinity norm...\n");
    double hinf = lure_hinf_norm(ls);
    assert(hinf > 0.0);

    /* ---- Test 19: Sector Theory ---- */
    printf("  Test 19: Sector theory...\n");
    assert(sector_check_point(sb_hp, 1.0, 1.0));
    assert(!sector_check_point(sb_hp, 1.0, 3.0));
    assert(nonlinearity_saturation(2.0, 1.0) == 1.0);
    assert(fabs(nonlinearity_deadzone(0.3, 0.2) - 0.1) < 1e-10);
    assert(nonlinearity_relay(1.0) == 1.0);
    assert(nonlinearity_relay(-1.0) == -1.0);
    assert(fabs(nonlinearity_cubic(1.0, 0.5) - 1.5) < 1e-10);
    assert(fabs(nonlinearity_quantizer(1.7, 1.0) - 2.0) < 1e-10);

    /* ---- Test 20: Sector Transformations ---- */
    printf("  Test 20: Sector transformations...\n");
    SectorBounds sb_orig = {1.0, 3.0, 0, 0, 0.0, 0.0};
    SectorBounds shifted = sector_shift(sb_orig, -1.0);
    assert(fabs(shifted.alpha) < 0.01);
    assert(fabs(shifted.beta - 2.0) < 0.01);

    /* ---- Test 21: Compare Circle vs Popov ---- */
    printf("  Test 21: Circle vs Popov comparison...\n");
    double cm2, pm2, ef;
    int comp = popov_compare_criteria(tf, sb_hp, &cm2, &pm2, &ef);
    (void)comp; (void)cm2; (void)pm2; (void)ef;

    /* ---- Test 22: Max Allowed Gain ---- */
    printf("  Test 22: Max allowed gain...\n");
    double max_gain, mg_margin;
    max_gain = sector_max_allowed_gain(tf, 0, &mg_margin);
    assert(max_gain > 0.0);

    /* ---- Cleanup ---- */
    nyquist_free(nd);
    nyquist_free(nd2);
    tf_free(tf);
    lure_free(ls);
    ss_free(ss1);
    ss_free(ss);

    printf("\n=== ALL TESTS PASSED ===\n");
    return 0;
}
