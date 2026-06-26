#include "intercon_core.h"
#include "intercon_small_gain.h"
#include "intercon_passivity.h"
#include "intercon_lyapunov.h"
#include "intercon_decentralized.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define EPS 1e-9

int main(void) {
    /* === Subsystem creation === */
    ICSubsystem* s1 = ic_subsys_create(2, 1, 1, false);
    assert(s1);
    assert(s1->n == 2);
    double A1[] = {-2.0, 0.0, 0.0, -3.0};
    ic_subsys_set_A(s1, A1);
    double B1[] = {1.0, 0.0};
    ic_subsys_set_B(s1, B1);
    double C1[] = {1.0, 0.0};
    ic_subsys_set_C(s1, C1);
    double D1[] = {0.1};
    ic_subsys_set_D(s1, D1);
    double x0[] = {0.5, -0.3};
    ic_subsys_set_state(s1, x0);
    assert(fabs(s1->A[0][0] + 2.0) < EPS);

    ICSubsystem* s2 = ic_subsys_create(1, 1, 1, false);
    assert(s2);
    double A2[] = {-1.0};
    ic_subsys_set_A(s2, A2);
    double B2[] = {1.0};
    ic_subsys_set_B(s2, B2);
    double C2[] = {2.0};
    ic_subsys_set_C(s2, C2);
    double D2[] = {0.0};
    ic_subsys_set_D(s2, D2);
    double x02[] = {0.2};
    ic_subsys_set_state(s2, x02);
    assert(s2->n == 1);

    /* Freq response */
    double mag, phase;
    ic_subsys_freq_response(s1, 1.0, &mag, &phase);
    assert(mag >= 0.0);

    /* Step */
    double u[] = {1.0};
    ic_subsys_step(s1, u, 0.01);
    assert(fabs(s1->x[0]) > 0.0);
    double y[1];
    ic_subsys_output(s1, u, y);

    /* === Interconnection === */
    ICInterconnection* ic = ic_intercon_create(3);
    assert(ic);
    ic_intercon_set(ic, 0, 1, 0.5);
    assert(fabs(ic_intercon_get(ic, 0, 1) - 0.5) < EPS);
    double spec_rad;
    ic_intercon_spectral_radius(ic, &spec_rad);
    assert(spec_rad >= 0.0);
    ic_intercon_free(ic);

    /* === Interconnected System === */
    ICSystem* sys = ic_system_create();
    assert(sys);
    int idx1 = ic_system_add_subsys(sys, s1);
    int idx2 = ic_system_add_subsys(sys, s2);
    assert(idx1 == 0 && idx2 == 1);
    assert(sys->n_subsys == 2);
    assert(sys->total_states == 3);
    bool conn = ic_system_connect(sys, IC_FEEDBACK_NEG, 0, 1, 0.5);
    assert(conn);

    double r[] = {0.0, 0.0};
    ic_system_step(sys, r, 0.01);
    assert(sys->t > 0.0);

    double comp_state[IC_MAX_STATES];
    ic_system_get_composite_state(sys, comp_state);

    /* Composite A */
    double A_comp[IC_MAX_STATES * IC_MAX_STATES];
    ic_build_composite_A(sys, A_comp, IC_MAX_STATES);

    ICStability stab = ic_composite_eigenvalue_stability(sys);
    assert(stab == IC_ASYM_STABLE || stab == IC_MARG_STABLE || stab == IC_UNSTABLE);

    /* === Small-Gain === */
    SmallGainResult* sg = ic_small_gain_test(s1, s2, IC_FEEDBACK_NEG);
    assert(sg);
    assert(sg->product >= 0.0);
    ic_small_gain_result_print(sg);
    ic_small_gain_result_free(sg);

    SmallGainResult* sg_multi = ic_small_gain_multiloop(sys);
    assert(sg_multi);
    ic_small_gain_result_free(sg_multi);

    double opt_scale;
    double scaled = ic_scaled_small_gain(s1, s2, &opt_scale);
    assert(scaled >= 0.0);

    /* Circle criterion */
    CircleCriterion* cc = ic_circle_criterion(s1, 0.1, 2.0);
    assert(cc);
    ic_circle_criterion_print(cc);
    ic_circle_criterion_free(cc);

    /* Popov criterion */
    bool popov = ic_popov_criterion(s1, 1.0);
    (void)popov;

    /* ISS small-gain */
    double gains[IC_MAX_SUBSYS][IC_MAX_SUBSYS] = {{0.1, 0.3}, {0.2, 0.1}};
    ISS_SmallGain* iss = ic_iss_small_gain(gains, 2);
    assert(iss);
    ic_iss_small_gain_print(iss);
    ic_iss_small_gain_free(iss);

    /* === Passivity === */
    PassivityResult* pr1 = ic_passivity_test(s1);
    assert(pr1);
    ic_passivity_result_print(pr1);
    ic_passivity_result_free(pr1);

    PassivityResult* pr2 = ic_passivity_test(s2);
    assert(pr2);
    ic_passivity_result_free(pr2);

    bool pf_stable = ic_passivity_feedback_stable(s1, s2);
    (void)pf_stable;

    double pidx = ic_passivity_index(s2);
    assert(!isnan(pidx));

    double ofb = ic_output_feedback_passivity(s1);
    assert(ofb >= 0.0);

    double P_storage[IC_MAX_DIM * IC_MAX_DIM];
    ic_storage_function(s1, P_storage, IC_MAX_DIM);
    double av_stor = ic_available_storage(s1);
    assert(av_stor >= 0.0);

    double req_sup = ic_required_supply(s1, s1->x);
    assert(req_sup >= 0.0);

    /* Dissipativity */
    double Q[] = {-1.0}, S[] = {0.5}, R[] = {2.0};
    SupplyRate* sr = ic_supply_rate_create(1, Q, S, R);
    assert(sr);
    DissipativityResult* dr = ic_dissipativity_test(s1, sr);
    assert(dr);
    ic_dissipativity_result_print(dr);
    ic_dissipativity_result_free(dr);
    ic_supply_rate_free(sr);

    /* === Composite Lyapunov === */
    double weights[] = {1.0, 1.5};
    CompositeLyapunov* cl = ic_composite_lyapunov(sys, weights);
    assert(cl);
    assert(cl->V_composite >= 0.0);
    ic_composite_lyapunov_print(cl);
    ic_composite_lyapunov_free(cl);

    /* === M-Matrix === */
    MMatrixResult* mr = ic_build_M_matrix(sys);
    assert(mr);
    assert(mr->n == 2);
    ic_M_matrix_print(mr);
    ic_M_matrix_free(mr);

    /* Gershgorin */
    double A_test[] = {-5.0, 1.0, 2.0, -4.0};
    GershgorinResult* gr = ic_gershgorin_test(A_test, 2);
    assert(gr);
    assert(gr->all_left_half);
    ic_gershgorin_print(gr);
    ic_gershgorin_free(gr);

    assert(ic_is_row_diag_dominant(A_test, 2));
    assert(!ic_is_row_diag_dominant((double[]){-1.0, 2.0, 2.0, -1.0}, 2));
    assert(!ic_is_generalized_diag_dominant((double[]){-1.0, 3.0, 3.0, -1.0}, 2));

    /* Comparison system */
    double M_test[] = {-2.0, 0.5, 0.3, -1.5};
    double V0[] = {1.0, 0.5}, Vf[2];
    ic_comparison_system_simulate(M_test, 2, V0, 1.0, 0.01, Vf);
    assert(Vf[0] < V0[0]); /* Should decay */

    /* Connective stability */
    bool conn_stable = ic_connective_stability_test(sys);
    (void)conn_stable;

    /* Stability margins */
    double mu = ic_stability_margin_mu(sys);
    assert(!isnan(mu));
    double gm = ic_gain_margin(s1, s2);
    assert(gm >= 0.0);
    double pm = ic_phase_margin(s1, s2);
    assert(pm >= 0.0);

    /* === Decentralized === */
    DecentralizedController* dc = ic_decentralized_controller_create(2, 1);
    assert(dc);
    double K_gain[] = {0.5, 0.3};
    ic_decentralized_controller_set_gain(dc, K_gain);
    bool dc_stable = ic_decentralized_stabilize_check(dc, s1);
    (void)dc_stable;
    ic_decentralized_controller_free(dc);

    DecentralizedStability* ds = ic_decentralized_stability_test(sys);
    assert(ds);
    assert(ds->interaction_bound >= 0.0);
    ic_decentralized_stability_print(ds);
    ic_decentralized_stability_free(ds);

    double ib = ic_interaction_norm_bound(sys);
    assert(ib >= 0.0);

    OverlappingDecomp* od = ic_overlapping_decompose(A_comp, 3);
    assert(od);
    ic_overlapping_decomp_print(od);
    ic_overlapping_decomp_free(od);

    int fm = ic_fixed_modes_count(sys);
    assert(fm >= 0);
    bool dc_ctrl = ic_is_decentralized_controllable(sys);
    (void)dc_ctrl;
    bool dc_obs = ic_is_decentralized_observable(sys);
    (void)dc_obs;

    SparseHinfController* sc = ic_sparse_hinf_design(sys, 1.0);
    assert(sc);
    ic_sparse_hinf_print(sc);
    ic_sparse_hinf_free(sc);

    /* Matrix utilities */
    double A_mat[] = {1,2,3,4}, B_mat[] = {5,6,7,8}, C_mat[4];
    ic_matrix_mult(A_mat, B_mat, C_mat, 2, 2, 2);
    double D_mat[4];
    ic_matrix_add(A_mat, B_mat, D_mat, 2, 2);
    double E_mat[4];
    ic_matrix_scale(A_mat, 2.0, E_mat, 2, 2);
    double n_inf = ic_matrix_norm_inf(A_mat, 2, 2);
    assert(n_inf > 0.0);
    double n2 = ic_matrix_norm_2(A_mat, 2, 2);
    assert(n2 > 0.0);

    /* Lyapunov solver */
    double A_lyap[] = {-2.0, 0.0, 0.0, -3.0};
    double Q_lyap[] = {1.0, 0.0, 0.0, 1.0};
    double P_lyap[4];
    ic_solve_lyapunov(A_lyap, Q_lyap, P_lyap, 2);

    /* M-matrix check */
    double M_test2[] = {-2.0, 0.0, 0.0, -3.0};
    assert(ic_is_M_matrix(M_test2, 2));

    /* Cleanup: s1, s2 are owned by sys */
    ic_system_free(sys);

    printf("All tests passed.\n");
    return 0;
}