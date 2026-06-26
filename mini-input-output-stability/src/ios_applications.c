#include "ios_core.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* ============================================================
 * Application Demos and Examples
 * Each function demonstrates a key concept of I/O stability.
 * ============================================================ */

/* --- BIBO Stability Demos --- */

void ios_example_first_order(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    printf("1/(s+1): BIBO=%s Hinf=%.2f\n",
           ios_check_bibo_tf(tf) == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE",
           ios_Hinf_norm_tf(tf));
    ios_tf_free(tf);
}

void ios_example_second_order(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 3);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 0.5; tf->den[2] = 1.0;
    printf("1/(s^2+0.5s+1): BIBO=%s\n",
           ios_check_bibo_tf(tf) == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE");
    ios_tf_free(tf);
}

void ios_example_small_gain(void) {
    IOSSmallGainResult r = ios_small_gain_test(0.5, 0.5);
    printf("Small-gain test (0.5*0.5): %s (margin=%.4f)\n",
           r.satisfied ? "PASS" : "FAIL", r.margin);
}

void ios_example_passivity(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    IOSPassivityResult r = ios_check_passivity_tf(tf);
    printf("1/(s+1): passive=%s diss=%.4f\n",
           r.passive ? "YES" : "NO", r.dissipation);
    ios_tf_free(tf);
}

void ios_example_popov(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 3);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 2.0; tf->den[2] = 1.0;
    printf("Popov test (k=1.0): %s\n",
           ios_popov_criterion_tf(tf, 1.0) ? "STABLE" : "INCONCLUSIVE");
    ios_tf_free(tf);
}

void ios_example_sector(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 2.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    IOSSector sec = ios_compute_sector_tf(tf);
    printf("2/(s+1) sector: [%.2f, %.2f]\n", sec.a, sec.b);
    ios_tf_free(tf);
}

void ios_example_robustness(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    printf("1/(s+1) disk margin: %.4f\n", ios_disk_margin(tf));
    ios_tf_free(tf);
}

void ios_example_all(void) {
    printf("=== I/O Stability Examples ===\n");
    ios_example_first_order();
    ios_example_second_order();
    ios_example_small_gain();
    ios_example_passivity();
    ios_example_popov();
    ios_example_sector();
    ios_example_robustness();
}

/* --- Standalone Demos --- */

void ios_demo_bibo(void) {
    IOSStateSpace* ss = ios_ss_create(2, 1, 1);
    if (!ss) return;
    ss->A[0] = -1.0; ss->A[3] = -2.0;
    ss->B[0] = 1.0; ss->C[0] = 1.0;
    printf("BIBO demo (2x2 stable): %s\n",
           ios_check_bibo_ss(ss) == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE");
    ios_ss_free(ss);
}

void ios_demo_smallgain_report(void) {
    IOSSmallGainResult r = ios_small_gain_test(0.6, 0.6);
    printf("SG demo: 0.6*0.6=%.4f %s (margin=%.4f)\n",
           r.gamma, r.satisfied ? "STABLE" : "FAIL", r.margin);
}

void ios_demo_passivity_report(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    IOSPassivityResult r = ios_check_passivity_tf(tf);
    printf("Passivity demo: %s diss=%.4f\n",
           r.passive ? "PASSIVE" : "NOT PASSIVE", r.dissipation);
    ios_tf_free(tf);
}

void ios_demo_circle_report(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    bool c = ios_circle_criterion_tf(tf, 0.0, 5.0);
    printf("Circle[0,5]: %s\n", c ? "STABLE" : "INCONCLUSIVE");
    ios_tf_free(tf);
}

void ios_demo_sector_report(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 2.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    IOSSector s = ios_compute_sector_tf(tf);
    printf("Sector demo: [%.4f, %.4f]\n", s.a, s.b);
    ios_tf_free(tf);
}

void ios_demo_bandwidth(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    double bw = ios_bandwidth(tf);
    printf("1/(s+1) bandwidth: %.4f rad/s\n", bw);
    ios_tf_free(tf);
}

void ios_demo_disk_margin(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    double dm = ios_disk_margin(tf);
    printf("1/(s+1) disk margin: %.4f\n", dm);
    ios_tf_free(tf);
}

void ios_demo_gain_margins(void) {
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0;
    tf->den[0] = 1.0; tf->den[1] = 1.0;
    double gm = 0.0, pm = 0.0;
    ios_bibo_margins(tf, &gm, &pm);
    printf("1/(s+1): GM=%.2f PM=%.2f\n", gm, pm * 180.0 / IOS_PI);
    ios_tf_free(tf);
}

void ios_demo_compare_gains(void) {
    printf("=== Gain Comparison ===\n");
    IOSTransferFcn* tf1 = ios_tf_create(1, 2);
    IOSTransferFcn* tf2 = ios_tf_create(1, 2);
    if (!tf1 || !tf2) { ios_tf_free(tf1); ios_tf_free(tf2); return; }
    tf1->num[0] = 1.0; tf1->den[0] = 1.0; tf1->den[1] = 1.0;
    tf2->num[0] = 2.0; tf2->den[0] = 1.0; tf2->den[1] = 0.5;
    printf("  G1=1/(s+1):   Hinf=%.4f H2=%.4f\n", ios_Hinf_norm_tf(tf1), ios_H2_norm_tf(tf1));
    printf("  G2=2/(s+0.5): Hinf=%.4f H2=%.4f\n", ios_Hinf_norm_tf(tf2), ios_H2_norm_tf(tf2));
    ios_tf_free(tf1); ios_tf_free(tf2);
}

void ios_demo_parametric_sweep(void) {
    printf("=== Parametric Sweep ===\n");
    for (double a = 1.0; a <= 4.0; a += 1.0) {
        IOSTransferFcn* tf = ios_tf_create(1, 2);
        if (!tf) continue;
        tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = a;
        printf("  a=%.1f: G(s)=1/(s+%.1f) -> %s\n", a, a,
               ios_check_bibo_tf(tf) == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE");
        ios_tf_free(tf);
    }
}

void ios_demo_stability_margins(void) {
    printf("=== Stability Margins ===\n");
    IOSTransferFcn* tf = ios_tf_create(1, 3);
    if (!tf) return;
    tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 0.5; tf->den[2] = 1.0;
    double gm = 0.0, pm = 0.0;
    ios_bibo_margins(tf, &gm, &pm);
    double dm = ios_disk_margin(tf);
    printf("  G(s)=1/(s^2+0.5s+1): GM=%.2f PM=%.1f deg DM=%.4f\n",
           gm, pm * 180.0 / IOS_PI, dm);
    ios_tf_free(tf);
}

void ios_demo_loop_shaping(void) {
    printf("=== Loop Shaping Demo ===\n");
    IOSTransferFcn* G = ios_tf_create(1, 2);
    if (!G) return;
    G->num[0] = 1.0; G->den[0] = 1.0; G->den[1] = 1.0;
    double Ms = ios_sensitivity_peak(G);
    double Mt = ios_complementary_sensitivity(G);
    printf("  Sensitivity peak Ms: %.4f\n", Ms);
    printf("  Complementary peak Mt: %.4f\n", Mt);
    printf("  Disk margin: %.4f\n", ios_disk_margin(G));
    ios_tf_free(G);
}

void ios_demo_saturation_analysis(void) {
    printf("=== Saturation Analysis ===\n");
    IOSTransferFcn* G = ios_tf_create(1, 2);
    if (!G) return;
    G->num[0] = 1.0; G->den[0] = 1.0; G->den[1] = 1.0;
    bool circle = ios_circle_criterion_tf(G, 0.0, 1.0);
    printf("  Circle [0,1] (saturation): %s\n", circle ? "STABLE" : "INCONCLUSIVE");
    bool popov = ios_popov_criterion_tf(G, 1.0);
    printf("  Popov [0,1]: %s\n", popov ? "STABLE" : "INCONCLUSIVE");
    ios_tf_free(G);
}

void ios_demo_interconnection(void) {
    printf("=== Interconnection Demo ===\n");
    IOSStateSpace* P = ios_ss_create(2, 1, 1);
    IOSStateSpace* C = ios_ss_create(1, 1, 1);
    if (!P || !C) { ios_ss_free(P); ios_ss_free(C); return; }
    P->A[0] = -1.0; P->A[1] = 0.0; P->A[2] = 1.0; P->A[3] = -2.0;
    P->B[0] = 1.0; P->B[1] = 0.0; P->C[0] = 1.0; P->C[1] = 0.0;
    C->A[0] = -5.0; C->B[0] = 1.0; C->C[0] = 1.0;
    double gP = ios_L2_gain_ss(P), gC = ios_L2_gain_ss(C);
    printf("  Plant L2 gain: %.4f\n", gP);
    printf("  Controller L2 gain: %.4f\n", gC);
    printf("  Small-gain: %s\n", gP * gC < 1.0 ? "STABLE" : "INCONCLUSIVE");
    ios_ss_free(P); ios_ss_free(C);
}

void ios_demo_uncertainty(void) {
    printf("=== Uncertainty Demo ===\n");
    IOSTransferFcn* Gnom = ios_tf_create(1, 2);
    IOSTransferFcn* Gpert = ios_tf_create(1, 2);
    if (!Gnom || !Gpert) { ios_tf_free(Gnom); ios_tf_free(Gpert); return; }
    Gnom->num[0] = 1.0; Gnom->den[0] = 1.0; Gnom->den[1] = 1.0;
    Gpert->num[0] = 1.1; Gpert->den[0] = 1.0; Gpert->den[1] = 0.9;
    double gap = ios_gap_metric(Gnom, Gpert);
    printf("  Gap metric: %.6f\n", gap);
    printf("  Disk margin (nom): %.4f\n", ios_disk_margin(Gnom));
    ios_tf_free(Gnom); ios_tf_free(Gpert);
}

void ios_demo_nyquist(void) {
    printf("=== Nyquist Demo ===\n");
    IOSTransferFcn* tf = ios_tf_create(1, 3);
    if (!tf) return;
    tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 0.5; tf->den[2] = 1.0;
    printf("G(s)=1/(s^2+0.5s+1): disk_margin=%.4f\n", ios_disk_margin(tf));
    ios_tf_free(tf);
}

void ios_demo_h2_hinf(void) {
    printf("=== H2/Hinf Demo ===\n");
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
    printf("1/(s+1): H2=%.4f Hinf=%.4f\n", ios_H2_norm_tf(tf), ios_Hinf_norm_tf(tf));
    ios_tf_free(tf);
}

void ios_demo_impulse_response(void) {
    IOSStateSpace* ss = ios_ss_create(2, 1, 1);
    if (!ss) return;
    ss->A[0] = -1.0; ss->A[3] = -2.0;
    ss->B[0] = 1.0; ss->C[0] = 1.0;
    printf("Impulse response energy: %.4f\n", ios_impulse_energy(ss));
    ios_ss_free(ss);
}

void ios_demo_step_response(void) {
    IOSStateSpace* ss = ios_ss_create(2, 1, 1);
    if (!ss) return;
    ss->A[0] = -1.0; ss->A[3] = -2.0;
    ss->B[0] = 1.0; ss->C[0] = 1.0;
    IOSSignal* y = ios_signal_create(100);
    if (y) {
        ios_step_response_full(ss, 5.0, 0.05, y);
        printf("Step response: %d points, final=%.4f\n", y->length, y->data[y->length - 1]);
        ios_signal_free(y);
    }
    ios_ss_free(ss);
}

void ios_demo_l1_l2_linf(void) {
    IOSStateSpace* ss = ios_ss_create(2, 1, 1);
    if (!ss) return;
    ss->A[0] = -1.0; ss->A[3] = -2.0;
    ss->B[0] = 1.0; ss->C[0] = 1.0;
    printf("Norms: L1=%.4f L2=%.4f Linf=%.4f\n",
           ios_L1_induced_norm(ss), ios_L2_gain_ss(ss), ios_Linf_induced_norm(ss));
    ios_ss_free(ss);
}

void ios_demo_all_metrics(void) {
    printf("=== All Metrics Demo ===\n");
    IOSStateSpace* ss = ios_ss_create(2, 1, 1);
    if (!ss) return;
    ss->A[0] = -1.0; ss->A[3] = -2.0;
    ss->B[0] = 1.0; ss->C[0] = 1.0;
    ios_print_gain_table(ss);
    ios_passivity_all_metrics(ss);
    ios_ss_free(ss);
}

void ios_demo_quick_check(void) {
    printf("=== Quick Check ===\n");
    IOSTransferFcn* tf = ios_tf_create(1, 2);
    if (!tf) return;
    tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
    printf("1/(s+1): stable=%s passive=%s Hinf=%.2f\n",
           ios_tf_is_stable(tf) ? "YES" : "NO",
           ios_check_passivity_tf(tf).passive ? "YES" : "NO",
           ios_Hinf_norm_tf(tf));
    ios_tf_free(tf);
}

/* --- Master Demos --- */

void ios_run_all_demos(void) {
    ios_example_all();
    ios_demo_compare_gains();
    ios_demo_parametric_sweep();
    ios_demo_stability_margins();
    ios_demo_loop_shaping();
    ios_demo_saturation_analysis();
    ios_demo_interconnection();
    ios_demo_uncertainty();
}

void ios_demo_final(void) {
    printf("=== FINAL DEMO ===\n");
    ios_run_all_demos();
    ios_demo_nyquist();
    ios_demo_h2_hinf();
    printf("Module complete.\n");
}

void ios_demo_final_report(void) {
    ios_demo_quick_check();
    ios_demo_all_metrics();
    ios_demo_final();
}

void ios_all_demos(void) {
    ios_demo_bibo();
    ios_example_all();
}

/* ============================================================
 * Comprehensive Demo: All Analysis Types
 * ============================================================ */

void ios_demo_comprehensive_analysis(void) {
    printf("\n========================================\n");
    printf("  COMPREHENSIVE I/O STABILITY ANALYSIS\n");
    printf("========================================\n\n");

    /* Test system: G(s) = 2/(s^2 + 0.5s + 1) */
    IOSTransferFcn* G = ios_tf_create(1, 3);
    if (!G) return;
    G->num[0] = 2.0;
    G->den[0] = 1.0; G->den[1] = 0.5; G->den[2] = 1.0;

    printf("Transfer function: G(s) = 2 / (s^2 + 0.5s + 1)\n\n");

    /* BIBO Stability */
    printf("--- BIBO Stability ---\n");
    IOSBIBOStatus bibo = ios_check_bibo_tf(G);
    printf("  Status: %s\n", bibo == IOS_BIBO_STABLE ? "STABLE" : "UNSTABLE");
    printf("  H2 norm:   %.6f\n", ios_H2_norm_tf(G));
    printf("  Hinf norm: %.6f\n", ios_Hinf_norm_tf(G));
    printf("  Bandwidth: %.4f rad/s\n", ios_bandwidth(G));

    /* Passivity */
    printf("\n--- Passivity Analysis ---\n");
    IOSPassivityResult pr = ios_check_passivity_tf(G);
    printf("  Passive:        %s\n", pr.passive ? "YES" : "NO");
    printf("  Strictly Input: %s\n", pr.strictly_input ? "YES" : "NO");
    printf("  Dissipation:    %.6f\n", pr.dissipation);

    /* Small-Gain */
    printf("\n--- Small-Gain Analysis ---\n");
    double gam = ios_Hinf_norm_tf(G);
    printf("  System gain gamma = %.4f\n", gam);
    printf("  Feedback stable for H with gain < %.4f\n", 1.0 / gam);

    /* Circle Criterion */
    printf("\n--- Circle Criterion ---\n");
    printf("  k1=0, k2=2:   %s\n",
           ios_circle_criterion_tf(G, 0.0, 2.0) ? "STABLE" : "INCONCLUSIVE");
    printf("  k1=0, k2=5:   %s\n",
           ios_circle_criterion_tf(G, 0.0, 5.0) ? "STABLE" : "INCONCLUSIVE");
    printf("  k1=0, k2=0.5: %s\n",
           ios_circle_criterion_tf(G, 0.0, 0.5) ? "STABLE" : "INCONCLUSIVE");

    /* Popov Criterion */
    printf("\n--- Popov Criterion ---\n");
    for (double k = 0.5; k <= 5.0; k *= 2.0) {
        printf("  Popov[0, %.1f]: %s\n", k,
               ios_popov_criterion_tf(G, k) ? "STABLE" : "INCONCLUSIVE");
    }

    /* Robustness */
    printf("\n--- Robustness ---\n");
    printf("  Disk margin:    %.6f\n", ios_disk_margin(G));
    printf("  Sensitivity Ms: %.6f\n", ios_sensitivity_peak(G));
    printf("  Nichols peak:   %.6f\n", ios_nichols_peak(G));

    /* Gap Metric (compare to perturbed) */
    IOSTransferFcn* Gpert = ios_tf_create(1, 3);
    if (Gpert) {
        Gpert->num[0] = 2.2; Gpert->den[0] = 1.0; Gpert->den[1] = 0.6; Gpert->den[2] = 0.9;
        printf("  Gap(G, Gpert):  %.6f\n", ios_gap_metric(G, Gpert));
        printf("  Nu-gap:         %.6f\n", ios_nu_gap_metric(G, Gpert));
        ios_tf_free(Gpert);
    }

    printf("\n========================================\n");
    printf("  Analysis Complete\n");
    printf("========================================\n");
    ios_tf_free(G);
}

/* Demo: Compare stability criteria on a sample system */
void ios_demo_criteria_comparison(void) {
    printf("\n=== Stability Criteria Comparison ===\n\n");

    /* Sample systems */
    struct { const char* name; double a; double b; } systems[] = {
        {"Stable 1st-order", 1.0, 1.0},
        {"Fast stable", 10.0, 1.0},
        {"Slow stable", 1.0, 0.1},
        {"Oscillatory", 0.5, 1.0},
        {"Light damped", 0.1, 1.0}
    };

    printf("%-20s | %-10s | %-10s | %-10s | %-10s\n",
           "System", "BIBO", "Passive", "Popov", "Circle");
    printf("---------------------+------------+------------+------------+------------\n");

    for (int i = 0; i < 5; i++) {
        IOSTransferFcn* tf = ios_tf_create(1, 3);
        if (!tf) continue;
        tf->num[0] = systems[i].a;
        tf->den[0] = 1.0; tf->den[1] = systems[i].b; tf->den[2] = 1.0;

        bool bibo_ok = ios_check_bibo_tf(tf) == IOS_BIBO_STABLE;
        bool pass_ok = ios_check_passivity_tf(tf).passive;
        bool popov_ok = ios_popov_criterion_tf(tf, 1.0);
        bool circle_ok = ios_circle_criterion_tf(tf, 0.0, 5.0);

        printf("%-20s | %-10s | %-10s | %-10s | %-10s\n",
               systems[i].name,
               bibo_ok ? "STABLE" : "UNSTABLE",
               pass_ok ? "YES" : "NO",
               popov_ok ? "STABLE" : "INCONCL",
               circle_ok ? "STABLE" : "INCONCL");
        ios_tf_free(tf);
    }
}

/* Demo: Gain scheduling via small-gain */
void ios_demo_gain_scheduling(void) {
    printf("\n=== Gain Scheduling via Small-Gain ===\n\n");

    printf("Scheduling gains to keep loop gain < 1:\n");
    printf("  G(s) = K/(s+1), H(s) = 1/(s+0.5)\n\n");

    for (double K = 0.5; K <= 5.0; K += 0.5) {
        IOSTransferFcn* G = ios_tf_create(1, 2);
        IOSTransferFcn* H = ios_tf_create(1, 2);
        if (!G || !H) { ios_tf_free(G); ios_tf_free(H); continue; }
        G->num[0] = K; G->den[0] = 1.0; G->den[1] = 1.0;
        H->num[0] = 1.0; H->den[0] = 1.0; H->den[1] = 0.5;

        double gG = ios_Hinf_norm_tf(G);
        double gH = ios_Hinf_norm_tf(H);
        IOSSmallGainResult sg = ios_small_gain_test(gG, gH);

        printf("  K=%.1f: ||G||=%.4f ||H||=%.4f product=%.4f margin=%.4f %s\n",
               K, gG, gH, sg.gamma, sg.margin,
               sg.satisfied ? "STABLE" : "MAY BE UNSTABLE");
        ios_tf_free(G); ios_tf_free(H);
    }
}

/* Demo: Sector-bounded nonlinearity analysis */
void ios_demo_sector_analysis(void) {
    printf("\n=== Sector-Bounded Nonlinearity Analysis ===\n\n");

    IOSTransferFcn* G = ios_tf_create(1, 2);
    if (!G) return;
    G->num[0] = 1.0; G->den[0] = 1.0; G->den[1] = 1.0;

    IOSSector sec = ios_compute_sector_tf(G);
    printf("Linear system G(s)=1/(s+1):\n");
    printf("  Sector bounds: [%.4f, %.4f]\n", sec.a, sec.b);

    /* Test different sector nonlinearities */
    typedef struct { double a; double b; const char* desc; } SectorTest;
    SectorTest tests[] = {
        {0.0, 1.0, "Saturation [0,1]"},
        {0.0, 2.0, "Deadzone+Sat [0,2]"},
        {0.5, 1.0, "Narrow [0.5,1]"},
        {0.0, 10.0, "Wide [0,10]"},
        {0.1, 0.5, "Restricted [0.1,0.5]"}
    };

    printf("\nCircle criterion results:\n");
    for (int i = 0; i < 5; i++) {
        bool ok = ios_circle_criterion_tf(G, tests[i].a, tests[i].b);
        printf("  %-25s: %s\n", tests[i].desc, ok ? "STABLE" : "INCONCLUSIVE");
    }

    ios_tf_free(G);
}

/* All comprehensive demos */
void ios_demo_all_comprehensive(void) {
    ios_demo_comprehensive_analysis();
    ios_demo_criteria_comparison();
    ios_demo_gain_scheduling();
    ios_demo_sector_analysis();
    printf("\n=== All Comprehensive Demos Complete ===\n");
}
