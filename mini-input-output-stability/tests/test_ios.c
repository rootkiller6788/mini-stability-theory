#include "ios_core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

int main(void) {
    printf("=== I/O Stability Test Suite ===\n");

    /* Test 1: Signal create/free */
    { IOSSignal* s = ios_signal_create(10);
      assert(s != NULL);
      assert(s->length == 0);
      assert(s->capacity >= 10);
      ios_signal_free(s); }

    /* Test 2: Signal append and length */
    { IOSSignal* s = ios_signal_create(10);
      assert(s != NULL);
      ios_signal_append(s, 1.0);
      ios_signal_append(s, 2.0);
      assert(s->length == 2);
      assert(fabs(s->data[0] - 1.0) < 1e-10);
      assert(fabs(s->data[1] - 2.0) < 1e-10);
      ios_signal_free(s); }

    /* Test 3: Signal norm L2 */
    { IOSSignal* s = ios_signal_create(10);
      ios_signal_append(s, 3.0);
      ios_signal_append(s, 4.0);
      double n2 = ios_signal_norm(s, IOS_NORM_L2);
      assert(fabs(n2 - 5.0) < 1e-10);
      ios_signal_free(s); }

    /* Test 4: Signal norm L1 */
    { IOSSignal* s = ios_signal_create(10);
      ios_signal_append(s, 3.0);
      ios_signal_append(s, -4.0);
      double n1 = ios_signal_norm(s, IOS_NORM_L1);
      assert(fabs(n1 - 7.0) < 1e-10);
      ios_signal_free(s); }

    /* Test 5: Signal statistics */
    { IOSSignal* s = ios_signal_create(10);
      for (int i = 0; i < 5; i++) ios_signal_append(s, (double)(i + 1));
      double mn = ios_signal_mean(s);
      assert(fabs(mn - 3.0) < 1e-10);
      double mx = ios_signal_max(s);
      assert(fabs(mx - 5.0) < 1e-10);
      double mi = ios_signal_min(s);
      assert(fabs(mi - 1.0) < 1e-10);
      double energy = ios_signal_energy(s);
      assert(energy > 0.0);
      double power = ios_signal_power(s);
      assert(power > 0.0);
      double rms_val = ios_signal_rms(s);
      assert(rms_val > 0.0);
      ios_signal_free(s); }

    /* Test 6: BIBO stability - stable 2x2 system */
    { IOSStateSpace* ss = ios_ss_create(2, 1, 1);
      assert(ss != NULL);
      ss->A[0] = -1.0; ss->A[3] = -2.0;
      ss->B[0] = 1.0; ss->B[1] = 0.0;
      ss->C[0] = 1.0; ss->C[1] = 0.0;
      assert(ios_check_bibo_ss(ss) == IOS_BIBO_STABLE);
      ios_ss_free(ss); }

    /* Test 7: BIBO stability - unstable system */
    { IOSStateSpace* ss = ios_ss_create(1, 1, 1);
      ss->A[0] = 1.0;
      assert(ios_check_bibo_ss(ss) == IOS_BIBO_UNSTABLE);
      ios_ss_free(ss); }

    /* Test 8: BIBO stability - stable TF */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
      assert(ios_check_bibo_tf(tf) == IOS_BIBO_STABLE);
      ios_tf_free(tf); }

    /* Test 9: BIBO stability - unstable TF */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = -1.0;
      assert(ios_check_bibo_tf(tf) == IOS_BIBO_UNSTABLE);
      ios_tf_free(tf); }

    /* Test 10: Small-gain theorem */
    { IOSSmallGainResult r = ios_small_gain_test(0.5, 0.5);
      assert(r.satisfied);
      assert(r.gamma < 1.0); }
    { IOSSmallGainResult r = ios_small_gain_test(0.8, 0.8);
      assert(r.satisfied); }
    { IOSSmallGainResult r = ios_small_gain_test(0.9, 0.9);
      assert(r.satisfied); }
    { IOSSmallGainResult r = ios_small_gain_test(0.8, 1.5);
      assert(!r.satisfied); }

    /* Test 11: Passivity check */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
      IOSPassivityResult r = ios_check_passivity_tf(tf);
      assert(r.passive);
      ios_tf_free(tf); }

    /* Test 12: Circle criterion */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
      bool c = ios_circle_criterion_tf(tf, 0.0, 5.0);
      assert(c);
      ios_tf_free(tf); }

    /* Test 13: Disk margin */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
      double dm = ios_disk_margin(tf);
      assert(dm > 0.0);
      ios_tf_free(tf); }

    /* Test 14: Sector computation */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
      IOSSector sec = ios_compute_sector_tf(tf);
      assert(sec.type == IOS_SECTOR_FINITE);
      ios_tf_free(tf); }

    /* Test 15: H2 norm */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
      double h2 = ios_H2_norm_tf(tf);
      assert(h2 >= 0.0);
      ios_tf_free(tf); }

    /* Test 16: Hinf norm */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 1.0;
      double hi = ios_Hinf_norm_tf(tf);
      assert(hi >= 0.0);
      ios_tf_free(tf); }

    /* Test 17: Popov criterion */
    { IOSTransferFcn* tf = ios_tf_create(1, 3);
      assert(tf != NULL);
      tf->num[0] = 1.0; tf->den[0] = 1.0; tf->den[1] = 2.0; tf->den[2] = 1.0;
      assert(ios_popov_criterion_tf(tf, 1.0));
      ios_tf_free(tf); }

    /* Test 18: L2 gain from state-space */
    { IOSStateSpace* ss = ios_ss_create(2, 1, 1);
      assert(ss != NULL);
      ss->A[0] = -1.0; ss->A[3] = -2.0;
      ss->B[0] = 1.0; ss->C[0] = 1.0;
      double L2 = ios_L2_gain_ss(ss);
      assert(L2 >= 0.0);
      ios_ss_free(ss); }

    /* Test 19: Exponential stability */
    { IOSStateSpace* ss = ios_ss_create(2, 1, 1);
      assert(ss != NULL);
      ss->A[0] = -1.0; ss->A[3] = -2.0;
      ss->B[0] = 1.0; ss->C[0] = 1.0;
      assert(ios_is_exponentially_stable(ss));
      ios_ss_free(ss); }

    /* Test 20: Signal variance */
    { IOSSignal* s = ios_signal_create(10);
      ios_signal_append(s, 1.0);
      ios_signal_append(s, 2.0);
      ios_signal_append(s, 3.0);
      double v = ios_signal_variance(s);
      assert(v > 0.0);
      ios_signal_free(s); }

    /* Test 21: Signal energy and power */
    { IOSSignal* s = ios_signal_create(10);
      ios_signal_append(s, 3.0);
      ios_signal_append(s, 4.0);
      double e = ios_signal_energy(s);
      assert(fabs(e - 25.0) < 1e-10);
      double p = ios_signal_power(s);
      assert(fabs(p - 12.5) < 1e-10);
      ios_signal_free(s); }

    /* Test 22: Gap metric */
    { IOSTransferFcn* tf1 = ios_tf_create(1, 2);
      IOSTransferFcn* tf2 = ios_tf_create(1, 2);
      assert(tf1 && tf2);
      tf1->num[0] = 1.0; tf1->den[0] = 1.0; tf1->den[1] = 1.0;
      tf2->num[0] = 1.0; tf2->den[0] = 1.0; tf2->den[1] = 1.0;
      double gap = ios_gap_metric(tf1, tf2);
      /* Same systems => gap should be near 0 */
      assert(gap < 0.1);
      ios_tf_free(tf1); ios_tf_free(tf2); }

    /* Test 23: System create/analyze */
    { IOSStateSpace* ss = ios_ss_create(2, 1, 1);
      assert(ss != NULL);
      ss->A[0] = -1.0; ss->A[3] = -2.0;
      ss->B[0] = 1.0; ss->C[0] = 1.0;
      IOSSystem* sys = ios_system_create(ss);
      assert(sys != NULL);
      ios_system_analyze(sys);
      assert(sys->analyzed);
      assert(sys->bibo == IOS_BIBO_STABLE);
      ios_system_free(sys);
      ios_ss_free(ss); }

    /* Test 24: Null pointer safety */
    ios_signal_free(NULL);
    ios_ss_free(NULL);
    ios_tf_free(NULL);
    ios_system_free(NULL);
    assert(ios_check_bibo_ss(NULL) == IOS_BIBO_UNSTABLE);
    assert(ios_check_bibo_tf(NULL) == IOS_BIBO_UNSTABLE);

    /* Test 25: Signal convolution */
    { IOSSignal* a = ios_signal_create(5);
      IOSSignal* b = ios_signal_create(5);
      ios_signal_append(a, 1.0); ios_signal_append(a, 2.0);
      ios_signal_append(b, 1.0); ios_signal_append(b, 1.0);
      IOSSignal* c = ios_signal_convolve(a, b);
      assert(c != NULL);
      assert(c->length >= 3);
      ios_signal_free(a); ios_signal_free(b); ios_signal_free(c); }

    /* Test 26: Transfer function operations */
    { IOSTransferFcn* tf = ios_tf_create(1, 2);
      assert(tf != NULL);
      tf->num[0] = 2.0; tf->den[0] = 1.0; tf->den[1] = 3.0;
      assert(ios_tf_order(tf) == 1);
      assert(ios_tf_is_proper(tf));
      double dc = ios_tf_dcgain(tf);
      assert(fabs(dc - 2.0/3.0) < 1e-10);
      ios_tf_free(tf); }

    printf("All 26 tests passed.\n");
    return 0;
}
