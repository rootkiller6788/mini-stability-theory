#ifndef IOS_BIBO_H
#define IOS_BIBO_H
#include "ios_types.h"
IOSSignal* ios_signal_create(int cap);
void ios_signal_free(IOSSignal* s);
void ios_signal_append(IOSSignal* s, double v);
double ios_signal_norm(const IOSSignal* s, IOSNormType type);
double ios_impulse_response(const IOSTransferFcn* tf, double t);
double ios_step_response(const IOSTransferFcn* tf, double t);
double ios_convolution(const IOSSignal* u, const IOSTransferFcn* tf, double t);

IOSBIBOStatus ios_check_bibo_tf(const IOSTransferFcn* tf);
IOSBIBOStatus ios_check_bibo_ss(const IOSStateSpace* sys);
bool ios_is_bibo_stable(const IOSSystem* sys);
double ios_bibo_gain(const IOSSystem* sys, IOSNormType type);

IOSStateSpace* ios_ss_create(int n, int m, int p);
void ios_ss_free(IOSStateSpace* sys);
void ios_ss_simulate(const IOSStateSpace* sys, const IOSSignal* u, IOSSignal* y);
IOSTransferFcn* ios_tf_create(int num_n, int den_n);
void ios_tf_free(IOSTransferFcn* tf);
double ios_tf_eval(const IOSTransferFcn* tf, double s);
double ios_Linf_induced_norm(const IOSStateSpace* sys);
double ios_L1_induced_norm(const IOSStateSpace* sys);
double ios_peak_gain(const IOSStateSpace* sys);
bool ios_is_exponentially_stable(const IOSStateSpace* sys);
double ios_settling_time(const IOSStateSpace* sys, double tol);
double ios_overshoot_estimate(const IOSStateSpace* sys);
double ios_rise_time_estimate(const IOSStateSpace* sys);
void ios_bibo_margins(const IOSTransferFcn* tf, double* gm, double* pm);
void ios_impulse_response_full(const IOSStateSpace* sys, double t_max, double dt, IOSSignal* h);
void ios_step_response_full(const IOSStateSpace* sys, double t_max, double dt, IOSSignal* y);
double ios_impulse_energy(const IOSStateSpace* sys);
double ios_L1_norm_impulse(const IOSTransferFcn* tf);
double ios_H2_norm_tf(const IOSTransferFcn* tf);
double ios_Hinf_norm_tf(const IOSTransferFcn* tf);
double ios_h2_norm_ss(const IOSStateSpace* sys);
double ios_impulse_response_tf(const IOSTransferFcn* tf, double t);
double ios_step_response_tf(const IOSTransferFcn* tf, double t);
double ios_impulse_response_ss(const IOSStateSpace* sys, double t);
double ios_ss_eigen_max_real(const IOSStateSpace* sys);
void ios_print_gain_table(const IOSStateSpace* sys);
double ios_dc_gain_ss(const IOSStateSpace* sys);
double ios_system_norm_ratio(const IOSStateSpace* sys);
int ios_bibo_margin_check(const IOSStateSpace* sys, double margin);
void ios_bibo_test_suite(void);
#endif
