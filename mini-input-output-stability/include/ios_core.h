#ifndef IOS_CORE_H
#define IOS_CORE_H
#include "ios_types.h"
#include "ios_bibo.h"
#include "ios_smallgain.h"
#include "ios_passivity.h"
#define IOS_PI 3.14159265358979323846
#define IOS_EPS 1e-12
#define IOS_INF 1e100
#define IOS_DEFAULT_SIGNAL_CAP 1000
#define IOS_MAX_ITER 200
#define IOS_FREQ_POINTS 500
IOSSystem* ios_system_create(IOSStateSpace* ss);
void ios_system_free(IOSSystem* sys);
double ios_system_gain_margin(const IOSSystem* s);
void ios_system_analyze(IOSSystem* s);
void ios_print_system(const IOSSystem* s);
double ios_signal_energy(const IOSSignal* s);
double ios_signal_power(const IOSSignal* s);
double ios_signal_max(const IOSSignal* s);
double ios_signal_min(const IOSSignal* s);
double ios_signal_mean(const IOSSignal* s);
double ios_signal_rms(const IOSSignal* s);
void ios_signal_fill(IOSSignal* s, double v);
void ios_signal_ramp(IOSSignal* s, double slope);
void ios_signal_sine(IOSSignal* s, double amp, double freq, double dt);
void ios_signal_step(IOSSignal* s, double amp, int start);
void ios_signal_random(IOSSignal* s, double amp);
double ios_signal_corr(const IOSSignal* a, const IOSSignal* b, int lag);
double ios_signal_variance(const IOSSignal* s);
void ios_bode_plot(const IOSTransferFcn* tf, double* w, double* mag, double* phase, int n);
double ios_bandwidth(const IOSTransferFcn* tf);
int ios_tf_order(const IOSTransferFcn* tf);
double ios_tf_dcgain(const IOSTransferFcn* tf);
bool ios_tf_is_proper(const IOSTransferFcn* tf);
bool ios_tf_is_stable(const IOSTransferFcn* tf);
void ios_print_tf(const IOSTransferFcn* tf);
IOSSignal* ios_signal_convolve(const IOSSignal* a, const IOSSignal* b);
int ios_find_crossings(const IOSSignal* s, double thresh, int* idx, int max_n);
void ios_signal_clear(IOSSignal* s);
void ios_signal_differentiate(const IOSSignal* s, IOSSignal* ds, double dt);
void ios_signal_integrate(const IOSSignal* s, IOSSignal* is, double dt);
double ios_frequency_response(const IOSTransferFcn* tf, double w, double* mag, double* phase);
void ios_all_demos(void);
#endif
