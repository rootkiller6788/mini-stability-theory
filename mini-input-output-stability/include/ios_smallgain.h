#ifndef IOS_SMALLGAIN_H
#define IOS_SMALLGAIN_H
#include "ios_types.h"
IOSSmallGainResult ios_small_gain_test(double gamma1, double gamma2);
bool ios_small_gain_stable(const IOSSystem* sys1, const IOSSystem* sys2);
double ios_L2_gain_ss(const IOSStateSpace* sys);
double ios_L2_gain_tf(const IOSTransferFcn* tf);

bool ios_circle_criterion(const IOSStateSpace* sys, double k1, double k2);
bool ios_circle_criterion_tf(const IOSTransferFcn* tf, double k1, double k2);
double ios_nyquist_distance(const IOSTransferFcn* tf, double w);
double ios_disk_margin(const IOSTransferFcn* tf);

typedef struct { double center; double radius; bool stable; } IOSCircleResult;
IOSCircleResult ios_circle_analysis(const IOSTransferFcn* tf, double k1, double k2);
double ios_sensitivity_peak(const IOSTransferFcn* L);
double ios_complementary_sensitivity(const IOSTransferFcn* L);
double ios_gap_metric(const IOSTransferFcn* G1, const IOSTransferFcn* G2);
double ios_structured_singular_value(const IOSStateSpace* sys);
bool ios_zames_falb_criterion(const IOSTransferFcn* tf, double k);
void ios_sg_report(double g1, double g2);
void ios_circle_report(const IOSTransferFcn* tf, double k1, double k2);
double ios_sensitivity_peak_tf(const IOSTransferFcn* L);
double ios_complementary_sensitivity_tf(const IOSTransferFcn* L);
void ios_smallgain_report_full(double g1, double g2);
double ios_nichols_peak(const IOSTransferFcn* tf);
double ios_nu_gap_metric(const IOSTransferFcn* G1, const IOSTransferFcn* G2);
void ios_circle_comprehensive_report(const IOSTransferFcn* G, double k1, double k2);
void ios_smallgain_comprehensive_report(const IOSTransferFcn* G1, const IOSTransferFcn* G2);
#endif
