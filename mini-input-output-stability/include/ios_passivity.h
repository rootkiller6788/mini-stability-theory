#ifndef IOS_PASSIVITY_H
#define IOS_PASSIVITY_H
#include "ios_types.h"
IOSPassivityResult ios_check_passivity_ss(const IOSStateSpace* sys);
IOSPassivityResult ios_check_passivity_tf(const IOSTransferFcn* tf);
bool ios_is_passive(const IOSSystem* sys);
bool ios_is_strictly_passive(const IOSSystem* sys, double epsilon);

bool ios_popov_criterion(const IOSStateSpace* sys, double k);
bool ios_popov_criterion_tf(const IOSTransferFcn* tf, double k);
double ios_popov_sector_bound(const IOSStateSpace* sys);

double ios_dissipation_inequality(const IOSStateSpace* sys, const IOSSignal* u, const IOSSignal* y);
double ios_storage_function(const IOSStateSpace* sys, const double* x);
bool ios_is_dissipative(const IOSStateSpace* sys, double supply_rate);

IOSSector ios_compute_sector_ss(const IOSStateSpace* sys);
IOSSector ios_compute_sector_tf(const IOSTransferFcn* tf);
bool ios_in_conic_sector(const IOSTransferFcn* tf, double a, double b);
double ios_passivity_index(const IOSTransferFcn* tf);
bool ios_strict_positive_real(const IOSTransferFcn* tf, double eps);
double ios_positive_real_condition(const IOSTransferFcn* tf);
bool ios_system_is_passive(const IOSStateSpace* sys);
void ios_passivity_report(const IOSStateSpace* sys);
void ios_sector_report(const IOSTransferFcn* tf);
void ios_passivity_all_metrics(const IOSStateSpace* sys);
double ios_excess_passivity(const IOSStateSpace* sys);
double ios_shortage_passivity(const IOSStateSpace* sys);
bool ios_feedback_passivity_stable(const IOSStateSpace* G, const IOSStateSpace* H);
double ios_dissipation_rate(const IOSStateSpace* sys, const double* x, const double* u);
#endif
