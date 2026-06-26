#ifndef IOS_TYPES_H
#define IOS_TYPES_H
#include <stdbool.h>
#include <stddef.h>
/* Input-Output Stability Theory -- Core Types
 * BIBO, L2, Lp stability, small-gain, passivity, circle criterion.
 * Ref: Desoer & Vidyasagar (1975), Khalil (2002), Zames (1966) */

typedef struct { double* data; int length; int capacity; } IOSSignal;
typedef struct { double* num; double* den; int n_num; int n_den; } IOSTransferFcn;
typedef struct { double* A; double* B; double* C; double* D; int n; int m; int p; } IOSStateSpace;

typedef enum { IOS_NORM_L1=0, IOS_NORM_L2=1, IOS_NORM_LINF=2 } IOSNormType;
typedef enum { IOS_BIBO_STABLE=0, IOS_BIBO_UNSTABLE=1, IOS_BIBO_CRITICAL=2 } IOSBIBOStatus;
typedef enum { IOS_SECTOR_NONE=0, IOS_SECTOR_FINITE=1, IOS_SECTOR_INFINITE=2, IOS_SECTOR_INCREMENTAL=3 } IOSSectorType;

typedef struct { double a; double b; IOSSectorType type; } IOSSector;
typedef struct { double gamma; bool satisfied; double margin; } IOSSmallGainResult;
typedef struct { bool passive; bool strictly_input; bool strictly_output; double dissipation; } IOSPassivityResult;
typedef struct { double gain_margin; double phase_margin; double delay_margin; bool stable; } IOSMarginResult;

typedef double (*IOSMapFunction)(const double* u, double* y, int n, void* params);

typedef struct {
    IOSStateSpace* sys; IOSMapFunction map; void* params;
    double L2_gain; double L1_gain; double Linf_gain;
    IOSBIBOStatus bibo; IOSSector sector;
    IOSSmallGainResult sg; IOSPassivityResult passivity;
    IOSMarginResult margins; bool analyzed;
} IOSSystem;
#endif
