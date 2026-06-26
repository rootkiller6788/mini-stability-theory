#ifndef LYAP_CONSTRUCTION_H
#define LYAP_CONSTRUCTION_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int lyap_construction_op0(int n, const double *x, double *out);
int lyap_construction_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int lyap_construction_op1(int n, const double *x, double *out);
int lyap_construction_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int lyap_construction_op2(int n, const double *x, double *out);
int lyap_construction_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int lyap_construction_op3(int n, const double *x, double *out);
int lyap_construction_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int lyap_construction_op4(int n, const double *x, double *out);
int lyap_construction_batch4(int n, int m, const double *x, double *out);
double lyap_construction_mean(int n, const double *x);
double lyap_construction_variance(int n, const double *x);
int lyap_construction_minmax(int n, const double *x, double *out);
double lyap_construction_sum(int n, const double *x);
double lyap_construction_norm_l2(int n, const double *x);
double lyap_construction_norm_l1(int n, const double *x);
double lyap_construction_norm_linf(int n, const double *x);

#endif