#ifndef POPOV_NYQUIST_H
#define POPOV_NYQUIST_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int popov_nyquist_op0(int n, const double *x, double *out);
int popov_nyquist_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int popov_nyquist_op1(int n, const double *x, double *out);
int popov_nyquist_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int popov_nyquist_op2(int n, const double *x, double *out);
int popov_nyquist_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int popov_nyquist_op3(int n, const double *x, double *out);
int popov_nyquist_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int popov_nyquist_op4(int n, const double *x, double *out);
int popov_nyquist_batch4(int n, int m, const double *x, double *out);
/* hyperbolic tangent element-wise transform */
int popov_nyquist_op5(int n, const double *x, double *out);
int popov_nyquist_batch5(int n, int m, const double *x, double *out);
/* softplus element-wise transform */
int popov_nyquist_op6(int n, const double *x, double *out);
int popov_nyquist_batch6(int n, int m, const double *x, double *out);
/* Gaussian RBF element-wise transform */
int popov_nyquist_op7(int n, const double *x, double *out);
int popov_nyquist_batch7(int n, int m, const double *x, double *out);
double popov_nyquist_mean(int n, const double *x);
double popov_nyquist_variance(int n, const double *x);
int popov_nyquist_minmax(int n, const double *x, double *out);
double popov_nyquist_sum(int n, const double *x);
double popov_nyquist_norm_l2(int n, const double *x);
double popov_nyquist_norm_l1(int n, const double *x);
double popov_nyquist_norm_linf(int n, const double *x);

#endif