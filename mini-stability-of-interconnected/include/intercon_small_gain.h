#ifndef INTERCON_SMALL_GAIN_H
#define INTERCON_SMALL_GAIN_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int intercon_small_gain_op0(int n, const double *x, double *out);
int intercon_small_gain_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int intercon_small_gain_op1(int n, const double *x, double *out);
int intercon_small_gain_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int intercon_small_gain_op2(int n, const double *x, double *out);
int intercon_small_gain_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int intercon_small_gain_op3(int n, const double *x, double *out);
int intercon_small_gain_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int intercon_small_gain_op4(int n, const double *x, double *out);
int intercon_small_gain_batch4(int n, int m, const double *x, double *out);
/* hyperbolic tangent element-wise transform */
int intercon_small_gain_op5(int n, const double *x, double *out);
int intercon_small_gain_batch5(int n, int m, const double *x, double *out);
/* softplus element-wise transform */
int intercon_small_gain_op6(int n, const double *x, double *out);
int intercon_small_gain_batch6(int n, int m, const double *x, double *out);
/* Gaussian RBF element-wise transform */
int intercon_small_gain_op7(int n, const double *x, double *out);
int intercon_small_gain_batch7(int n, int m, const double *x, double *out);
/* swish activation element-wise transform */
int intercon_small_gain_op8(int n, const double *x, double *out);
int intercon_small_gain_batch8(int n, int m, const double *x, double *out);
/* GELU approx element-wise transform */
int intercon_small_gain_op9(int n, const double *x, double *out);
int intercon_small_gain_batch9(int n, int m, const double *x, double *out);
/* log-cosh element-wise transform */
int intercon_small_gain_op10(int n, const double *x, double *out);
int intercon_small_gain_batch10(int n, int m, const double *x, double *out);
/* softsign element-wise transform */
int intercon_small_gain_op11(int n, const double *x, double *out);
int intercon_small_gain_batch11(int n, int m, const double *x, double *out);
double intercon_small_gain_mean(int n, const double *x);
double intercon_small_gain_variance(int n, const double *x);
int intercon_small_gain_minmax(int n, const double *x, double *out);
double intercon_small_gain_sum(int n, const double *x);
double intercon_small_gain_norm_l2(int n, const double *x);
double intercon_small_gain_norm_l1(int n, const double *x);
double intercon_small_gain_norm_linf(int n, const double *x);

#endif