#ifndef FTS_ANALYSIS_H
#define FTS_ANALYSIS_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int fts_analysis_op0(int n, const double *x, double *out);
int fts_analysis_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int fts_analysis_op1(int n, const double *x, double *out);
int fts_analysis_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int fts_analysis_op2(int n, const double *x, double *out);
int fts_analysis_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int fts_analysis_op3(int n, const double *x, double *out);
int fts_analysis_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int fts_analysis_op4(int n, const double *x, double *out);
int fts_analysis_batch4(int n, int m, const double *x, double *out);
/* hyperbolic tangent element-wise transform */
int fts_analysis_op5(int n, const double *x, double *out);
int fts_analysis_batch5(int n, int m, const double *x, double *out);
/* softplus element-wise transform */
int fts_analysis_op6(int n, const double *x, double *out);
int fts_analysis_batch6(int n, int m, const double *x, double *out);
/* Gaussian RBF element-wise transform */
int fts_analysis_op7(int n, const double *x, double *out);
int fts_analysis_batch7(int n, int m, const double *x, double *out);
/* swish activation element-wise transform */
int fts_analysis_op8(int n, const double *x, double *out);
int fts_analysis_batch8(int n, int m, const double *x, double *out);
double fts_analysis_mean(int n, const double *x);
double fts_analysis_variance(int n, const double *x);
int fts_analysis_minmax(int n, const double *x, double *out);
double fts_analysis_sum(int n, const double *x);
double fts_analysis_norm_l2(int n, const double *x);
double fts_analysis_norm_l1(int n, const double *x);
double fts_analysis_norm_linf(int n, const double *x);

#endif