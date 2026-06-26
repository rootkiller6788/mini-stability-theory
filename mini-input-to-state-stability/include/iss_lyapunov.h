#ifndef ISS_LYAPUNOV_H
#define ISS_LYAPUNOV_H
#include <stddef.h>
#include <math.h>

/* absolute value element-wise transform */
int iss_lyapunov_op0(int n, const double *x, double *out);
int iss_lyapunov_batch0(int n, int m, const double *x, double *out);
/* square element-wise transform */
int iss_lyapunov_op1(int n, const double *x, double *out);
int iss_lyapunov_batch1(int n, int m, const double *x, double *out);
/* cube element-wise transform */
int iss_lyapunov_op2(int n, const double *x, double *out);
int iss_lyapunov_batch2(int n, int m, const double *x, double *out);
/* logistic sigmoid element-wise transform */
int iss_lyapunov_op3(int n, const double *x, double *out);
int iss_lyapunov_batch3(int n, int m, const double *x, double *out);
/* ReLU activation element-wise transform */
int iss_lyapunov_op4(int n, const double *x, double *out);
int iss_lyapunov_batch4(int n, int m, const double *x, double *out);
/* hyperbolic tangent element-wise transform */
int iss_lyapunov_op5(int n, const double *x, double *out);
int iss_lyapunov_batch5(int n, int m, const double *x, double *out);
/* softplus element-wise transform */
int iss_lyapunov_op6(int n, const double *x, double *out);
int iss_lyapunov_batch6(int n, int m, const double *x, double *out);
/* Gaussian RBF element-wise transform */
int iss_lyapunov_op7(int n, const double *x, double *out);
int iss_lyapunov_batch7(int n, int m, const double *x, double *out);
/* swish activation element-wise transform */
int iss_lyapunov_op8(int n, const double *x, double *out);
int iss_lyapunov_batch8(int n, int m, const double *x, double *out);
/* GELU approx element-wise transform */
int iss_lyapunov_op9(int n, const double *x, double *out);
int iss_lyapunov_batch9(int n, int m, const double *x, double *out);
/* log-cosh element-wise transform */
int iss_lyapunov_op10(int n, const double *x, double *out);
int iss_lyapunov_batch10(int n, int m, const double *x, double *out);
/* softsign element-wise transform */
int iss_lyapunov_op11(int n, const double *x, double *out);
int iss_lyapunov_batch11(int n, int m, const double *x, double *out);
double iss_lyapunov_mean(int n, const double *x);
double iss_lyapunov_variance(int n, const double *x);
int iss_lyapunov_minmax(int n, const double *x, double *out);
double iss_lyapunov_sum(int n, const double *x);
double iss_lyapunov_norm_l2(int n, const double *x);
double iss_lyapunov_norm_l1(int n, const double *x);
double iss_lyapunov_norm_linf(int n, const double *x);

#endif