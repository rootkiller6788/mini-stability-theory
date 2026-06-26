#ifndef GST_CORE_H
#define GST_CORE_H
#include <stdbool.h>
#include <stddef.h>

/* gst_core.h - Core math types for stability analysis.
 * L1: Vector, Matrix, ODEFunc, LyapunovFunc, SystemTrajectory
 * L2: state-space trajectory, Lyapunov function evaluation
 * L3: dense matrix (row-major), vector, function pointer types
 * Ref: Khalil (2002) Nonlinear Systems, LaSalle (1960) */

typedef struct { double *data; int n; } Vector;
typedef struct { double *data; int rows, cols; } Matrix;
typedef void (*ODEFunc)(const double *x, double *dxdt, int n, double t);
typedef double (*LyapunovFunc)(const double *x, int n);
typedef struct {
    double *x_history, *t_history;
    int n_steps, n_dims; double dt, t_final;
} SystemTrajectory;

Vector* vec_create(int n);
void vec_free(Vector *v);
double vec_norm(const Vector *v);
double vec_dot(const Vector *a, const Vector *b);
Vector* vec_add(const Vector *a, const Vector *b);
Vector* vec_sub(const Vector *a, const Vector *b);
Vector* vec_scale(const Vector *v, double s);
Vector* vec_copy(const Vector *v);
void vec_axpy(double a, const Vector *x, Vector *y);
void vec_to_array(const Vector *v, double *arr);
void array_to_vec(const double *arr, int n, Vector *v);

Matrix* mat_create(int rows, int cols);
void mat_free(Matrix *m);
Matrix* mat_identity(int n);
Matrix* mat_copy(const Matrix *m);
void mat_set(Matrix *m, int i, int j, double v);
double mat_get(const Matrix *m, int i, int j);
Matrix* mat_mul(const Matrix *a, const Matrix *b);
Matrix* mat_transpose(const Matrix *m);
Matrix* mat_add(const Matrix *a, const Matrix *b);
Matrix* mat_scale(const Matrix *m, double s);
double mat_trace(const Matrix *m);
double mat_det(const Matrix *m);
Matrix* mat_inverse(const Matrix *m);
double mat_frob_norm(const Matrix *m);
void mat_print(const Matrix *m, const char *name);

int mat_eigen_sym(const Matrix *A, double *eigvals);
Matrix* mat_lyapunov_solve(const Matrix *A, const Matrix *Q);

void rk4_step(ODEFunc f, double *x, int n, double t, double dt);
int rkf45_step(ODEFunc f, double *x, int n, double *t, double *dt,
               double tol, double dt_min, double dt_max);
int rk4_solve(ODEFunc f, const double* x0, int n, double t0, double tf,
              double dt, double** traj_out, int* n_steps_out);
SystemTrajectory* trajectory_simulate(ODEFunc f, const double *x0, int n,
                                       double dt, double t_final);
void trajectory_free(SystemTrajectory *traj);

/* Grid-based state-space sampling for invariant set computation */
int sample_grid(const double* bmin, const double* bmax, int n,
                int grid_n, double* samples, int max_samples);
int sample_random_uniform(const double* center, double radius, int n,
                          int nsamples, double* samples);

double lyapunov_derivative(ODEFunc f, LyapunovFunc V,
                            const double *x, int n, double eps);
void numerical_gradient(LyapunovFunc V, const double *x, int n,
                         double eps, double *grad);

bool is_hurwitz(const double *A_data, int n);
double quadratic_form(const double *P_data, const double *x, int n);
bool is_positive_definite(const double *P_data, int n);

#endif