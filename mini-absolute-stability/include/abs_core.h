#ifndef ABS_CORE_H
#define ABS_CORE_H

/*
 * abs_core.h — Core data structures and linear algebra for absolute stability theory
 *
 * Absolute stability studies Lur'e systems:
 *   dx/dt = A*x + b*u     (linear forward path)
 *   y     = c^T * x        (output)
 *   u     = -phi(t, y)     (nonlinear feedback)
 *
 * where phi belongs to sector [alpha, beta]:
 *   alpha*y^2 <= y*phi(y) <= beta*y^2   for all y
 *
 * Transfer function:  G(s) = c^T * (sI - A)^{-1} * b
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ── Mathematical Constants ───────────────────────────────────────── */

#define ABS_PI        3.14159265358979323846   /* pi */
#define ABS_EPS       1e-12                    /* default tolerance */
#define ABS_EPS_SQRT  1e-6                     /* sqrt tolerance */
#define ABS_MAX_DIM   16                       /* max system dimension */
#define ABS_MAX_ITER  1000                     /* max iterations for iterative solvers */
#define ABS_TWO_PI    (2.0 * ABS_PI)           /* 2*pi */
#define ABS_HALF_PI   (ABS_PI / 2.0)           /* pi/2 */

/* ── Sector Nonlinearity Types ────────────────────────────────────── */

/** Nonlinearity class: sector, saturation, deadzone, hysteresis, etc. */
typedef enum {
    ABS_PHI_GENERIC      = 0,  /* generic sector-bounded */
    ABS_PHI_SATURATION   = 1,  /* sat(u) = clip(u, -L, L) */
    ABS_PHI_DEADZONE     = 2,  /* deadzone of width delta */
    ABS_PHI_RELAY        = 3,  /* ideal relay: sgn(u) */
    ABS_PHI_HYSTERESIS   = 4,  /* hysteresis loop */
    ABS_PHI_CUBIC        = 5,  /* phi(y) = y^3 */
    ABS_PHI_SINUSOIDAL   = 6,  /* phi(y) = sin(y) */
    ABS_PHI_TANH         = 7   /* phi(y) = tanh(y) */
} AbsPhiType;

/* ── Core Data Structures ─────────────────────────────────────────── */

/**
 * AbsLureSystem — Lur'e problem description
 *
 * Linear part:  dx/dt = A*x + b*u,  y = c^T*x
 * Nonlinearity: u = -phi(y),  phi in sector [alpha, beta]
 *
 * Memory layout: A is stored row-major as A[n*n] array.
 * b and c are vectors of length n.
 */
typedef struct {
    double *A;        /* [n*n] state matrix, row-major */
    double *b;        /* [n]   input vector */
    double *c;        /* [n]   output vector */
    double alpha;     /* lower sector bound */
    double beta;      /* upper sector bound */
    int    n;         /* system dimension */
} AbsLureSystem;

/**
 * AbsSector — sector condition parameters
 *
 * A nonlinearity phi(y) belongs to sector [alpha, beta] if:
 *   (phi(y) - alpha*y) * (beta*y - phi(y)) >= 0   for all y
 */
typedef struct {
    double alpha;     /* lower slope bound */
    double beta;      /* upper slope bound */
    double center;    /* (alpha+beta)/2, center of sector */
    double radius;    /* (beta-alpha)/2, half-width of sector */
} AbsSector;

/**
 * AbsComplex — complex number for frequency-domain analysis
 */
typedef struct {
    double re;        /* real part */
    double im;        /* imaginary part */
} AbsComplex;

/**
 * AbsFrechet — Frechet derivative structure for sensitivity analysis
 * of Lyapunov solutions with respect to A matrix perturbations.
 */
typedef struct {
    double *dP_dA;    /* [n*n * n*n] Frechet derivative matrix (row-major) */
    double norm;      /* Frobenius norm of derivative */
    int    n;         /* dimension */
} AbsFrechet;

/* ── Lur'e System Lifecycle ──────────────────────────────────────── */

/** Allocate and initialize a Lur'e system of dimension n.
 *  A, b, c are copied into internally allocated storage.
 *  Returns NULL on allocation failure or invalid parameters. */
AbsLureSystem* abs_lure_create(const double *A, const double *b,
                                const double *c, int n,
                                double alpha, double beta);

/** Deallocate a Lur'e system and all its internal arrays. */
void abs_lure_free(AbsLureSystem *sys);

/** Create a copy of a Lur'e system. Returns NULL on failure. */
AbsLureSystem* abs_lure_clone(const AbsLureSystem *sys);

/** Validate Lur'e system: check A Hurwitz, sector valid, dimensions. */
bool abs_lure_validate(const AbsLureSystem *sys);

/* ── Sector Operations ────────────────────────────────────────────── */

/** Initialize a sector with given bounds. Swaps if alpha > beta. */
AbsSector abs_sector_make(double alpha, double beta);

/** Check if value y and phi(y) satisfy sector condition. */
bool abs_sector_check(AbsSector sec, double y, double phi_y);

/** Compute the Nyquist disk center for circle criterion:
 *  center = -(alpha+beta) / (2*alpha*beta)         (for 0 < alpha < beta)
 *  center = -1/beta                                (for alpha = 0)
 *  Returns true if well-defined, false if degenerate. */
bool abs_sector_disk_center(AbsSector sec, double *center_real);

/** Compute Nyquist disk radius for circle criterion:
 *  radius = (beta-alpha) / (2*|alpha*beta|)        (for alpha != 0)
 *  radius = 1/(2*beta)                             (for alpha = 0) */
double abs_sector_disk_radius(AbsSector sec);

/** Compute the Hurwitz sector [k_min, k_max] for a given A matrix:
 *  All k in (k_min, k_max) make A-k*b*c^T Hurwitz stable. */
bool abs_sector_hurwitz_range(const double *A, const double *b,
                               const double *c, int n,
                               double *k_min, double *k_max);

/* ── Linear Algebra Utilities ─────────────────────────────────────── */

/** Matrix-vector multiply: y = alpha*A*x + beta*y  (A is n×n row-major) */
void abs_linalg_gemv(int n, double alpha_coeff, const double *A,
                     const double *x, double beta_coeff, double *y);

/** Matrix-matrix multiply: C = alpha*A*B + beta*C  (all n×n row-major) */
void abs_linalg_gemm(int n, double alpha_coeff, const double *A,
                     const double *B, double beta_coeff, double *C);

/** Vector 2-norm */
double abs_linalg_norm2(int n, const double *x);

/** Dot product */
double abs_linalg_dot(int n, const double *a, const double *b);

/** Matrix trace */
double abs_linalg_trace(int n, const double *A);

/** Matrix transpose in-place: B = A^T (A and B must not alias) */
void abs_linalg_transpose(int n, const double *A, double *B);

/** Check if matrix is symmetric: max|A - A^T| < tol */
bool abs_linalg_is_symmetric(int n, const double *A, double tol);

/** Check if matrix is positive definite via Cholesky attempt.
 *  Returns true if Cholesky factorization succeeds. */
bool abs_linalg_is_positive_definite(int n, const double *A);

/** Cholesky decomposition: A = L*L^T, lower triangular L stored in A.
 *  Returns 0 on success, -1 if A is not positive definite. */
int abs_linalg_cholesky(int n, double *A);

/** Solve L*x = b  (L is n×n lower triangular, row-major).
 *  x and b may alias. */
void abs_linalg_forward_sub(int n, const double *L, const double *b, double *x);

/** Solve U*x = b  (U is n×n upper triangular, row-major).
 *  x and b may alias. */
void abs_linalg_backward_sub(int n, const double *U, const double *b, double *x);

/** Copy n×n matrix: dst = src */
void abs_linalg_copy_matrix(int n, const double *src, double *dst);

/** Set n×n matrix to identity */
void abs_linalg_identity(int n, double *A);

/** Scale n×n matrix: A = alpha * A */
void abs_linalg_scale(int n, double alpha, double *A);

/** Compute eigenvalues of real symmetric n×n matrix via Jacobi method.
 *  A is modified (destroyed). eigenvalues[n] receives result.
 *  Returns number of iterations used, or -1 on failure. */
int abs_linalg_symeig(int n, double *A, double *eigenvalues);

/** Compute the transfer function G(s) = c^T*(sI-A)^{-1}*b
 *  at complex frequency s. Result stored in G. */
bool abs_linalg_transfer(const AbsLureSystem *sys, AbsComplex s, AbsComplex *G);

/** Frequency response G(jw) = c^T*(jw*I - A)^{-1}*b */
bool abs_linalg_freqresp(const AbsLureSystem *sys, double omega,
                          AbsComplex *G);

#endif /* ABS_CORE_H */
