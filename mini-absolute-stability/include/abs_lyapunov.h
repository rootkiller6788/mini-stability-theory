#ifndef ABS_LYAPUNOV_H
#define ABS_LYAPUNOV_H

/*
 * abs_lyapunov.h — Lyapunov equation solvers for absolute stability
 *
 * The continuous-time Lyapunov equation:
 *   A^T * P + P * A = -Q
 *
 * If A is Hurwitz (all eigenvalues have negative real parts) and
 * Q = Q^T > 0, then P = P^T > 0 is the unique solution.
 *
 * V(x) = x^T * P * x is a Lyapunov function proving stability.
 *
 * For absolute stability, we solve:
 *   (A - alpha*b*c^T)^T * P + P * (A - alpha*b*c^T) = -Q
 * This P is used to construct the Lur'e-Postnikov Lyapunov function.
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>

/* ── Lyapunov Solution Method ────────────────────────────────────── */

typedef enum {
    ABS_LYAP_BARTELS_STEWART = 0,  /* Bartels-Stewart: real Schur + Sylvester */
    ABS_LYAP_DIRECT_KRONECKER = 1, /* Kronecker product direct solve */
    ABS_LYAP_MATRIX_SIGN     = 2,  /* Matrix sign function iteration */
    ABS_LYAP_ALTERNATING     = 3   /* Alternating direction implicit (ADI) */
} AbsLyapMethod;

/* ── Lyapunov Solution Result ────────────────────────────────────── */

typedef struct {
    double *P;            /* [n*n] solution matrix, row-major */
    int     n;            /* dimension */
    double  residual;     /* ||A^T*P + P*A + Q||_F / ||Q||_F */
    int     iterations;   /* iterations used (for iterative methods) */
    double  cond_P;       /* condition number estimate of P */
    double  min_eig_P;    /* minimum eigenvalue of P */
    double  max_eig_P;    /* maximum eigenvalue of P */
} AbsLyapResult;

/* ── Functions ───────────────────────────────────────────────────── */

/** Solve continuous-time Lyapunov equation: A^T*P + P*A = -Q
 *
 * @param A      [n*n] system matrix (must be Hurwitz for unique solution)
 * @param Q      [n*n] right-hand side (typically positive definite)
 * @param n      dimension
 * @param method solution algorithm to use
 * @return       allocated result with P, or NULL on failure.
 *               Caller must abs_lyap_result_free(). */
AbsLyapResult* abs_lyap_solve(const double *A, const double *Q, int n,
                               AbsLyapMethod method);

/** Free a Lyapunov result structure. */
void abs_lyap_result_free(AbsLyapResult *res);

/** Bartels-Stewart algorithm (most numerically stable).
 *  Uses real Schur decomposition: A = U*T*U^T, then solves
 *  the transformed equation T^T*Y + Y*T = -U^T*Q*U via
 *  forward substitution on the quasi-triangular T.
 *
 *  Complexity: O(n^3). Preferred for n <= ABS_MAX_DIM. */
AbsLyapResult* abs_lyap_bartels_stewart(const double *A, const double *Q, int n);

/** Direct Kronecker product solution.
 *  Forms the n^2 × n^2 system:
 *    (I ⊗ A^T + A^T ⊗ I) * vec(P) = -vec(Q)
 *  and solves via Gaussian elimination.
 *  Complexity: O(n^6). Only practical for small n. */
AbsLyapResult* abs_lyap_kronecker(const double *A, const double *Q, int n);

/** Matrix sign function iteration.
 *  Uses: P = (1/2) * sign([A^T, Q; 0, -A])_12
 *  Iterative, converges quadratically.
 *  Complexity: O(n^3 * iter). */
AbsLyapResult* abs_lyap_matrix_sign(const double *A, const double *Q, int n);

/** Compute a Lyapunov function P for the Lur'e system linear part.
 *  Solves A^T*P + P*A = -I and returns P. */
AbsLyapResult* abs_lyap_for_lure(const double *A, int n);

/** Compute the sector-shifted Lyapunov solution:
 *  (A - alpha*b*c^T)^T * P + P * (A - alpha*b*c^T) = -Q
 *  Used to construct the Lur'e Lyapunov function:
 *  V(x) = x^T*P*x + integral_0^{c^T*x} phi(sigma) dsigma */
AbsLyapResult* abs_lyap_sector_shifted(const double *A, const double *b,
                                        const double *c, int n,
                                        double alpha, const double *Q);

/** Check the Lyapunov inequality condition for absolute stability:
 *  Given P = P^T > 0, verify:
 *    x^T * (A^T*P + P*A) * x + 2*x^T*P*b * (-phi(c^T*x)) < 0
 *  for a set of test points in a ball around the origin.
 *
 *  Returns the minimum negative definiteness margin found. */
double abs_lyap_verify_stability(const double *A, const double *b,
                                  const double *c, int n,
                                  const double *P, double alpha, double beta,
                                  int num_test_points);

/** Compute the condition number of a Lyapunov solution.
 *  cond = 2*||A||_F / sep(A^T, -A)
 *  where sep is the separation of A^T and -A. */
double abs_lyap_condition(const double *A, int n);

/** Estimate the separation sep(A^T, -A) = min_{X≠0} ||A^T*X + X*A||_F / ||X||_F
 *  via power iteration on the Sylvester operator. */
double abs_lyap_separation(const double *A, int n);

/** Compute the real Schur decomposition of a matrix.
 *  A = Q*T*Q^T where Q is orthogonal and T is quasi-upper-triangular
 *  (1×1 or 2×2 blocks on diagonal for real eigenvalues or complex pairs).
 *  Input A is destroyed. Q and T are output.
 *  Returns 0 on success, -1 on failure. */
int abs_lyap_schur_decompose(int n, double *A, double *Q, double *T);

#endif /* ABS_LYAPUNOV_H */
