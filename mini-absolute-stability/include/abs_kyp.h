#ifndef ABS_KYP_H
#define ABS_KYP_H

/*
 * abs_kyp.h — Kalman-Yakubovich-Popov (KYP) Lemma
 *
 * The KYP lemma bridges frequency-domain inequalities (FDI) and
 * linear matrix inequalities (LMI) in state-space.
 *
 * Central equivalence:
 *   The FDI  G(jw) + G(jw)^H >= 0  for all real w
 *   is equivalent to the existence of P = P^T > 0 such that:
 *     [ A^T*P + P*A ,  P*B - C^T ]
 *     [ B^T*P - C   ,  -(D + D^T) ]  <= 0
 *
 * For absolute stability, the KYP lemma provides a constructive
 * method to find the Lyapunov matrix P from frequency-domain
 * conditions (circle/Popov criteria).
 */

#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include "abs_core.h"

/* ── KYP Data ─────────────────────────────────────────────────────── */

/** State-space quadruple (A, B, C, D) for KYP analysis.
 *  Standard form: dx/dt = A*x + B*u, y = C*x + D*u */
typedef struct {
    double *A;     /* [n*n] state matrix */
    double *B;     /* [n*m] input matrix */
    double *C;     /* [p*n] output matrix */
    double *D;     /* [p*m] feedthrough matrix */
    int     n;     /* state dimension */
    int     m;     /* input dimension */
    int     p;     /* output dimension */
} AbsKYPStateSpace;

/** KYP problem specification: find P > 0 such that LMI(P) <= 0 */
typedef struct {
    double *M11;   /* [n*n] coefficient for A^T*P + P*A term */
    double *M12;   /* [n*m] coefficient for P*B term */
    double *M22;   /* [m*m] coefficient for D+D^T term */
    int     n;     /* state dimension */
    int     m;     /* input dimension */
} AbsKYPProblem;

/** KYP solution result */
typedef struct {
    double *P;             /* [n*n] Lyapunov matrix solution */
    double  min_eig_P;     /* minimum eigenvalue of P (must be > 0) */
    double  max_eig_P;     /* maximum eigenvalue of P */
    double  lmi_residual;  /* ||LMI(P)||_F / max(1, ||M||) */
    bool    is_feasible;   /* true if a feasible P found */
    int     iterations;    /* iterations taken by solver */
    int     n;
} AbsKYPResult;

/* ── Core KYP Functions ───────────────────────────────────────────── */

/** Solve the strict KYP lemma (positive-real lemma):
 *
 *  Find P = P^T > 0 such that:
 *    [ A^T*P + P*A ,  P*B - C^T ]  < 0
 *    [ B^T*P - C   ,  -(D + D^T) ]
 *
 *  Equivalent to: transfer function H(s) = C*(sI-A)^{-1}*B + D
 *  is strictly positive real (SPR).
 *
 *  Uses LMI approach with interior-point method (simplified).
 *
 *  @param A, B, C, D  state-space matrices
 *  @param n, m, p     dimensions
 *  @return allocated result, or NULL on failure */
AbsKYPResult* abs_kyp_solve_strict(const double *A, const double *B,
                                    const double *C, const double *D,
                                    int n, int m, int p);

/** Solve the non-strict KYP lemma:
 *  Same as above but with <= 0 (positive real, not strictly).
 *  More numerically challenging due to boundary solutions. */
AbsKYPResult* abs_kyp_solve_nonstrict(const double *A, const double *B,
                                       const double *C, const double *D,
                                       int n, int m, int p);

/** Free KYP result. */
void abs_kyp_result_free(AbsKYPResult *res);

/* ── KYP for Absolute Stability ───────────────────────────────────── */

/** KYP-based circle criterion:
 *  Construct the KYP problem from the circle condition and solve.
 *  Returns the Lyapunov matrix P confirming absolute stability.
 *  Combines the frequency-domain test with state-space verification. */
AbsKYPResult* abs_kyp_circle_to_lmi(const AbsLureSystem *sys);

/** KYP-based Popov criterion:
 *  For time-invariant phi in [0, k], constructs the multiplier
 *  M(s) = 1 + q*s and solves the KYP LMI.
 *  Returns P and q if feasible. */
AbsKYPResult* abs_kyp_popov_to_lmi(const AbsLureSystem *sys, double q);

/** KYP-based verification: given P, check if it satisfies the LMI
 *  for the given Lur'e system and sector. */
bool abs_kyp_verify_solution(const AbsLureSystem *sys,
                              const double *P, double tol);

/* ── Positive Real Lemma ──────────────────────────────────────────── */

/** Check if a transfer function H(s) = C*(sI-A)^{-1}*B + D
 *  is strictly positive real (SPR).
 *  Returns true if both:
 *    1. A is Hurwitz
 *    2. The KYP LMI is feasible with P > 0 */
bool abs_kyp_is_spr(const double *A, const double *B,
                     const double *C, const double *D,
                     int n, int m, int p);

/** Check if H(s) is positive real (PR, non-strict).
 *  Returns true if H(jw) + H(jw)^H >= 0 for all w. */
bool abs_kyp_is_pr(const double *A, const double *B,
                    const double *C, const double *D,
                    int n, int m, int p);

/** Bounded-real lemma: check if ||H(s)||_inf < gamma.
 *  Equivalent to the existence of P > 0 such that:
 *    [ A^T*P + P*A + C^T*C ,  P*B + C^T*D           ]
 *    [ B^T*P + D^T*C       ,  D^T*D - gamma^2*I    ] < 0 */
bool abs_kyp_bounded_real(const double *A, const double *B,
                           const double *C, const double *D,
                           int n, int m, int p, double gamma);

/* ── LMI Utilities ────────────────────────────────────────────────── */

/** Solve a simple Lyapunov-type LMI:
 *    A^T*P + P*A + Q < 0
 *  by solving the Lyapunov equation A^T*P + P*A = -(Q + eps*I)
 *  and checking positive definiteness.
 *  Returns P in row-major format (caller allocates n*n doubles). */
bool abs_kyp_solve_lyap_lmi(const double *A, const double *Q,
                             int n, double *P);

/** Form the KYP LMI matrix for a given P:
 *    LMI(P) = [ A^T*P + P*A,  P*B - C^T ]
 *             [ B^T*P - C,    -(D+D^T)  ]
 *  Stored row-major in lmi[(n+m)*(n+m)].
 *  Returns the maximum eigenvalue (<= 0 means feasible). */
double abs_kyp_form_lmi(const double *A, const double *B,
                         const double *C, const double *D,
                         int n, int m,
                         const double *P, double *lmi);

/** Check if an (n+m)×(n+m) symmetric matrix is negative semidefinite.
 *  Uses Cholesky on -M and returns true if successful. */
bool abs_kyp_is_negative_semidef(int dim, const double *M, double tol);

#endif /* ABS_KYP_H */
