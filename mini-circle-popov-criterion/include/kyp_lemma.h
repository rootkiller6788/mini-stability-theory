#ifndef KYP_LEMMA_H
#define KYP_LEMMA_H

/*
 * kyp_lemma.h — Kalman-Yakubovich-Popov Lemma
 *
 * The KYP lemma establishes the fundamental equivalence between
 * frequency-domain inequalities (FDIs) and linear matrix inequalities (LMIs)
 * in state-space form. It is the bridge between classical control (Nyquist,
 * Popov, circle criterion) and modern control (Lyapunov, LMI, H-infinity).
 *
 * Statement (Positive Real Lemma):
 *   Given a minimal state-space realization (A, B, C, D) of G(s), the
 *   following are equivalent:
 *
 *   (FDI)  G(jw) + G(-jw)^T + 2*epsilon*I >= 0  for all w
 *
 *   (LMI)  There exists P = P^T > 0 such that:
 *          [ A^T P + P A    P B - C^T ]
 *          [ B^T P - C      -(D + D^T) ]  <=  0
 *
 *   With the strict version (KYP proper):
 *          A^T P + P A = -L L^T - epsilon*P
 *          P B - C^T   = -L W
 *          D + D^T     =  W^T W
 *
 * The KYP lemma provides:
 *   1. Frequency-domain -> state-space: find P from G(jw) condition.
 *   2. State-space -> frequency-domain: verify G(jw) from LMI solution.
 *   3. Passivity verification: G(s) is positive real iff KYP holds.
 *
 * References:
 *   - Kalman, "Lyapunov functions for the problem of Lur'e", 1963
 *   - Yakubovich, "Solution of certain matrix inequalities", 1962
 *   - Popov, "Hyperstability of Control Systems", 1973
 *   - Boyd et al., "Linear Matrix Inequalities in System and Control", 1994
 */

#include "popov_criterion.h"

/* ---- State-Space Representation ---- */
/* Continuous-time LTI system:
 *   dx/dt = A*x + B*u
 *   y     = C*x + D*u
 *
 * A is stored row-major: A(i,j) = A[i*n + j]
 * B is a column vector of length n
 * C is a row vector of length n
 * D is a scalar (SISO systems)
 */
typedef struct {
    double *A;      /* n x n state matrix (row-major)     */
    double *B;      /* n x 1 input vector                  */
    double *C;      /* 1 x n output vector                 */
    double  D;      /* feedthrough scalar                  */
    int     n;      /* state dimension                     */
} StateSpace;

/* ---- KYP Solution ---- */
typedef struct {
    double *P;          /* n x n Lyapunov matrix (row-major)  */
    double *L;          /* n x 1 factor vector                */
    double  W;          /* feedthrough factor                 */
    int     n;          /* state dimension                    */
    bool    is_feasible; /* true if KYP LMI is feasible        */
    int     iterations;  /* number of iterations used          */
    double  residual;    /* final LMI residual norm            */
} KYPSolution;

/* ---- State-Space Lifecycle ---- */

/* Create a state-space system of dimension n.
 * All matrices are zero-initialized. */
StateSpace* ss_create(int n);

/* Free a state-space system. */
void ss_free(StateSpace* ss);

/* Set the state-space matrices from arrays.
 * A is n x n, B is n x 1, C is 1 x n, D is scalar. */
void ss_set_matrices(StateSpace* ss,
                     const double *A, const double *B,
                     const double *C, double D);

/* ---- KYP Solver ---- */

/* Solve the KYP LMI for a given state-space system with tolerance epsilon.
 *
 * Approach:
 *   For small systems (n <= 10), solve the Lyapunov equation
 *   A^T P + P A + Q = 0 using the Bartels-Stewart algorithm, then
 *   verify the remaining KYP conditions.
 *
 *   For larger systems, use an iterative gradient-descent method.
 *
 * Parameters:
 *   ss      - state-space realization
 *   epsilon - numerical tolerance (positive small number)
 *
 * Returns a KYPSolution with P matrix and feasibility flag. */
KYPSolution* kyp_solve(const StateSpace* ss, double epsilon);

/* Free a KYP solution. */
void kyp_free(KYPSolution* sol);

/* ---- KYP Verification ---- */

/* Check if a state-space system satisfies the KYP (positive real) conditions.
 *
 * Verifies the frequency-domain inequality:
 *   Re[G(jw)] >= -epsilon  for all w
 *
 * by checking at a grid of frequency points. */
bool kyp_feasibility_check(const StateSpace* ss, double epsilon);

/* Convert a state-space representation to a transfer function.
 * G(s) = C*(sI - A)^(-1)*B + D
 *
 * Uses Leverrier's algorithm (Faddeev-Leverrier) to compute
 * the characteristic polynomial and the numerator coefficients.
 *
 * The resulting transfer function has den of degree n and num of
 * degree at most n (typically n-1 for strictly proper, n for D≠0). */
void kyp_to_transfer_function(const StateSpace* ss, TransferFunction* tf);

/* ---- Passivity Verification ---- */

/* Verify if G(s) is strictly positive real (SPR).
 * A system is SPR if:
 *   1. G(s) is analytic in Re[s] >= 0
 *   2. Re[G(jw)] > 0 for all w
 *   3. lim_{w->inf} w^2 * Re[G(jw)] > 0 (for strictly proper)
 *
 * Returns: 1 if SPR, 0 if not, -1 on error. */
int kyp_verify_passivity(const StateSpace* ss);

/* ---- Lyapunov Equation Solver ---- */

/* Solve the continuous-time Lyapunov equation:
 *   A^T * X + X * A = -Q
 *
 * where Q is assumed positive semi-definite. X is returned in P.
 * Uses the Bartels-Stewart algorithm with real Schur decomposition.
 * For small systems (n <= 6), uses a direct method. */
int kyp_solve_lyapunov(const double *A, double *X, const double *Q, int n);

/* ---- Spectral Factorization ---- */

/* Perform spectral factorization: given a positive real transfer
 * function G(s), find a state-space realization (A, B, C, D) such
 * that G(s) + G(-s) = H(s) * H(-s) for some H(s). */
int kyp_spectral_factor(const TransferFunction* tf,
                         StateSpace* ss_out);

/* ---- Constants ---- */
#define KYP_EPS        1e-6    /* default numerical tolerance       */
#define KYP_MAX_ITER   200     /* max iterations for LMI solver     */
#define KYP_SMALL_N      6     /* threshold for direct methods      */

#endif /* KYP_LEMMA_H */
