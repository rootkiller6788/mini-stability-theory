#ifndef ISS_INTERCONNECTED_H
#define ISS_INTERCONNECTED_H
#include "iss_core.h"
#include "iss_gain.h"

/* ==============================================================
 * iss_interconnected.h ? ISS for Interconnected Systems
 *
 * Cascade: dx1=f1(x1,x2), dx2=f2(x2)
 * If x1-subsystem is ISS wrt x2, and x2-subsystem is GAS,
 * then the cascade is GAS.
 *
 * Feedback: dx1=f1(x1,x2), dx2=f2(x2,x1)
 * Requires small-gain condition for ISS.
 *
 * References:
 *   Sontag & Teel (1995) Changing Supply Functions
 *   Jiang, Teel, Praly (1994) Small-Gain Theorem for ISS
 * ============================================================== */

/* Cascade system */
typedef struct {
    ISS_System* driving;      /* S2: GAS subsystem */
    ISS_System* driven;       /* S1: ISS wrt x2 */
    bool is_ISS_cascade;
    ISS_Gain cascade_gain;
} ISS_Cascade;

/* Feedback interconnection */
typedef struct {
    ISS_System* system1;
    ISS_System* system2;
    double coupling_strength;
    bool is_ISS_feedback;
    bool small_gain_satisfied;
} ISS_Feedback;

/* Build cascade connection */
ISS_Cascade* iss_cascade_create(ISS_System* driven, ISS_System* driving);
void iss_cascade_free(ISS_Cascade* cascade);
bool iss_cascade_is_ISS(const ISS_Cascade* cascade);

/* Build feedback connection */
ISS_Feedback* iss_feedback_create(ISS_System* s1, ISS_System* s2, double coupling);
void iss_feedback_free(ISS_Feedback* fb);
bool iss_feedback_is_ISS(const ISS_Feedback* fb);

/* Nonlinear small-gain for feedback */
bool iss_feedback_small_gain(const ISS_Feedback* fb, double tol);

/* Lyapunov-based ISS for cascade */
ISS_Certificate iss_cascade_certificate(const ISS_Cascade* cascade);
double iss_cascade_gain(const ISS_Cascade* c);
int iss_cascade_count(void);
bool iss_feedback_has_finite_gain(const ISS_Feedback* fb);

/* Composit Lyapunov for feedback */
ISS_Certificate iss_feedback_certificate(const ISS_Feedback* fb);

/* Cascade verification with ISS gains */
bool iss_cascade_verify_ISS(ISS_Cascade* c, double tol);
bool iss_parallel_is_ISS(ISS_System* s1, ISS_System* s2);
ISS_System** iss_cascade_chain_create(ISS_System** systems, int n_systems);
void iss_cascade_chain_free(ISS_System** chain);
bool iss_feedback_lyapunov_check(ISS_Feedback* fb, double tolerance);
bool iss_interconnected_trajectory_check(ISS_Cascade* c,
    const ISS_State* x1_0, const ISS_State* x2_0,
    double dt, int steps, double u_max);
ISS_Certificate iss_interconnection_certificate(ISS_System** systems,
    int n_sys, const int* connections, int n_conn);

#endif




























































