#ifndef INTERCON_PASSIVITY_H
#define INTERCON_PASSIVITY_H
#include "intercon_core.h"

/* ==============================================================
 * Passivity-Based Stability (Willems 1972, Hill & Moylan 1976)
 *
 * A system is passive if there exists a storage function V(x) >= 0
 * such that V(x(t)) - V(x(0)) <= int_0^t u^T(s)*y(s) ds
 *
 * Passivity theorem: feedback interconnection of two passive
 * systems is stable (under detectability conditions).
 * ============================================================== */

typedef enum {
    PASSIVE_NONE         = 0,
    PASSIVE_LOSSY        = 1,  /* input strictly passive */
    PASSIVE_OUTPUT_STRICT = 2, /* output strictly passive */
    PASSIVE_STRICT       = 3,  /* strictly passive (both) */
    PASSIVE_NONEXPANSIVE = 4,  /* ||y|| <= ||u|| */
    PASSIVE_ACTIVE       = 5   /* not passive */
} PassivityType;

typedef struct {
    PassivityType type;
    double        excess;     /* passivity excess/shortage */
    double        shortage;   /* negative = excess */
    bool          is_passive;
    double        dissipation; /* dissipation rate */
} PassivityResult;

/* ==============================================================
 * Passivity Tests
 * ============================================================== */

/* KYP Lemma: LTI system is passive iff exists P > 0 such that
 * [A^T*P + P*A    P*B - C^T]  <= 0
 * [B^T*P - C      -(D + D^T)] */
PassivityResult* ic_passivity_test(const ICSubsystem* s);
void             ic_passivity_result_free(PassivityResult* pr);
void             ic_passivity_result_print(const PassivityResult* pr);

/* ==============================================================
 * Passivity-Based Interconnection Stability
 * ============================================================== */

/* Check if feedback interconnection of two passive systems is stable */
bool ic_passivity_feedback_stable(const ICSubsystem* s1,
                                    const ICSubsystem* s2);

/* Passivity index: find the maximum nu such that system + nu*I is passive */
double ic_passivity_index(const ICSubsystem* s);

/* Output feedback passivity: system is output-feedback passive if
 * there exists output feedback u = -K*y + v making system passive */
double ic_output_feedback_passivity(const ICSubsystem* s);

/* ==============================================================
 * Storage Function Construction
 * ============================================================== */

/* Compute storage function P for passive LTI system via KYP */
void ic_storage_function(const ICSubsystem* s, double* P, int max_dim);

/* Compute available storage (max energy extractable) */
double ic_available_storage(const ICSubsystem* s);

/* Compute required supply (min energy to reach state) */
double ic_required_supply(const ICSubsystem* s, const double* x);

/* ==============================================================
 * Dissipativity (Generalized Passivity)
 *
 * System is (Q,S,R)-dissipative if:
 * int_0^t w(u(s), y(s)) ds >= V(x(t)) - V(x(0))
 * with supply rate w(u,y) = y^T*Q*y + 2*y^T*S*u + u^T*R*u
 * ============================================================== */

typedef struct {
    double Q[IC_MAX_DIM][IC_MAX_DIM];  /* output weighting */
    double S[IC_MAX_DIM][IC_MAX_DIM];  /* cross weighting */
    double R[IC_MAX_DIM][IC_MAX_DIM];  /* input weighting */
    int    dim;
} SupplyRate;

typedef struct {
    bool is_dissipative;
    double P[IC_MAX_DIM][IC_MAX_DIM]; /* storage matrix */
    int   n;
} DissipativityResult;

SupplyRate* ic_supply_rate_create(int dim, const double* Q, const double* S,
                                    const double* R);
void        ic_supply_rate_free(SupplyRate* sr);
DissipativityResult* ic_dissipativity_test(const ICSubsystem* s,
                                             const SupplyRate* sr);
void        ic_dissipativity_result_free(DissipativityResult* dr);
void        ic_dissipativity_result_print(const DissipativityResult* dr);

/* Check if interconnection of (Q1,S1,R1) and (Q2,S2,R2) is stable */
bool ic_dissipativity_intercon_stable(const SupplyRate* sr1,
                                        const SupplyRate* sr2);

#endif /* INTERCON_PASSIVITY_H */
double ic_passivity_l2_gain_bound(const ICSubsystem* s);
bool   ic_is_finite_gain_stable(const ICSubsystem* s);
bool   ic_mixed_passivity_small_gain(const ICSubsystem* s1, const ICSubsystem* s2);
bool   ic_incremental_passivity_test(const ICSubsystem* s);

bool ic_output_strict_passivity_margin(const ICSubsystem* s, double* rho);
double ic_passivity_sector_max(const ICSubsystem* s);
bool ic_feedback_preserves_passivity(const ICSubsystem* s1, const ICSubsystem* s2);
