# Input-Output Stability Theory

## Zames and Sandberg Framework

Input-output stability analyzes systems as operators mapping between normed signal spaces (L2, L-infinity, etc.). A system H is L2-stable if u in L2 implies Hu in L2.

## Small-Gain Theorem

The feedback interconnection of H1 and H2 is L2-stable if gamma(H1)*gamma(H2) < 1, where gamma is the L2 operator gain (the H-infinity norm for LTI systems).

## Passivity

A system is passive if (u, Hu) >= 0 for all u (net energy absorbed). Passivity is preserved under negative feedback interconnection. Strict passivity ensures asymptotic stability.

## Conic Sectors

For Lur'e systems, the circle and Popov criteria provide frequency-domain tests for absolute stability. The KYP lemma connects frequency-domain inequalities to state-space LMIs.
