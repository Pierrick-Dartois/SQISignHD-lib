#ifndef _BIEXT_H_
#define _BIEXT_H_

//#include <sqisign_namespace.h>
#include <ec.h>

typedef struct weil_params
{
    uint64_t e;     // When points have order 2^e
    uint64_t *n;    // When points have order n (for the generic case only)
    unsigned int nwords; // Number of words of the multiprecision integer n
    ec_point_t P;   // x(P)
    ec_point_t Q;   // x(Q)
    ec_point_t PQ;  // x(P+Q) = (PQX/PQZ : 1)
    fp2_t ixP;      // PZ/PX
    fp2_t ixQ;      // QZ/QX
    fp2_t ixPQ;     // z(P+Q)/x(P+Q) (for the generic case only)
    ec_point_t A24; // ((A+2)/4 : 1)
} weil_params_t;

// For two bases <P, Q> and <R, S> store:
// x(P - R), x(P - S), x(R - Q), x(S - Q)
typedef struct weil_dlog_diff_points
{
    ec_point_t PmR; // x(P - R)
    ec_point_t PmS; // x(P - S)
    ec_point_t RmQ; // x(R - Q)
    ec_point_t SmQ; // x(S - Q)
} weil_dlog_diff_points_t;

typedef struct weil_dlog_params
{
    uint64_t e;                   // Points have order 2^e
    ec_basis_t PQ;                // x(P), x(Q), x(P-Q)
    ec_basis_t RS;                // x(R), x(S), x(R-S)
    weil_dlog_diff_points_t diff; // x(P - R), x(P - S), x(R - Q), x(S - Q)
    fp2_t ixP;                    // PZ/PX
    fp2_t ixQ;                    // QZ/QX
    fp2_t ixR;                    // RZ/RX
    fp2_t ixS;                    // SZ/SX
    ec_point_t A24;               // ((A+2)/4 : 1)
} weil_dlog_params_t;

// Computes e = e_{2^e}(P, Q) using biextension ladder
void weil_2e(fp2_t *r, uint64_t e, ec_point_t *P, ec_point_t *Q, ec_point_t *PQ, ec_curve_t *E);

// Computes e = e_n(P, Q) using biextension ladder
void weil(fp2_t *r, uint64_t *n, unsigned int nwords, ec_point_t *P, ec_point_t *Q, ec_point_t *PQ, ec_curve_t *E);

// Given two bases <P, Q> and <R, S> computes scalars
// such that R = [r1]P + [r2]Q, S = [s1]P + [s2]Q
void ec_dlog_2_weil(digit_t *r1,
                    digit_t *r2,
                    digit_t *s1,
                    digit_t *s2,
                    ec_basis_t *PQ,
                    ec_basis_t *RS,
                    ec_curve_t *curve,
                    int e,
                    unsigned int nwords);

// Given a basis (P, Q) of E[l^e] and R\in E[l^e] computes scalars
// such that R = [r1]P + [r2]Q.
void ec_single_dlog_le_weil(digit_t *r1,
               digit_t *r2,
               ec_basis_t *PQ,
               ec_point_t *R,
               ec_curve_t *curve,
               int l,
               int e, 
               unsigned int nwords);

#endif
