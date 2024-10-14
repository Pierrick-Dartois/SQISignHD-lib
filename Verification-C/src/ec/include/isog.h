#ifndef _ISOG_H_
#define _ISOG_H_
//#include <sqisign_namespace.h>
#include <ec.h>

/* KPS structure for isogenies of degree 2 or 4 */
// Common subexpressions for a 2 or 4-isogeny.
typedef struct
{
    ec_point_t K;
} ec_kps2_t;
typedef struct
{
    ec_point_t K[3];
} ec_kps4_t;
// Common subexpressions for a (2d+1)-isogeny
// Meant to encode kernel points x([i]P) for i\in{1,...,d} for an (2d+1)-isogeny
typedef struct
{
    ec_point_t *K;
} ec_kps_t;
// Structure for a 2-isogeny chain
typedef struct
{
    ec_curve_t domain;
    ec_curve_t codomain;
    ec_kps2_t *kps;
    ec_point_t *A24;// Also contains A24 of domain and codomain
    bool *is_singular;// is_singular[i] is True if the i-th isogeny of the chain has kernel (0,0) 
    unsigned int len;
} ec_2_isog_chain_t;
// Structure for a (2d+1)-isogeny chain
typedef struct
{
    ec_curve_t domain;
    ec_curve_t codomain;
    ec_kps_t *kps;
    ec_point_t *A24; // Also contains A24 of domain and codomain
    unsigned int len;
    unsigned int d;
} ec_odd_isog_chain_t;


void xisog_2(ec_kps2_t *kps, ec_point_t *B, const ec_point_t P); // degree-2 isogeny construction
void xisog_2_singular(ec_kps2_t *kps, ec_point_t *B24, ec_point_t A24);

void xisog_4(ec_kps4_t *kps, ec_point_t *B, const ec_point_t P); // degree-4 isogeny construction
void xisog_4_singular(ec_kps4_t *kps, ec_point_t *B24, const ec_point_t P, ec_point_t A24);

void xeval_2(ec_point_t *R, ec_point_t *const Q, const int lenQ, const ec_kps2_t *kps);
void xeval_2_singular(ec_point_t *R, const ec_point_t *Q, const int lenQ, const ec_kps2_t *kps);

void xeval_4(ec_point_t *R, const ec_point_t *Q, const int lenQ, const ec_kps4_t *kps);
void xeval_4_singular(ec_point_t *R,
                      const ec_point_t *Q,
                      const int lenQ,
                      const ec_point_t P,
                      const ec_kps4_t *kps);

void xisog_3(ec_kps_t *kps, ec_point_t *B24, const ec_point_t P);
void xeval_3(ec_point_t *R, ec_point_t *const Q, const int lenQ, const ec_kps_t *kps);

void xisog_odd(ec_kps_t *kps, ec_point_t *B24, const ec_point_t P, const ec_point_t *A24, ec_point_t P2, unsigned int d);
void xeval_odd(ec_point_t *R, ec_point_t *const Q, const int lenQ, const ec_kps_t *kps, unsigned int d);

// Strategy-based 4-isogeny chain
/*
static void ec_eval_even_strategy(ec_curve_t *image,
                                  ec_point_t *points,
                                  unsigned short points_len,
                                  ec_point_t *A24,
                                  const ec_point_t *kernel,
                                  const int isog_len);
*/

void ec_2_isog_chain(ec_2_isog_chain_t *chain,const ec_point_t *kernel, const ec_curve_t *domain, unsigned int len,
    const unsigned int *strategy);

void ec_eval_2_isog_chain(ec_point_t *Q, const ec_point_t *P, const ec_2_isog_chain_t *chain);

void ec_odd_isog_chain(ec_odd_isog_chain_t *chain,const ec_point_t *kernel, const ec_curve_t *domain, unsigned int l, 
    unsigned int len, const unsigned int *strategy);

void ec_eval_odd_isog_chain(ec_point_t *Q, const ec_point_t *P, const ec_odd_isog_chain_t *chain);

void del_2_isog_chain(ec_2_isog_chain_t *chain);
void del_odd_isog_chain(ec_odd_isog_chain_t *chain);

#endif
