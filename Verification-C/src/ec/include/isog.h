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
    unsigned int d;
} ec_kps_t;

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

void xisog_odd(ec_kps_t *kps, ec_point_t *B24, const ec_point_t P, ec_point_t A24, unsigned int d);
void xeval_odd(ec_point_t *R, ec_point_t *const Q, const int lenQ, const ec_kps_t *kps);

// Strategy-based 4-isogeny chain
static void ec_eval_even_strategy(ec_curve_t *image,
                                  ec_point_t *points,
                                  unsigned short points_len,
                                  ec_point_t *A24,
                                  const ec_point_t *kernel,
                                  const int isog_len);

#endif
