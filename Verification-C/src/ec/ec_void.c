#include "include/curve_extras.h"
//#include "tedwards.h"
//#include "ec_params.h"
//#include <assert.h>
//The latter being replaced by
#include "../common/include/constants.h"
#include "../gf/include/fp2.h"
#include "include/ec.h"

bool ec_is_zero(ec_point_t const* P)
{
    return fp2_is_zero(&P->z);
}

void ec_set_zero(ec_point_t *P)
{
    fp2_set(&(P->x), 1);
    fp2_set(&(P->z), 0);
    return;
}

void ec_init(ec_point_t* P)
{ // Initialize point as identity element (1:0)
    ec_set_zero(P);
}


void copy_curve(ec_curve_t* E1, ec_curve_t const* E2)
{
    fp2_copy(&(E1->A), &(E2->A));
    fp2_copy(&(E1->C), &(E2->C));
}

void xDBL(ec_point_t* Q, ec_point_t const* P, ec_point_t const* AC)
{
    // This version computes the coefficient values A+2C and 4C on-the-fly 
    // The curve coefficients are passed via AC = (A:C)
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_add(&t3, &AC->z, &AC->z);  
    fp2_mul(&t1, &t1, &t3);
    fp2_add(&t1, &t1, &t1);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_add(&t0, &t3, &AC->x);
    fp2_mul(&t0, &t0, &t2);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void xDBLv2(ec_point_t* Q, ec_point_t const* P, ec_point_t const* A24)
{
    // This version receives the coefficient value A24 = (A+2C:4C) 
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    fp2_mul(&t1, &t1, &A24->z);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_mul(&t0, &t2, &A24->x);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

void xADD(ec_point_t* R, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ)
{
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_add(&t2, &Q->x, &Q->z);
    fp2_sub(&t3, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t3);
    fp2_mul(&t1, &t1, &t2);
    fp2_add(&t2, &t0, &t1);
    fp2_sub(&t3, &t0, &t1);
    fp2_sqr(&t2, &t2);
    fp2_sqr(&t3, &t3);
    fp2_mul(&t2, &PQ->z, &t2);
    fp2_mul(&R->z, &PQ->x, &t3);
    fp2_copy(&R->x, &t2);
}

void xDBLADD(ec_point_t* R, ec_point_t* S, ec_point_t const* P, ec_point_t const* Q, ec_point_t const* PQ, ec_point_t const* A24)
{
    // Requires precomputation of A24 = (A+2C:4C)
    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&R->x, &t0);
    fp2_sub(&t2, &Q->x, &Q->z);
    fp2_add(&S->x, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t2);
    fp2_sqr(&R->z, &t1);
    fp2_mul(&t1, &t1, &S->x);
    fp2_sub(&t2, &R->x, &R->z);
    fp2_mul(&R->z, &R->z, &A24->z);
    fp2_mul(&R->x, &R->x, &R->z);
    fp2_mul(&S->x, &A24->x, &t2);
    fp2_sub(&S->z, &t0, &t1);
    fp2_add(&R->z, &R->z, &S->x);
    fp2_add(&S->x, &t0, &t1);
    fp2_mul(&R->z, &R->z, &t2);
    fp2_sqr(&S->z, &S->z);
    fp2_sqr(&S->x, &S->x);
    fp2_mul(&S->z, &S->z, &PQ->x);
    fp2_mul(&S->x, &S->x, &PQ->z);
}

bool is_point_equal(const ec_point_t* P, const ec_point_t* Q)
{ // Evaluate if two points in Montgomery coordinates (X:Z) are equal
  // Returns 1 (true) if P=Q, 0 (false) otherwise
    fp2_t t0, t1;

    if ((fp2_is_zero(&P->x) && fp2_is_zero(&P->z)) || (fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z))) {
        return fp2_is_zero(&P->x) && fp2_is_zero(&P->z) && fp2_is_zero(&Q->x) && fp2_is_zero(&Q->z);
    }

    fp2_mul(&t0, &P->x, &Q->z);
    fp2_mul(&t1, &Q->x, &P->z);
    fp2_sub(&t0, &t0, &t1);

    return fp2_is_zero(&t0);
}

void swap_points(ec_point_t* P, ec_point_t* Q, const digit_t option)
{ // Swap points
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;

    for (int i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->x.re[i] ^ Q->x.re[i]);
        P->x.re[i] = temp ^ P->x.re[i];
        Q->x.re[i] = temp ^ Q->x.re[i];
        temp = option & (P->x.im[i] ^ Q->x.im[i]);
        P->x.im[i] = temp ^ P->x.im[i];
        Q->x.im[i] = temp ^ Q->x.im[i];
        temp = option & (P->z.re[i] ^ Q->z.re[i]);
        P->z.re[i] = temp ^ P->z.re[i];
        Q->z.re[i] = temp ^ Q->z.re[i];
        temp = option & (P->z.im[i] ^ Q->z.im[i]);
        P->z.im[i] = temp ^ P->z.im[i];
        Q->z.im[i] = temp ^ Q->z.im[i];
    }
}

void copy_point(ec_point_t* P, ec_point_t const* Q)
{
    fp2_copy(&(P->x), &(Q->x));
    fp2_copy(&(P->z), &(Q->z));
}

void ec_normalize(ec_point_t* P){
    fp2_inv(&P->z);
    fp2_mul(&P->x, &P->x, &P->z);
    fp_mont_setone(P->z.re);
    fp_set(P->z.im, 0);
}

void ec_neg(ec_point_t* res, const ec_point_t* P){
    // DOES NOTHING
    copy_point(res, P);
}

void point_to_curve(ec_curve_t* curve, ec_point_t const * AC){
    fp2_copy(&curve->A,&AC->x);
    fp2_copy(&curve->C,&AC->z);
}

void AC_to_A24(ec_point_t* A24, ec_point_t const * AC){
    // Computation of A24=(A+2C:4C)
    fp2_add(&A24->x, &AC->z, &AC->z);    
    fp2_add(&A24->z, &A24->x, &A24->x);
    fp2_add(&A24->x, &A24->x, &AC->x);
}

void curve_to_A24(ec_point_t* A24, ec_curve_t const * curve){
    // Computation of A24=(A+2C:4C)
    fp2_add(&A24->x, &curve->C, &curve->C);
    fp2_add(&A24->z, &A24->x, &A24->x);
    fp2_add(&A24->x, &A24->x, &curve->A);
}


void xMUL(ec_point_t* Q, ec_point_t const* P, digit_t const* k, ec_curve_t const* curve)
{
    ec_point_t R0, R1, A24;
    digit_t mask;
    unsigned int bit = 0, prevbit = 0, swap;

    fp2_add(&A24.x, &curve->C, &curve->C);    // Precomputation of A24=(A+2C:4C)
    fp2_add(&A24.z, &A24.x, &A24.x);
    fp2_add(&A24.x, &A24.x, &curve->A);

    // R0 <- (1:0), R1 <- P
    ec_init(&R0);
    fp2_copy(&R1.x, &P->x);
    fp2_copy(&R1.z, &P->z);

    // Main loop
    // number of bits BITS = NWORDS_FIELD*RADIX
    for (int i = NWORDS_FIELD*RADIX-1; i >= 0; i--) {                                          
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;                         
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, &R0, &R1, P, &A24);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(&R0, &R1, mask);

    fp2_copy(&Q->x, &R0.x);
    fp2_copy(&Q->z, &R0.z);
}

void xMULv2(ec_point_t* Q, ec_point_t const* P, digit_t const* k, const int kbits, ec_point_t const* A24)
{
    // This version receives the coefficient value A24 = (A+2C:4C) 
    ec_point_t R0, R1;
    digit_t mask;
    unsigned int bit = 0, prevbit = 0, swap;

    // R0 <- (1:0), R1 <- P
    ec_init(&R0);
    fp2_copy(&R1.x, &P->x);
    fp2_copy(&R1.z, &P->z);

    // Main loop
    for (int i = kbits-1; i >= 0; i--) {                              
        bit = (k[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;                      
        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(&R0, &R1, mask);
        xDBLADD(&R0, &R1, &R0, &R1, P, A24);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(&R0, &R1, mask);

    fp2_copy(&Q->x, &R0.x);
    fp2_copy(&Q->z, &R0.z);
}

static void mp_add(digit_t* c, const digit_t* a, const digit_t* b, const unsigned int nwords)
{ // Multiprecision addition
    unsigned int i, carry = 0;

    for (i = 0; i < nwords; i++) {
        ADDC(c[i], carry, a[i], b[i], carry);
    }
}

static void mp_sub(digit_t* c, digit_t const* a, digit_t const* b, const unsigned int nwords)
{ // Multiprecision subtraction, assuming a > b
    unsigned int i, borrow = 0;

    for (i = 0; i < nwords; i++) {
        SUBC(c[i], borrow, a[i], b[i], borrow);
    }
}

void select_ct(digit_t* c, const digit_t* a, const digit_t* b, const digit_t mask, const int nwords)
{ // Select c <- a if mask = 0, select c <- b if mask = 1...1

    for (int i = 0; i < nwords; i++) {
        c[i] = ((a[i] ^ b[i]) & mask) ^ a[i];
    }
}

void swap_ct(digit_t* a, digit_t* b, const digit_t option, const int nwords)
{ // Swap entries
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then a <- b and b <- a
    digit_t temp;

    for (int i = 0; i < nwords; i++) {
        temp = option & (a[i] ^ b[i]);
        a[i] = temp ^ a[i];
        b[i] = temp ^ b[i];
    }
}

// Compute S = k*P + l*Q, with PQ = P+Q
void xDBLMUL(ec_point_t* S, ec_point_t const* P, digit_t const* k, ec_point_t const* Q, digit_t const* l, ec_point_t const* PQ, ec_curve_t const* curve)
{
    int i;
    digit_t evens, mevens, bitk0, bitl0, maskk, maskl, temp, bs1_ip1, bs2_ip1, bs1_i, bs2_i, h;
    digit_t sigma[2] = {0}, pre_sigma = 0;
    digit_t k_t[NWORDS_FIELD], l_t[NWORDS_FIELD], one[NWORDS_FIELD] = {0}, r[2*NWORDS_FIELD*RADIX] = {0};            
    ec_point_t A24, DIFF1a, DIFF1b, DIFF2a, DIFF2b, R[3] = {0}, T[3];

    // Derive sigma according to parity
    bitk0 = (k[0] & 1);
    bitl0 = (l[0] & 1);
    maskk = 0 - bitk0;               // Parity masks: 0 if even, otherwise 1...1
    maskl = 0 - bitl0;
    sigma[0] = (bitk0 ^ 1);
    sigma[1] = (bitl0 ^ 1);
    evens = sigma[0] + sigma[1];     // Count number of even scalars
    mevens = 0 - (evens & 1);        // Mask mevens <- 0 if # even scalars = 0 or 2, otherwise mevens = 1...1

    // If k and l are both even or both odd, pick sigma = (0,1)
    sigma[0] = (sigma[0] & mevens);
    sigma[1] = (sigma[1] & mevens) | (1 & ~mevens);

    // Convert even scalars to odd
    one[0] = 1;
    mp_sub(k_t, k, one, NWORDS_FIELD);
    mp_sub(l_t, l, one, NWORDS_FIELD);
    select_ct(k_t, k_t, k, maskk, NWORDS_FIELD);                                        
    select_ct(l_t, l_t, l, maskl, NWORDS_FIELD);

    // Scalar recoding
    for (i = 0; i < NWORDS_FIELD*RADIX; i++) {
        // If sigma[0] = 1 swap k_t and l_t
        maskk = 0 - (sigma[0] ^ pre_sigma);
        swap_ct(k_t, l_t, maskk, NWORDS_FIELD);
        
        if (i == NWORDS_FIELD*RADIX-1) {
            bs1_ip1 = 0;
            bs2_ip1 = 0;
        } else {
            bs1_ip1 = mp_shiftr(k_t, 1, NWORDS_FIELD);
            bs2_ip1 = mp_shiftr(l_t, 1, NWORDS_FIELD);
        }
        bs1_i = k_t[0] & 1;
        bs2_i = l_t[0] & 1;

        r[2*i]   = bs1_i ^ bs1_ip1;
        r[2*i+1] = bs2_i ^ bs2_ip1;

        // Revert sigma if second bit, r_(2i+1), is 1 
        pre_sigma = sigma[0];
        maskk = 0 - r[2*i+1];
        select_ct(&temp, &sigma[0], &sigma[1], maskk, 1);
        select_ct(&sigma[1], &sigma[1], &sigma[0], maskk, 1);
        sigma[0] = temp;
    }

    // Point initialization
    ec_init(&R[0]);
    maskk = 0 - sigma[0];
    select_ct((digit_t*)&R[1], (digit_t*)P, (digit_t*)Q, maskk, 4*NWORDS_FIELD);
    select_ct((digit_t*)&R[2], (digit_t*)Q, (digit_t*)P, maskk, 4*NWORDS_FIELD);
    fp2_copy(&DIFF1a.x, &R[1].x);
    fp2_copy(&DIFF1a.z, &R[1].z);
    fp2_copy(&DIFF1b.x, &R[2].x);
    fp2_copy(&DIFF1b.z, &R[2].z);

    // Initialize DIFF2a <- P+Q, DIFF2b <- P-Q
    xADD(&R[2], &R[1], &R[2], PQ);
    fp2_copy(&DIFF2a.x, &R[2].x);
    fp2_copy(&DIFF2a.z, &R[2].z);
    fp2_copy(&DIFF2b.x, &PQ->x);
    fp2_copy(&DIFF2b.z, &PQ->z);

    fp2_add(&A24.x, &curve->C, &curve->C);    // Precomputation of A24=(A+2C:4C)
    fp2_add(&A24.z, &A24.x, &A24.x);
    fp2_add(&A24.x, &A24.x, &curve->A);

    // Main loop
    for (i = NWORDS_FIELD*RADIX-1; i>=0; i--) {
        h = r[2*i] + r[2*i+1];    // in {0, 1, 2}
        maskk = 0 - (h & 1);
        select_ct((digit_t*)&T[0], (digit_t*)&R[0], (digit_t*)&R[1], maskk, 4*NWORDS_FIELD);
        maskk = 0 - (h >> 1);
        select_ct((digit_t*)&T[0], (digit_t*)&T[0], (digit_t*)&R[2], maskk, 4*NWORDS_FIELD);
        xDBLv2(&T[0], &T[0], &A24);

        maskk = 0 - r[2*i+1];     // in {0, 1}
        select_ct((digit_t*)&T[1], (digit_t*)&R[0], (digit_t*)&R[1], maskk, 4*NWORDS_FIELD);
        select_ct((digit_t*)&T[2], (digit_t*)&R[1], (digit_t*)&R[2], maskk, 4*NWORDS_FIELD);
        swap_points(&DIFF1a, &DIFF1b, maskk);
        xADD(&T[1], &T[1], &T[2], &DIFF1a);
        xADD(&T[2], &R[0], &R[2], &DIFF2a);
        
        // If hw (mod 2) = 1 then swap DIFF2a and DIFF2b
        maskk = 0 - (h & 1);
        swap_points(&DIFF2a, &DIFF2b, maskk);

        // R <- T
        memcpy((digit_t*)&R[0], (digit_t*)&T[0], NWORDS_FIELD*RADIX*4/8);
        memcpy((digit_t*)&R[1], (digit_t*)&T[1], NWORDS_FIELD*RADIX*4/8);
        memcpy((digit_t*)&R[2], (digit_t*)&T[2], NWORDS_FIELD*RADIX*4/8);
    }

    // Output R[evens]
    select_ct((digit_t*)S, (digit_t*)&R[0], (digit_t*)&R[1], mevens, 4*NWORDS_FIELD);
    maskk = 0 - (bitk0 & bitl0);
    select_ct((digit_t*)S, (digit_t*)S, (digit_t*)&R[2], maskk, 4*NWORDS_FIELD);
}

// Computes P+m*Q
void ec_ladder3pt(ec_point_t *R, digit_t const* m, ec_point_t const *P, ec_point_t const *Q, ec_point_t const *PQ, ec_curve_t const *A)
{
    // Curve constant in the form A24=(A+2C:4C)
    ec_point_t A24;
    fp2_add(&A24.z, &A->C, &A->C);
    fp2_add(&A24.x, &A->A, &A24.z);
    fp2_add(&A24.z, &A24.z, &A24.z);

    ec_point_t X0, X1, X2;
    copy_point(&X0, Q);
    copy_point(&X1, P);
    copy_point(&X2, PQ);

    int i,j;
    uint64_t t;
    for (i = 0; i < NWORDS_FIELD; i++)
    {
        t = 1;
        for (j = 0 ; j < 64; j++)
        {
            swap_points(&X1, &X2, -((t & m[i]) == 0));
            xDBLADD(&X0, &X1, &X0, &X1, &X2, &A24);
            swap_points(&X1, &X2, -((t & m[i]) == 0));
            t <<= 1;
        };
    };
    copy_point(R, &X1);
}

void ec_j_inv(fp2_t* j_inv, const ec_curve_t* curve){
    /* j-invariant computation for montgommery coefficient A24=(A+2C:4C) */
    fp2_t t0, t1;
    
    fp2_tomont(j_inv,j_inv);// We have to put the j-invariant in Montgomery form (why?)
    fp2_sqr(&t1, &curve->C);//t1=C^2
    fp2_sqr(j_inv, &curve->A);//j_inv=A^2
    fp2_add(&t0, &t1, &t1);//t0=2C^2
    fp2_sub(&t0, j_inv, &t0);//t0=A^2-2C^2
    fp2_sub(&t0, &t0, &t1);//t0=A^2-3C^2
    fp2_sub(j_inv, &t0, &t1);//j_inv=A^2-4C^2
    fp2_sqr(&t1, &t1);//t1=C^4
    fp2_mul(j_inv, j_inv, &t1);//j_inv=(A^2-4C^2)*C^4
    fp2_add(&t0, &t0, &t0);//t0=2(A^2-3C^2)
    fp2_add(&t0, &t0, &t0);//t0=4(A^2-3C^2)
    fp2_sqr(&t1, &t0);//t1=16(A^2-3C^2)^2
    fp2_mul(&t0, &t0, &t1);//t0=64(A^2-3C^2)^3
    fp2_add(&t0, &t0, &t0);//t0=128(A^2-3C^2)^3
    fp2_add(&t0, &t0, &t0);//t0=256(A^2-3C^2)^3
    fp2_inv(j_inv);//j_inv=1/((A^2-4C^2)*C^4)
    fp2_mul(j_inv, &t0, j_inv);//j_inv=256(A^2-3C^2)^3/((A^2-4C^2)*C^4)
    fp2_frommont(j_inv,j_inv);
}
