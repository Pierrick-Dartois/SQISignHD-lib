#include "ec.h"
#include "fp2.h"

#include "gf_constants.h"
#define LEN_NQR_TABLE (sizeof(NQR_TABLE) / sizeof(*NQR_TABLE))

#include <assert.h>
static_assert(LEN_NQR_TABLE == (sizeof(Z_NQR_TABLE) / sizeof(*Z_NQR_TABLE)),
              "tables of nonsquares should have the same length");

// New methods for basis generation using Entangled / ApresSQI like
// methods. Finds a point of full order with square checking, then
// finds the point of desired order by clearing cofactors. Not need
// to double all the way down to ensure the correct order.
//
// This also allows a faster method for completing a torsion basis
// as if we know a point P is above (0 : 0) or not, we can directly
// compute the point Q orthogonal to this one using
// ec_curve_to_point_2f_above_montgomery or
// ec_curve_to_point_2f_not_above_montgomery

static void
difference_point(ec_point_t *PQ, const ec_point_t *P, const ec_point_t *Q, const ec_curve_t *curve)
{
    // Given P,Q in projective x-only, computes a deterministic choice for (P-Q)
    // Based on Proposition 3 of https://eprint.iacr.org/2017/518.pdf

    fp2_t Bxx, Bxz, Bzz, t0, t1;

    fp2_mul(&t0, &P->x, &Q->x);
    fp2_mul(&t1, &P->z, &Q->z);
    fp2_sub(&Bxx, &t0, &t1);
    fp2_sqr(&Bxx, &Bxx);
    fp2_mul(&Bxx, &Bxx, &curve->C); // C*(P.x*Q.x-P.z*Q.z)^2
    fp2_add(&Bxz, &t0, &t1);
    fp2_mul(&t0, &P->x, &Q->z);
    fp2_mul(&t1, &P->z, &Q->x);
    fp2_add(&Bzz, &t0, &t1);
    fp2_mul(&Bxz, &Bxz, &Bzz); // (P.x*Q.x+P.z*Q.z)(P.x*Q.z+P.z*Q.x)
    fp2_sub(&Bzz, &t0, &t1);
    fp2_sqr(&Bzz, &Bzz);
    fp2_mul(&Bzz, &Bzz, &curve->C); // C*(P.x*Q.z-P.z*Q.x)^2
    fp2_mul(&Bxz, &Bxz, &curve->C); // C*(P.x*Q.x+P.z*Q.z)(P.x*Q.z+P.z*Q.x)
    fp2_mul(&t0, &t0, &t1);
    fp2_mul(&t0, &t0, &curve->A);
    fp2_add(&t0, &t0, &t0);
    fp2_add(&Bxz, &Bxz, &t0); // C*(P.x*Q.x+P.z*Q.z)(P.x*Q.z+P.z*Q.x) + 2*A*P.x*Q.z*P.z*Q.x

    // To ensure that the denominator is a fourth power in Fp, we normalize by
    // C*C_bar^2*(P.z)_bar^2*(Q.z)_bar^2
    fp_copy(&t0.re, &curve->C.re);
    fp_neg(&t0.im, &curve->C.im);
    fp2_sqr(&t0, &t0);
    fp2_mul(&t0, &t0, &curve->C);
    fp_copy(&t1.re, &P->z.re);
    fp_neg(&t1.im, &P->z.im);
    fp2_sqr(&t1, &t1);
    fp2_mul(&t0, &t0, &t1);
    fp_copy(&t1.re, &Q->z.re);
    fp_neg(&t1.im, &Q->z.im);
    fp2_sqr(&t1, &t1);
    fp2_mul(&t0, &t0, &t1);
    fp2_mul(&Bxx, &Bxx, &t0);
    fp2_mul(&Bxz, &Bxz, &t0);
    fp2_mul(&Bzz, &Bzz, &t0);

    // Solving quadratic equation
    fp2_sqr(&t0, &Bxz);
    fp2_mul(&t1, &Bxx, &Bzz);
    fp2_sub(&t0, &t0, &t1);
    fp2_sqrt(&t0);
    fp2_add(&PQ->x, &Bxz, &t0);
    fp2_copy(&PQ->z, &Bzz);
}

static void compute_xonly_basis(ec_basis_t *PQ, const ec_curve_t *curve, const ec_point_t *P, const ec_point_t *Q){
    copy_point(&PQ->P, P);
    copy_point(&PQ->Q, Q);
    difference_point(&PQ->PmQ, P, Q, curve);
}

// Given an x-coordinate, determines if this is a valid
// point on the curve. Assumes C=1.
static uint32_t
is_on_curve(fp2_t *x, const ec_curve_t *curve)
{
    assert(fp2_is_one(&curve->C));
    fp2_t t0;

    fp_t one;
    fp_set_one(&one);

    fp2_add(&t0, x, &curve->A);   // x + (A/C)
    fp2_mul(&t0, &t0, x);         // x^2 + (A/C)*x
    fp_add(&t0.re, &t0.re, &one); // x^2 + (A/C)*x + 1
    fp2_mul(&t0, &t0, x);         // x^3 + (A/C)*x^2 + x

    return fp2_is_square(&t0);
}

// Given an x-coordinate, computes the y coordinate. Assumes C=1.
static void
y_coordinate(fp2_t *y, const fp2_t *x, const ec_curve_t *curve)
{
    assert(fp2_is_one(&curve->C));

    fp_t one;
    fp_set_one(&one);

    fp2_add(y, x, &curve->A);   // x + (A/C)
    fp2_mul(y, y, x);         // x^2 + (A/C)*x
    fp_add(&(y->re), &(y->re), &one); // x^2 + (A/C)*x + 1
    fp2_mul(y, y, x);         // x^3 + (A/C)*x^2 + x

    fp2_sqrt(y);
}

static void
normalize_point(ec_point_t *P)
{
    fp2_t z_inv;

    fp2_copy(&z_inv,&(P->z));
    fp2_inv(&z_inv);
    fp2_mul(&(P->x),&(P->x),&z_inv);
    fp2_set_one(&(P->z));
}

static void
compute_alpha(fp2_t *alpha, ec_curve_t *curve)
{
    // Compute a root of x^2 + Ax + 1
    fp2_t d, four;

    // d = sqrt(A^2 - 4)
    fp2_set_small(&four, 4);
    fp2_sqr(&d, &curve->A);
    fp2_sub(&d, &d, &four);
    fp2_sqrt(&d);

    // alpha = (-A + d) / 2
    fp2_sub(alpha, &d, &curve->A);
    fp2_half(alpha, alpha);
}

/// Finds a point of order k * 2^n where n is the largest power of two
/// dividing (p+1).
/// The x-coordinate is picked such that the point (0 : 0) is always
/// the point of order two below the point.
static uint8_t
ec_curve_to_point_2f_above_montgomery(ec_point_t *P, ec_curve_t *curve)
{
    // Compute alpha of x^2 + (A/C)x + 1
    fp2_t alpha;
    compute_alpha(&alpha, curve);

    fp_t one;
    fp_set_one(&one);

    uint_fast8_t hint = 0;
    for (;; ++hint) {
        fp2_t z1, z2;

        // collect z2-value from table, we have 20 chances
        // and expect to be correct 50% of the time.
        if (hint < LEN_NQR_TABLE) {
            z2=Z_NQR_TABLE[hint];
        }

        // Fallback method for when we're unlucky
        else {
            if (hint == LEN_NQR_TABLE) {
                fp_set_one(&z1.im);
                fp_set_one(&z2.im);
                fp_set_small(&z1.re, hint - 2);
                fp_set_small(&z2.re, hint - 1);
            }

            // Look for z2 = i + hint with z2 a square and
            // z2 - 1 not a square.
            for (;; ++hint) {
                // Set z2 = i + hint and z1 = z2 - 1
                // TODO: we could swap z1 and z2 on failure
                // and save one addition
                fp_add(&z1.re, &z1.re, &one);
                fp_add(&z2.re, &z2.re, &one);

                // Now check whether z2 is a square and z1 is not
                if (fp2_is_square(&z2) & !fp2_is_square(&z1)){
                    break;
                }

                assert(hint < UINT_FAST8_MAX);
            }
        }

        // Compute x-coordinate
        fp2_t t0, x;

        // Find a point on curve with x a NQR
        fp2_mul(&x, &z2, &alpha); // x = z2 * alpha
        if (is_on_curve(&x, curve)) {
            fp2_copy(&P->x, &x);
            fp2_set_one(&P->z);
            break;
        }

        assert(hint < UINT_FAST8_MAX);
    }

    return hint;
}

/// Finds a point of order k * 2^n where n is the largest power of two
/// dividing (p+1) using a hint such that z2 = i + hint above the point
/// (0 : 0).
static void
ec_curve_to_point_2f_above_montgomery_from_hint(ec_point_t *P, ec_curve_t *curve, uint_fast8_t hint)
{
    // Compute a root of x^2 + (A/C)x + 1
    fp2_t alpha;
    compute_alpha(&alpha, curve);

    // Compute the x coordinate from the hint and alpha
    // With 1/2^20 chance we can use the table look up
    fp2_t z1, z2;
    if (hint < LEN_NQR_TABLE) {
        z2=Z_NQR_TABLE[hint];
    }
    // Otherwise we create this using the form i + hint
    else {
        fp_set_one(&z2.im);
        fp_set_small(&z2.re, hint);
    }

    fp2_t x;
    fp2_mul(&x, &z2, &alpha);

    // Set the point
    fp2_copy(&P->x, &x);
    fp2_set_one(&P->z);
}

/// Finds a point of order k * 2^n where n is the largest power of two
/// dividing (p+1).
/// The x-coordinate is picked such that the point (0 : 0) is never the
/// point of order two.
static uint_fast8_t
ec_curve_to_point_2f_not_above_montgomery(ec_point_t *P, const ec_curve_t *curve)
{
    fp2_t x, t0;

    fp_t one;
    fp_set_one(&one);

    uint_fast8_t hint = 0;
    for (;; ++hint) {
        // For each guess of an x, we expect it to be a point 1/2
        // the time, so our table look up will work with failure 2^20
        if (hint < LEN_NQR_TABLE) {
            x=NQR_TABLE[hint];
        }

        // Fallback method in case we do not find a value!
        // For the cases where we are unlucky, we try points of the form
        // x = hint + i
        else {
            // When we first hit this loop, set x to be i + (hint - 1)
            // NOTE: we do hint -1 as we add one before checking a square
            //       in the below loop
            if (hint == LEN_NQR_TABLE) {
                fp_set_one(&x.im);
                fp_set_small(&x.re, hint - 1);
            }

            // Now we find a t which is a NQR of the form i + hint
            for (;; ++hint) {
                // Increase the real part by one until a NQR is found
                // TODO: could be made faster by adding one rather
                // than setting each time, but this is OK for now.
                fp_add(&x.re, &x.re, &one);
                if (!fp2_is_square(&x)){
                    break;
                }

                assert(hint < UINT_FAST8_MAX);
            }
        }

        // Now we have x which is a NQR -- is it on the curve?
        if (is_on_curve(&x, curve)) {
            fp2_copy(&P->x, &x);
            fp2_set_one(&P->z);
            break;
        }
        assert(hint < UINT_FAST8_MAX);
    }

    return hint;
}

/// Finds a point of order k * 2^n where n is the largest power of two
/// dividing (p+1) using a hint such that z2 = i + hint not above
/// the point (0 : 0).
static void
ec_curve_to_point_2f_not_above_montgomery_from_hint(ec_point_t *P,
                                                    const ec_curve_t *curve,
                                                    uint_fast8_t hint)
{
    fp2_t x;

    // If we got lucky (1/2^20) then we just grab an x-value
    // from the table
    if (hint < LEN_NQR_TABLE) {
        x=NQR_TABLE[hint];
    }
    // Otherwise, we find points of the form
    // i + hint
    else {
        fp_set_one(&x.im);
        fp_set_small(&x.re, hint);
    }

    fp2_copy(&P->x, &x);
    fp2_set_one(&P->z);
}

// Helper function which given a point of order k*2^n with n maximal
// and k odd, computes a point of order 2^f
static inline void
clear_cofactor_for_maximal_even_order(ec_point_t *P, ec_curve_t *curve, int f)
{
    // clear out the odd cofactor to get a point of order 2^n
    ec_mul(P, p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, P, curve);

    // clear the power of two to get a point of order 2^f
    for (int i = 0; i < POWER_OF_2 - f; i++) {
        xDBL_A24_normalized(P, P, &curve->A24);
    }
}

// Computes a basis E[2^f] = <P, Q> where the point Q is above (0 : 0)
void
ec_curve_to_basis_2f(ec_basis_t *PQ2, ec_curve_t *curve, int f)
{
    // Normalise (A/C : 1) and ((A + 2)/4 : 1)
    ec_normalize_curve_and_A24(curve);

    // Compute the points P, Q
    ec_point_t P, Q;
    ec_curve_to_point_2f_not_above_montgomery(&P, curve);
    ec_curve_to_point_2f_above_montgomery(&Q, curve);

    // clear out the odd cofactor to get a point of order 2^f
    clear_cofactor_for_maximal_even_order(&P, curve, f);
    clear_cofactor_for_maximal_even_order(&Q, curve, f);

    // Copy points and compute the difference point
    compute_xonly_basis(PQ2, curve, &P, &Q);
}

// Computes a basis E[2^f] = <P, Q> where the point Q is above (0 : 0)
// and stores hints as an array for faster recomputation at a later point
void
ec_curve_to_basis_2f_to_hint(ec_basis_t *PQ2, uint8_t hint[2], ec_curve_t *curve, int f)
{
    // Normalise (A/C : 1) and ((A + 2)/4 : 1)
    ec_normalize_curve_and_A24(curve);

    // Compute the points P, Q
    ec_point_t P, Q;
    hint[0] = ec_curve_to_point_2f_not_above_montgomery(&P, curve);
    hint[1] = ec_curve_to_point_2f_above_montgomery(&Q, curve);

    // clear out the odd cofactor to get a point of order 2^f
    clear_cofactor_for_maximal_even_order(&P, curve, f);
    clear_cofactor_for_maximal_even_order(&Q, curve, f);

    // Copy points and compute the difference point
    compute_xonly_basis(PQ2, curve, &P, &Q);
}

void
ec_curve_to_basis_2f_from_hint(ec_basis_t *PQ2, ec_curve_t *curve, int f, const uint8_t hint[2])
{
    // Normalise (A/C : 1) and ((A + 2)/4 : 1)
    ec_normalize_curve_and_A24(curve);

    // Compute the points P, Q
    ec_point_t P, Q;
    ec_curve_to_point_2f_not_above_montgomery_from_hint(&P, curve, hint[0]);
    ec_curve_to_point_2f_above_montgomery_from_hint(&Q, curve, hint[1]);

    // clear out the odd cofactor to get a point of order 2^f
    clear_cofactor_for_maximal_even_order(&P, curve, f);
    clear_cofactor_for_maximal_even_order(&Q, curve, f);

    // Copy points and compute the difference point
    compute_xonly_basis(PQ2, curve, &P, &Q);
}

// Methods for 3-torsion basis

bool ec_curve_to_3_torsion_point(ec_point_t* P3,ec_point_t* R3,ec_curve_t* curve, const ec_point_t* A3){
    // P3: point of 3-torsion
    // R3: basis point of 3^**-torsion (if found)
    // Returns True if R3 has been found
    // !!! Assumes curve is normalized (C=1) !!!

    ec_point_t P, PTPL;
    fp2_t x;
    uint_fast8_t hint = 0;

    for(;;++hint){
        assert(hint < UINT_FAST8_MAX);
        // Finding a point of the form x=i+hint, z=1
        for(;;++hint){
            assert(hint < UINT_FAST8_MAX);
            fp_set_one(&x.im);
            fp_set_small(&x.re, hint);

            if (is_on_curve(&x, curve)) {
                fp2_copy(&P.x, &x);
                fp2_set_one(&P.z);
                break;
            }
        }

        ec_mul(&P, p_cofactor_for_3g, P_COFACTOR_FOR_3G_BITLENGTH, &P, curve);
        if(ec_is_zero(&P)){
            continue;
        }

        fp2_copy(&(R3->x),&(P.x));
        fp2_copy(&(R3->z),&(P.z));

        int i=0;
        for(;;i++){
            xTPL(&PTPL, &P, A3);
            if(ec_is_zero(&PTPL)){
                break;
            }
            else{
                fp2_copy(&(P.x),&(PTPL.x));
                fp2_copy(&(P.z),&(PTPL.z));
            }
        }
        
        fp2_copy(&(P3->x),&(P.x));
        fp2_copy(&(P3->z),&(P.z));
        if(i==POWER_OF_3){
            return 1;
        }
        else{
            return 0;
        }
    }
}

void ec_curve_tangent_at_point(fp2_t *lambda, fp2_t *mu, const ec_point_t *P, const ec_curve_t * curve){
    // !!! Assumes curve is normalized (C=1) !!!
    
    // Compute y_P
    fp2_t y;
    ec_point_t P1;
    copy_point(&P1,P);
    normalize_point(&P1);
    y_coordinate(&y,&(P1.x),curve);

    // Compute lambda
    fp2_t t0, t1, t2, dbl_y_inv;

    fp2_add(&dbl_y_inv,&y,&y); //2y
    fp2_inv(&dbl_y_inv); //1/(2y)

    fp2_set_one(&t0); //1
    fp2_sqr(&t2,&(P1.x)); //x^2
    fp2_sub(mu,&t0,&t2); //1-x^2
    fp2_mul(mu,mu,&(P1.x)); //x-x^3
    fp2_mul(mu,mu,&dbl_y_inv); //mu=(x-x^3)/(2y)

    fp2_add(&t1,&t2,&t2); //2x^2
    fp2_add(&t2,&t1,&t2); //3x^2
    fp2_mul(&t1,&(curve->A),&(P1.x)); //Ax
    fp2_add(&t1,&t1,&t1); //2Ax
    fp2_add(lambda,&t0,&t1); //2Ax+1
    fp2_add(lambda,lambda,&t2); //3x^2+2Ax+1
    fp2_mul(lambda,lambda,&dbl_y_inv); //lambda=(3x^2+2Ax+1)/(2y)
}

void ec_curve_to_basis_3(ec_basis_t* PQ3, const ec_curve_t* curve){

    fp2_t x, lambda, mu, y, t0;
    ec_point_t P, Q, Q3, P3, R3, S3, A3;
    ec_curve_t E;
    bool found_P;

    // Normalizing the curve
    copy_curve(&E,curve);
    ec_normalize_curve_and_A24(&E);

    // Curve coefficient in the form A3 = (A+2C:A-2C)
    fp2_sub(&A3.z, &(E.A24.x), &(E.A24.z));
    fp2_copy(&A3.x, &(E.A24.x));

    // Computing P3 of 3-torsion and maybe the first basis point P
    found_P=ec_curve_to_3_torsion_point(&P3,&P,&E,&A3);
    ec_curve_tangent_at_point(&lambda, &mu, &P3, &E);

    //Computing P if needed
    uint_fast8_t hint = 0;

    if(found_P==0){
        for(;;++hint){
            assert(hint < UINT_FAST8_MAX);
            // Finding a point of the form x=i+hint, z=1
            for(;;++hint){
                assert(hint < UINT_FAST8_MAX);
                fp_set_one(&x.im);
                fp_set_small(&x.re, hint);

                if (is_on_curve(&x, &E)) {
                    fp2_copy(&P.x, &x);
                    fp2_set_one(&P.z);
                    break;
                }
            }

            ec_mul(&P, p_cofactor_for_3g, P_COFACTOR_FOR_3G_BITLENGTH, &P, &E);
        
            // Testing if P\in [3]E
            normalize_point(&P);
            y_coordinate(&y,&(P.x),&E);
            fp2_mul(&t0,&lambda,&(P.x)); //lambda*x
            fp2_add(&t0,&t0,&mu); //lambda*x+mu
            fp2_sub(&t0,&y,&t0); //y-(lambda*x+mu)

            if(!fp2_is_cube(&t0)){
                break;
            }
        }

        // R3=3^(POWER_OF_3-1)*P
        fp2_copy(&R3.x, &P.x);
        fp2_copy(&R3.z, &P.z);
        for(int i=0;i<POWER_OF_3-1;i++){
            xTPL(&R3,&R3,&A3);
        }
    }
    else{
        // R3=3^(POWER_OF_3-1)*P=P3 (found_P==1)
        fp2_copy(&R3.x, &P3.x);
        fp2_copy(&R3.z, &P3.z);
    }

    //Computing Q
    hint = 0;
    for(;;++hint){
        assert(hint < UINT_FAST8_MAX);
        // Finding a point of the form x=i+hint, z=1
        for(;;++hint){
            assert(hint < UINT_FAST8_MAX);
            fp_set_one(&x.im);
            fp_set_small(&x.re, hint);

            if (is_on_curve(&x, &E)) {
                fp2_copy(&Q.x, &x);
                fp2_set_one(&Q.z);
                break;
            }
        }

        ec_mul(&Q, p_cofactor_for_3g, P_COFACTOR_FOR_3G_BITLENGTH, &Q, &E);
        
        // Testing if Q\in [3]E
        normalize_point(&Q);
        y_coordinate(&y,&(Q.x),&E);
        fp2_mul(&t0,&lambda,&(Q.x)); //lambda*x
        fp2_add(&t0,&t0,&mu); //lambda*x+mu
        fp2_sub(&t0,&y,&t0); //y-(lambda*x+mu)

        if(fp2_is_cube(&t0)){
            continue;
        }

        // Testing if P and Q are independent
        // S3=3^(POWER_OF_3-1)*Q
        fp2_copy(&S3.x, &Q.x);
        fp2_copy(&S3.z, &Q.z);
        for(int i=0;i<POWER_OF_3-1;i++){
            xTPL(&S3,&S3,&A3);
        }

        if(!ec_is_equal(&R3, &S3)){
            break;
        }
    }

    // Compute P-Q
    // TODO: reuse computed sqrts for P and Q
    difference_point(&PQ3->PmQ, &P, &Q, &E);
    copy_point(&PQ3->P, &P);
    copy_point(&PQ3->Q, &Q);
}
