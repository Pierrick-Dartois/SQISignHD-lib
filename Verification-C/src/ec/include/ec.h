/** @file
 *
 * @authors Luca De Feo, Francisco RH
 *
 * @brief Elliptic curve stuff
 */

#ifndef EC_H
#define EC_H
//#include <sqisign_namespace.h>
//#include <ec_params.h>
#include <fp2.h>
#include <tools.h>
#include <stdio.h>

/** @defgroup ec Elliptic curves
 * @{
 */

/** @defgroup ec_t Data structures
 * @{
 */

/** @brief Projective point
 *
 * @typedef ec_point_t
 *
 * @struct ec_point_t
 *
 * A projective point in (X:Z) or (X:Y:Z) coordinates (tbd).
 */
typedef struct ec_point_t
{
    fp2_t x;
    fp2_t z;
} ec_point_t;

/** @brief Projective point
 *
 * @typedef jac_point_t
 *
 * @struct jac_point_t
 *
 * A projective point in(X:Y:Z) coordinates
 */
typedef struct jac_point_t
{
    fp2_t x;
    fp2_t y;
    fp2_t z;
} jac_point_t;

/** @brief A basis of a torsion subgroup
 *
 * @typedef ec_basis_t
 *
 * @struct ec_basis_t
 *
 * A pair of points (or a triplet, tbd) forming a basis of a torsion subgroup.
 */
typedef struct ec_basis_t
{
    ec_point_t P;
    ec_point_t Q;
    ec_point_t PmQ;
} ec_basis_t;

/** @brief An elliptic curve
 *
 * @typedef ec_curve_t
 *
 * @struct ec_curve_t
 *
 * An elliptic curve in projective Montgomery form
 */
typedef struct ec_curve_t
{
    fp2_t A;
    fp2_t C;                             ///< cannot be 0
    ec_point_t A24;                      // the point (A+2 : 4C)
    bool is_A24_computed_and_normalized; // says if A24 has been computed and normalized
} ec_curve_t;

/** @brief An isogeny of degree a power of 2
 *
 * @typedef ec_isog_even_t
 *
 * @struct ec_isog_even_t
 */
typedef struct ec_isog_even_t
{
    ec_curve_t curve;      ///< The domain curve
    ec_point_t kernel;     ///< A kernel generator
    unsigned short length; ///< The length as a 2-isogeny walk
} ec_isog_even_t;

/** @brief Isomorphism of Montgomery curves
 *
 * @typedef ec_isom_t
 *
 * @struct ec_isom_t
 *
 * The isomorphism is given by the map maps (X:Z) ↦ ( (Nx X - Nz Z) : (D Z) )
 */
typedef struct ec_isom_t
{
    fp2_t Nx;
    fp2_t Nz;
    fp2_t D;
} ec_isom_t;

// end ec_t
/** @}
 */

/** @defgroup ec_curve_t Curves and isomorphisms
 * @{
 */

// Initalisation for curves and points
void ec_curve_init(ec_curve_t *E);
void ec_point_init(ec_point_t *P);

// Copying points, bases and curves
static inline void
copy_point(ec_point_t *P, const ec_point_t *Q)
{
    fp2_copy(&P->x, &Q->x);
    fp2_copy(&P->z, &Q->z);
}

static inline void
copy_basis(ec_basis_t *B1, const ec_basis_t *B0)
{
    copy_point(&B1->P, &B0->P);
    copy_point(&B1->Q, &B0->Q);
    copy_point(&B1->PmQ, &B0->PmQ);
}

void copy_curve(ec_curve_t *E1, const ec_curve_t *E2);

// Functions for working with the A24 point and normalisation

/**
 * @brief Reduce (A : C) to (A/C : 1) in place
 *
 * @param E a curve
 */
void ec_normalize_curve(ec_curve_t *E);

/**
 * @brief Reduce (A + 2 : 4C) to ((A+2)/4C : 1) in place
 *
 * @param E a curve
 */
void ec_curve_normalize_A24(ec_curve_t *E);

/**
 * @brief Normalise both (A : C) and (A + 2 : 4C) as above, in place
 *
 * @param E a curve
 */
void ec_normalize_curve_and_A24(ec_curve_t *E);

/**
 * @brief Given a curve E, compute (A+2 : 4C)
 *
 * @param A24 the value (A+2 : 4C) to return into
 * @param E a curve
 */
static inline void
AC_to_A24(ec_point_t *A24, const ec_curve_t *E)
{
    // Maybe we already have this computed
    if (E->is_A24_computed_and_normalized) {
        copy_point(A24, &E->A24);
        return;
    }

    // A24 = (A+2C : 4C)
    fp2_add(&A24->z, &E->C, &E->C);
    fp2_add(&A24->x, &E->A, &A24->z);
    fp2_add(&A24->z, &A24->z, &A24->z);
}

/**
 * @brief Given a curve the point (A+2 : 4C) compute the curve coefficients (A : C)
 *
 * @param E a curve to compute
 * @param A24 the value (A+2 : 4C)
 */
static inline void
A24_to_AC(ec_curve_t *E, const ec_point_t *A24)
{
    // (A:C) = ((A+2C)*2-4C : 4C)
    fp2_add(&E->A, &A24->x, &A24->x);
    fp2_sub(&E->A, &E->A, &A24->z);
    fp2_add(&E->A, &E->A, &E->A);
    fp2_copy(&E->C, &A24->z);
}

/**
 * @brief j-invariant.
 *
 * @param j_inv computed j_invariant
 * @param curve input curve
 */
void ec_j_inv(fp2_t *j_inv, const ec_curve_t *curve);

/**
 * @brief Isomorphism of elliptic curve
 *
 * @param isom computed isomorphism
 * @param from domain curve
 * @param to image curve
 */
void ec_isomorphism(ec_isom_t *isom, const ec_curve_t *from, const ec_curve_t *to);

/**
 * @brief In-place evaluation of an isomorphism
 *
 * @param P a point
 * @param isom an isomorphism
 */
void ec_iso_eval(ec_point_t *P, ec_isom_t *isom);

/** @}
 */
/** @defgroup ec_point_t Point operations
 * @{
 */

/**
 * @brief Point equality
 *
 * @param P a point
 * @param Q a point
 * @return 0xFFFFFFFF if equal, zero otherwise
 */
uint32_t ec_is_equal(const ec_point_t *P, const ec_point_t *Q);

/**
 * @brief Point equality
 *
 * @param P a point
 * @return 0xFFFFFFFF if point at infinity, zero otherwise
 */
uint32_t ec_is_zero(const ec_point_t *P);

/**
 * @brief Reduce Z-coordinate of point in place
 *
 * @param P a point
 */
void ec_normalize_point(ec_point_t *P);

// TODO: these are only used internally in the ec module, they could be kept elsewhere?
void xADD(ec_point_t *R, const ec_point_t *P, const ec_point_t *Q, const ec_point_t *PQ);
void xDBL_A24(ec_point_t *Q, const ec_point_t *P, const ec_point_t *A24);
void xDBL_A24_normalized(ec_point_t *Q, const ec_point_t *P, const ec_point_t *A24);
void xMUL(ec_point_t *Q, const ec_point_t *P, const digit_t *k, const int kbits, const ec_curve_t *curve);
void xMUL_A24(ec_point_t *Q,
                    const ec_point_t *P,
                    const digit_t *k,
                    const int kbits,
                    const ec_point_t *A24);
void xDBLMUL(ec_point_t *S,
        const ec_point_t *P,
        const digit_t *k,
        const ec_point_t *Q,
        const digit_t *l,
        const ec_point_t *PQ,
        const int kbits,
        const ec_curve_t *curve);
void xDBLADD(ec_point_t *R,
        ec_point_t *S,
        const ec_point_t *P,
        const ec_point_t *Q,
        const ec_point_t *PQ,
        const ec_point_t *A24);
void xDBLADD_normalized(ec_point_t *R,
                   ec_point_t *S,
                   const ec_point_t *P,
                   const ec_point_t *Q,
                   const ec_point_t *PQ,
                   const ec_point_t *A24);

void xTPL(ec_point_t* Q, const ec_point_t* P, const ec_point_t* A3);

/**
 * @brief Point doubling
 *
 * @param res computed double of P
 * @param P a point
 * @param curve an elliptic curve
 */
void ec_dbl(ec_point_t *res, const ec_point_t *P, const ec_curve_t *curve);

/**
 * @brief Point iterated doubling
 *
 * @param res computed double of P
 * @param P a point
 * @param n the number of double
 */
void ec_dbl_iter(ec_point_t *res, int n, const ec_point_t *P, ec_curve_t *curve);

/**
 * @brief Iterated doubling for a basis P, Q, PmQ
 *
 * @param res the computed iterated double of basis B
 * @param n the number of doubles
 * @param B the basis to double
 * @param curve the parent curve of the basis
 */
void ec_dbl_iter_basis(ec_basis_t *res, int n, const ec_basis_t *P, ec_curve_t *curve);

/**
 * @brief Point multiplication
 *
 * @param res computed scalar * P
 * @param curve the curve
 * @param scalar an unsigned multi-precision integer
 * @param P a point
 */
void ec_mul(ec_point_t *res,
            const digit_t *scalar,
            const int kbits,
            const ec_point_t *P,
            ec_curve_t *curve);

/**
 * @brief Combination P+m*Q
 *
 * @param R computed P + m * Q
 * @param curve the curve
 * @param m an unsigned multi-precision integer
 * @param P a point
 * @param Q a point
 * @param PQ the difference P-Q
 */
void ec_ladder3pt(ec_point_t *R,
                  digit_t *const m,
                  const ec_point_t *P,
                  const ec_point_t *Q,
                  const ec_point_t *PQ,
                  const ec_curve_t *A);

/**
 * @brief Linear combination of points of a basis
 *
 * @param res computed scalarP * P + scalarQ * Q
 * @param scalarP an unsigned multi-precision integer
 * @param scalarQ an unsigned multi-precision integer
 * @param kbits number of bits of the scalars, or n for points of order 2^n
 * @param PQ a torsion basis consisting of points P and Q
 * @param curve the curve
 */
void ec_biscalar_mul(ec_point_t *res,
                     const digit_t *scalarP,
                     const digit_t *scalarQ,
                     const int kbits,
                     const ec_basis_t *PQ,
                     const ec_curve_t *curve);

/** @defgroup ec_dlog_t Torsion basis computations
 * @{
 */

/**
 * @brief Generate a Montgomery curve and a 2^f-torsion basis
 *
 * The algorithm is deterministc
 */
void ec_curve_to_basis_2f(ec_basis_t *PQ2, ec_curve_t *curve, int f);
void ec_curve_to_basis_2f_to_hint(ec_basis_t *PQ2, uint8_t hint[2], ec_curve_t *curve, int f);
void ec_curve_to_basis_2f_from_hint(ec_basis_t *PQ2,
                                    ec_curve_t *curve,
                                    int f,
                                    const uint8_t hint[2]);

/** 
 * @brief Generate a 3^g-torsion basis 
 **/
void ec_curve_to_basis_3(ec_basis_t* PQ3, const ec_curve_t* curve);

/** @defgroup ec_isog_t Isogenies
 * @{
 */

/**
 * @brief Evaluate isogeny of even degree on list of points
 *
 * @param image computed image curve
 * @param phi isogeny
 * @param points a list of points to evaluate the isogeny on, modified in place
 * @param length of the list points
 */
void ec_eval_even(ec_curve_t *image,
                  ec_isog_even_t *phi,
                  ec_point_t *points,
                  unsigned short length);

// naive implementation for very small 2^n chain
void ec_eval_small_chain(ec_curve_t *image,
                         const ec_point_t *kernel,
                         int len,
                         ec_point_t *points,
                         int len_points);

// Jacobian point init and copying
void jac_init(jac_point_t *P);
void copy_jac_point(jac_point_t *P, const jac_point_t *Q);
uint32_t is_jac_equal(const jac_point_t *P, const jac_point_t *Q);

// Convert from Jacobian to x-only (just drop the Y-coordinate)
void jac_to_xz(ec_point_t *P, const jac_point_t *xyP);

// Jacobian arithmetic
void jac_neg(jac_point_t *Q, const jac_point_t *P);
void ADD(jac_point_t *R, const jac_point_t *P, const jac_point_t *Q, const ec_curve_t *AC);
void DBL(jac_point_t *Q, const jac_point_t *P, const ec_curve_t *AC);

// Given a basis in x-only, lift to a pair of Jacobian points
void lift_basis_normalized(jac_point_t *P, jac_point_t *Q, ec_basis_t *B, ec_curve_t *E);
void lift_basis(jac_point_t *P, jac_point_t *Q, ec_basis_t *B, ec_curve_t *E);

// Given a point in x-only, lift to a Jacobian point
void lift_point_normalized(jac_point_t *P, const ec_point_t *Q, const ec_curve_t *E);
void lift_point(jac_point_t *P, ec_point_t *Q, ec_curve_t *E);

// Check if a non-normalized point is on a non-normalized curve
uint32_t ec_is_on_curve_normalized(const ec_point_t *P, const ec_curve_t *E);
uint32_t ec_is_on_curve(ec_point_t *P, ec_curve_t *E);


/*
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Test functions for printing and order checking, only used in debug mode
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

// Check if a point (X : Z) has order exactly 2^f
static int
test_point_order_twof(const ec_point_t *P, const ec_curve_t *E, int t)
{
    ec_point_t test;
    ec_curve_t curve;
    test = *P;
    copy_curve(&curve, E);

    if (ec_is_zero(&test))
        return 0;
    // Scale point by 2^(t-1)
    ec_dbl_iter(&test, t - 1, &test, &curve);
    // If it's zero now, it doesnt have order 2^t
    if (ec_is_zero(&test))
        return 0;
    // Ensure [2^t] P = 0
    ec_dbl(&test, &test, &curve);
    return ec_is_zero(&test);
}

// Check if a Jacobian point (X : Y : Z) has order exactly 2^f
static int
test_jac_order_twof(const jac_point_t *P, const ec_curve_t *E, int t)
{
    jac_point_t test;
    test = *P;
    if (fp2_is_zero(&test.z))
        return 0;
    for (int i = 0; i < t - 1; i++) {
        DBL(&test, &test, E);
    }
    if (fp2_is_zero(&test.z))
        return 0;
    DBL(&test, &test, E);
    return (fp2_is_zero(&test.z));
}

// Prints the x-coordinate of the point (X : 1)
static void
ec_point_print(char *name, ec_point_t P)
{
    fp2_t a;
    if (fp2_is_zero(&P.z)) {
        printf("%s = INF\n", name);
    } else {
        fp2_copy(&a, &P.z);
        fp2_inv(&a);
        fp2_mul(&a, &a, &P.x);
        fp2_print(name, &a);
    }
}

// Prints the Montgomery coefficient A
static void
ec_curve_print(char *name, ec_curve_t E)
{
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, &a);
}

#endif
