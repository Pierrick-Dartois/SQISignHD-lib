/** @file
 *
 * @authors Antonin Leroux
 *
 * @brief The HD-isogenies algorithm required by the signature
 *
 */

#ifndef HD_H
#define HD_H

#include <sqisign_namespace.h>
#include <ec.h>
#include <torsion_constants.h>
#include <stdio.h>

/*************************** Structures *****************************/

/** @brief Type for couple point with XZ coordinates
 * @typedef theta_couple_point_t
 *
 * @struct theta_couple_point
 *
 * Structure for the couple point on an elliptic product
 * using XZ coordinates
 */
typedef struct theta_couple_point
{
    ec_point_t P1;
    ec_point_t P2;
} theta_couple_point_t;

/** @brief Type for couple point with XYZ coordinates
 * @typedef theta_couple_jac_point_t
 *
 * @struct theta_couple_jac_point
 *
 * Structure for the couple point on an elliptic product
 * using XYZ coordinates
 */
typedef struct theta_couple_jac_point
{
    jac_point_t P1;
    jac_point_t P2;
} theta_couple_jac_point_t;

/** @brief Type for couple curve *
 * @typedef theta_couple_curve_t
 *
 * @struct theta_couple_curve
 *
 * the  theta_couple_curve structure
 */
typedef struct theta_couple_curve
{
    ec_curve_t E1;
    ec_curve_t E2;
} theta_couple_curve_t;

/** @brief Type for theta point *
 * @typedef theta_point_t
 *
 * @struct theta_point
 *
 * the  theta_point structure used
 */
typedef struct theta_point
{
    fp2_t x;
    fp2_t y;
    fp2_t z;
    fp2_t t;
} theta_point_t;

/** @brief Type for theta structure *
 * @typedef theta_structure_t
 *
 * @struct theta_structure
 *
 * the  theta_structure structure used
 */
typedef struct theta_structure
{
    theta_point_t null_point;
    bool precomputation;

    // Eight precomputed values used for doubling and
    // (2,2)-isogenies.
    fp2_t XYZ0;
    fp2_t YZT0;
    fp2_t XZT0;
    fp2_t XYT0;

    fp2_t xyz0;
    fp2_t yzt0;
    fp2_t xzt0;
    fp2_t xyt0;
} theta_structure_t;

/** @brief A 2x2 matrix used for action by translation
 * @typedef translation_matrix_t
 *
 * @struct translation_matrix
 *
 * Structure to hold 4 fp2_t elements representing a 2x2 matrix used when computing
 * a compatible theta structure during gluing.
 */
typedef struct translation_matrix
{
    fp2_t g00;
    fp2_t g01;
    fp2_t g10;
    fp2_t g11;
} translation_matrix_t;

/** @brief A 4x4 matrix used for basis changes
 * @typedef basis_change_matrix_t
 *
 * @struct basis_change_matrix
 *
 * Structure to hold 16 elements representing a 4x4 matrix used for changing
 * the basis of a theta point.
 */
typedef struct basis_change_matrix
{
    fp2_t m00;
    fp2_t m01;
    fp2_t m02;
    fp2_t m03;

    fp2_t m10;
    fp2_t m11;
    fp2_t m12;
    fp2_t m13;

    fp2_t m20;
    fp2_t m21;
    fp2_t m22;
    fp2_t m23;

    fp2_t m30;
    fp2_t m31;
    fp2_t m32;
    fp2_t m33;
} basis_change_matrix_t;

/** @brief Type for gluing (2,2) theta isogeny *
 * @typedef theta_gluing_t
 *
 * @struct theta_gluing
 *
 * the  theta_gluing structure
 */
typedef struct theta_gluing
{

    theta_couple_curve_t domain;
    theta_couple_point_t K1_4;
    theta_couple_point_t K2_4;
    theta_couple_jac_point_t xyK1_4;
    theta_couple_jac_point_t xyK2_4;
    theta_point_t T1_8;
    theta_point_t T2_8;
    fp2_t M[4][4];
    theta_point_t precomputation;
    theta_point_t codomain;
    theta_point_t np_domain;
    theta_point_t prod_null_point;
    int zero_idx;
} theta_gluing_t;

/** @brief Type for standard (2,2) theta isogeny *
 * @typedef theta_isogeny_t
 *
 * @struct theta_isogeny
 *
 * the  theta_isogeny structure
 */
typedef struct theta_isogeny
{
    theta_point_t T1_8;
    theta_point_t T2_8;
    bool hadamard_bool_1;
    bool hadamard_bool_2;
    theta_structure_t domain;
    theta_point_t precomputation;
    theta_structure_t codomain;
} theta_isogeny_t;

/** @brief Type for splitting isomorphism *
 * @typedef theta_splitting_t
 *
 * @struct theta_splitting
 *
 * the theta_splitting structure
 */
typedef struct theta_splitting
{
    fp2_t M[4][4];
    //theta_structure_t B;
    theta_point_t domain;
    theta_couple_curve_t codomain;
    theta_point_t precomputation;
    theta_point_t prod_null_point;
} theta_splitting_t;

/** @brief Type for chain of (2,2) theta isogenies with preimage by 4 multiplication above *
 * @typedef theta_chain_t
 *
 * @struct theta_chain
 *
 * the  theta_chain structure
 */
typedef struct theta_chain
{
    theta_couple_point_t T1;
    theta_couple_point_t T2;
    int length;
    theta_couple_curve_t domain;
    theta_couple_curve_t codomain;
    theta_gluing_t first_step;
    theta_splitting_t last_step;
    theta_isogeny_t *steps;
} theta_chain_t;

static void
theta_chain_init(theta_chain_t *chain, size_t num_steps)
{
    memset(chain, 0, sizeof(*chain));
    chain->steps = calloc(num_steps, sizeof(*chain->steps));
}

static void
theta_chain_finalize(theta_chain_t *chain)
{
    free(chain->steps);
    chain->steps = NULL;
}

/*************************** Functions *****************************/

/**
 * @brief Compute the double of the theta couple point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_point
 * @param in the theta couple point in the elliptic product
 * @param E12 an elliptic product
 * in = (P1,P2)
 * out = [2] (P1,P2)
 *
 */
void double_couple_point(theta_couple_point_t *out,
                         const theta_couple_point_t *in,
                         const theta_couple_curve_t *E1E2);

/**
 * @brief Compute the iterated double of the theta couple point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_point
 * @param n : the number of iteration
 * @param E12 an elliptic product
 * @param in the theta couple point in the elliptic product
 * in = (P1,P2)
 * out = [2^n] (P1,P2)
 *
 */
void double_couple_point_iter(theta_couple_point_t *out,
                              unsigned n,
                              const theta_couple_point_t *in,
                              const theta_couple_curve_t *E1E2);

/**
 * @brief Compute the addition of two points in (X : Y : Z) coordinates on the elliptic product E12
 *
 * @param out Output: the theta_couple_jac_point
 * @param T1 the theta couple jac point in the elliptic product
 * @param T2 the theta couple jac point in the elliptic product
 * @param E1E2 an elliptic product
 * in  = (P1, P2), (Q1, Q2)
 * out = (P1 + Q1, P2 + Q2)
 *
 **/
void add_couple_jac_points(theta_couple_jac_point_t *out,
                           const theta_couple_jac_point_t *T1,
                           const theta_couple_jac_point_t *T2,
                           const theta_couple_curve_t *E1E2);

/**
 * @brief Compute the double of the theta couple point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_point
 * @param in the theta couple point in the elliptic product
 * @param E12 an elliptic product
 * in = (P1,P2)
 * out = [2] (P1,P2)
 *
 */
void double_couple_jac_point(theta_couple_jac_point_t *out,
                             const theta_couple_jac_point_t *in,
                             const theta_couple_curve_t *E1E2);

/**
 * @brief Compute the iterated double of the theta couple jac point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_jac_point
 * @param n : the number of iteration
 * @param in the theta couple jac point in the elliptic product
 * @param E12 an elliptic product
 * in  = (P1,P2)
 * out = [2^n] (P1,P2)
 *
 */
void double_couple_jac_point_iter(theta_couple_jac_point_t *out,
                                  unsigned n,
                                  const theta_couple_jac_point_t *in,
                                  const theta_couple_curve_t *E1E2);

/**
 * @brief A forgetful function which returns (X : Z) points given a pair of (X : Y : Z) points
 *
 * @param P Output: the theta_couple_point
 * @param xyP : the theta_couple_jac_point
 **/
void couple_jac_to_xz(theta_couple_point_t *P, const theta_couple_jac_point_t *xyP);

/**
 * @brief Compute  a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 * Returns 0 if the codomain fails to split and 1 otherwise.
 * @param out Output: the theta_chain
 * @param n : the length of the isogeny chain
 * @param E12 an elliptic curve product
 * @param T1 a couple point on E12[2^(n+2)]
 * @param T2 a couple point on E12[2^(n+2)]
 * @param T1m2 a couple point on E12[2^(n+2)] equal to T1-T2
 * @param eight_above boolean indicating if we give the points in E12[2^n] or E12[2^(n+2)]
 *
 * out : E1xE2 -> E3xE4 of kernel [4](T1,T2)
 * uses the precomputed strategy for this length
 *
 */
int theta_chain_compute(theta_chain_t *out,
                        unsigned n,
                        theta_couple_curve_t *E12,
                        const theta_couple_point_t *T1,
                        const theta_couple_point_t *T2,
                        const theta_couple_point_t *T1m2,
                        bool eight_above);

/**
 * @brief Compute a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 * and evaluate at a list of points of the form (P1,0) or (0,P2). Returns 0 if the codomain fails to
 * split and 1 otherwise.
 *
 * @param out Output: the theta_chain
 * @param n : the length of the isogeny chain
 * @param E12 an elliptic curve product
 * @param T1 a couple point on E12[2^(n+2)]
 * @param T2 a couple point on E12[2^(n+2)]
 * @param T1m2 a couple point on E12[2^(n+2)] equal to T1-T2
 * @param eight_above boolean indicating if we give the points in E12[2^n] or E12[2^(n+2)]
 * @param E34 Output: the codomain curve
 * @param P12 Input/Output: pointer to points to be pushed through the isogeny (in-place)
 * @param numP: length of the list of points given in P12 (can be zero)
 *
 * uses the precomputed strategy for this length
 *
 */

int theta_chain_compute_and_eval(unsigned n,
                                 /*const*/ theta_couple_curve_t *E12,
                                 const theta_couple_point_t *T1,
                                 const theta_couple_point_t *T2,
                                 const theta_couple_point_t *T1m2,
                                 bool eight_above,
                                 theta_couple_curve_t *E34,
                                 theta_couple_point_t *P12,
                                 size_t numP);

/**
 * @brief Evaluate a (2,2) isogeny chain in dimension 2 between elliptic products in the theta_model
 *
 * @param out Output: the image point
 * @param phi : the (2,2) isogeny chain of domain E12
 * @param P12 a couple point on E12 where one of the two point is zero,
 *
 * phi : E1xE2 -> E3xE4 of kernel
 * P12 in E1xE2 with P12 = (P1,0) or (0,P2)
 * out = phi(P12) in E3xE4
 *
 */
void theta_chain_eval_special_case(theta_couple_point_t *out,
                                   const theta_chain_t *phi,
                                   const theta_couple_point_t *P12);

#endif
