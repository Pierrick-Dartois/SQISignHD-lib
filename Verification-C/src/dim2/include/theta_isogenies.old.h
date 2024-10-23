#include <sqisign_namespace.h>
#include <ec.h>
#include <fp2.h>
#include "theta_structure.h"
#include <hd.h>
#include <hd_splitting_transforms.h>

/** @file
 *
 * @authors Antonin Leroux
 *
 * @brief the theta isogeny header
 */

#ifndef THETA_ISOGENY_H
#define THETA_ISOGENY_H

/*************************** Functions *****************************/

/**
 * @brief Evaluate a gluing isogeny from an elliptic product on a basis
 *
 * @param image1 Output: the theta_point of the image of the first couple of points
 * @param image2 Output : the theta point of the image of the second couple of points
 * @param P : a couple point in E12
 * @param Q : a couple point in E12
 * @param PmQ : a couple point in E12 corresponding to P-Q
 * @param a : ibz
 * @param b : ibz
 * @param E12 : an elliptic product
 * @param phi : a gluing isogeny E1 x E2 -> A
 *
 * The integers a,b are such that the phi.K1_4 = a P + b Q
 * out : phi( P ), Phi (Q)
 *
 */
void gluing_eval_basis(theta_point_t *image1,
                       theta_point_t *image2,
                       const theta_couple_jac_point_t *xyT1,
                       const theta_couple_jac_point_t *xyT2,
                       const theta_couple_curve_t *E12,
                       const theta_gluing_t *phi);

/**
 * @brief Compute  a (2,2) isogeny in dimension 2 in the theta_model
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_8 a point in A[8]
 * @param T2_8 a point in A[8]
 * @param hadamard_bool_1 a boolean
 * @param hadamard_bool_2 a boolean
 *
 * out : A -> B of kernel [4](T1_8,T2_8)
 * hadamard_bool_1 controls if the domain is in standard or dual coordinates
 * hadamard_bool_2 controls if the codomain is in standard or dual coordinates
 *
 */
void theta_isogeny_compute(theta_isogeny_t *out,
                           const theta_structure_t *A,
                           const theta_point_t *T1_8,
                           const theta_point_t *T2_8,
                           bool hadamard_bool_1,
                           bool hadamard_bool_2);

/**
 * @brief Evaluate a theta isogeny
 *
 * @param out Output: the evaluating point
 * @param phi a theta isogeny
 * @param P a point in the domain of phi
 *
 * out = phi(P)
 *
 */
void theta_isogeny_eval(theta_point_t *out, const theta_isogeny_t *phi, const theta_point_t *P);

/**
 * @brief Compute the splitting isomorphism from a theta structure to the product theta structure,
 * returns false if the given theta structure is not isomorphic to an elliptic product
 *
 * @return a boolean indicating if A is isomorphic to an elliptic product
 * @param out: the splitting isomorphism
 * @param A : the theta_structure in consideration
 *
 * out : A -> B where B is theta product associated to ExF an elliptic product
 */
int splitting_compute(theta_splitting_t *out, const theta_structure_t *A);

#endif
