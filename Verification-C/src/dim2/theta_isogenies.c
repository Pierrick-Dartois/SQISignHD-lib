#include <theta_isogenies.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <tools.h>

static inline void
choose_index_theta_point(fp2_t *res, int ind, const theta_point_t *T)
{
    const fp2_t *src;
    switch (ind % 4) {
        case 0:
            src = &T->x;
            break;
        case 1:
            src = &T->y;
            break;
        case 2:
            src = &T->z;
            break;
        case 3:
            src = &T->t;
            break;
        default:
            assert(0);
    }
    fp2_copy(res, src);
}

void
apply_isomorphism(theta_point_t *res, const basis_change_matrix_t *M, const theta_point_t *P)
{
    fp2_t x1;
    theta_point_t temp;

    fp2_mul(&temp.x, &P->x, &M->m00);
    fp2_mul(&x1, &P->y, &M->m01);
    fp2_add(&temp.x, &temp.x, &x1);
    fp2_mul(&x1, &P->z, &M->m02);
    fp2_add(&temp.x, &temp.x, &x1);
    fp2_mul(&x1, &P->t, &M->m03);
    fp2_add(&temp.x, &temp.x, &x1);

    fp2_mul(&temp.y, &P->x, &M->m10);
    fp2_mul(&x1, &P->y, &M->m11);
    fp2_add(&temp.y, &temp.y, &x1);
    fp2_mul(&x1, &P->z, &M->m12);
    fp2_add(&temp.y, &temp.y, &x1);
    fp2_mul(&x1, &P->t, &M->m13);
    fp2_add(&temp.y, &temp.y, &x1);

    fp2_mul(&temp.z, &P->x, &M->m20);
    fp2_mul(&x1, &P->y, &M->m21);
    fp2_add(&temp.z, &temp.z, &x1);
    fp2_mul(&x1, &P->z, &M->m22);
    fp2_add(&temp.z, &temp.z, &x1);
    fp2_mul(&x1, &P->t, &M->m23);
    fp2_add(&temp.z, &temp.z, &x1);

    fp2_mul(&temp.t, &P->x, &M->m30);
    fp2_mul(&x1, &P->y, &M->m31);
    fp2_add(&temp.t, &temp.t, &x1);
    fp2_mul(&x1, &P->z, &M->m32);
    fp2_add(&temp.t, &temp.t, &x1);
    fp2_mul(&x1, &P->t, &M->m33);
    fp2_add(&temp.t, &temp.t, &x1);

    fp2_copy(&res->x, &temp.x);
    fp2_copy(&res->y, &temp.y);
    fp2_copy(&res->z, &temp.z);
    fp2_copy(&res->t, &temp.t);
}

// compute the theta_point corresponding to the couple of point T on an elliptic product
void
base_change(theta_point_t *out, const theta_gluing_t *phi, const theta_couple_point_t *T)
{
    uint32_t c1, c2;
    theta_point_t null_point;
    fp2_t x1, x2, one;
    fp2_set_one(&one);

    // If P1 = (0 : 0) set x1 = 1 else
    // x1 = x(P1)
    c1 = fp2_is_zero(&T->P1.x);
    c2 = fp2_is_zero(&T->P1.z);
    fp2_select(&x1, &T->P1.x, &one, c1 & c2);

    // If P1 = (0 : 0) set x2 = 1 else
    // x2 = x(P2)
    c1 = fp2_is_zero(&T->P2.x);
    c2 = fp2_is_zero(&T->P2.z);
    fp2_select(&x2, &T->P2.x, &one, c1 & c2);

    // null_point = (a : b : c : d)
    // a = P1.x P2.x, b = P1.x P2.z, c = P1.z P2.x, d = P1.z P2.z
    fp2_mul(&null_point.x, &x1, &x2);
    fp2_mul(&null_point.y, &x1, &T->P2.z);
    fp2_mul(&null_point.z, &x2, &T->P1.z);
    fp2_mul(&null_point.t, &T->P1.z, &T->P2.z);

    // Apply the basis change
    apply_isomorphism(out, &phi->M, &null_point);
}

static void
action_by_translation_z_and_det(fp2_t *z_inv,
                                fp2_t *det_inv,
                                const ec_point_t *P4,
                                const ec_point_t *P2)
{
    // Store the Z-coordinate to invert
    fp2_copy(z_inv, &P4->z);

    // Then collect detij = xij wij - uij zij
    fp2_t tmp;
    fp2_mul(det_inv, &P4->x, &P2->z);
    fp2_mul(&tmp, &P4->z, &P2->x);
    fp2_sub(det_inv, det_inv, &tmp);
}

static void
action_by_translation_compute_matrix(translation_matrix_t *G,
                                     const ec_point_t *P4,
                                     const ec_point_t *P2,
                                     const fp2_t *z_inv,
                                     const fp2_t *det_inv)
{
    fp2_t tmp;

    // Gi.g10 = uij xij /detij - xij/zij
    fp2_mul(&tmp, &P4->x, z_inv);
    fp2_mul(&G->g10, &P4->x, &P2->x);
    fp2_mul(&G->g10, &G->g10, det_inv);
    fp2_sub(&G->g10, &G->g10, &tmp);

    // Gi.g11 = uij zij * detij
    fp2_mul(&G->g11, &P2->x, det_inv);
    fp2_mul(&G->g11, &G->g11, &P4->z);

    // Gi.g00 = -Gi.g11
    fp2_neg(&G->g00, &G->g11);

    // Gi.g01 = - wij zij detij
    fp2_mul(&G->g01, &P2->z, det_inv);
    fp2_mul(&G->g01, &G->g01, &P4->z);
    fp2_neg(&G->g01, &G->g01);
}

// Returns 1 if the basis is as expected and 0 otherwise
// We only expect this to fail for malformed signatures, so
// do not require this to run in constant time.
static int
verify_two_torsion(const theta_couple_point_t *K1_2,
                   const theta_couple_point_t *K2_2,
                   const theta_couple_curve_t *E12)
{
    // First check if any point in K1_2 or K2_2 is zero, if they are then the points did not have
    // order 8 when we started gluing
    if (ec_is_zero(&K1_2->P1) || ec_is_zero(&K1_2->P2) || ec_is_zero(&K2_2->P1) ||
        ec_is_zero(&K2_2->P2)) {
        return 0;
    }

    // Now ensure that P1, Q1 and P2, Q2 are independent. For points of order two this means
    // that they're not the same
    if (ec_is_equal(&K1_2->P1, &K2_2->P1) || ec_is_equal(&K1_2->P2, &K2_2->P2)) {
        return 0;
    }

    // Finally, double points to ensure all points have order exactly 0
    theta_couple_point_t O1, O2;
    double_couple_point(&O1, K1_2, E12);
    double_couple_point(&O2, K2_2, E12);
    // If this check fails then the points had order 2*f for some f, and the kernel is malformed.
    if (!(ec_is_zero(&O1.P1) && ec_is_zero(&O1.P2) && ec_is_zero(&O2.P1) && ec_is_zero(&O2.P2))) {
        return 0;
    }

    return 1;
}

// Computes the action by translation for four points
// (P1, P2) and (Q1, Q2) on E1 x E2 simultaneously to
// save on inversions.
// Returns 0 if any of Pi or Qi does not have order 2
// and 1 otherwise
static int
action_by_translation(translation_matrix_t *Gi,
                      const theta_couple_point_t *K1_4,
                      const theta_couple_point_t *K2_4,
                      const theta_couple_curve_t *E12)
{
    // Compute points of order 2 from Ki_4
    theta_couple_point_t K1_2, K2_2;
    double_couple_point(&K1_2, K1_4, E12);
    double_couple_point(&K2_2, K2_4, E12);

    if (!verify_two_torsion(&K1_2, &K2_2, E12)) {
        return 0;
    }

    // We need to invert four Z coordinates and
    // four determinants which we do with batched
    // inversion
    fp2_t inverses[8];
    action_by_translation_z_and_det(&inverses[0], &inverses[4], &K1_4->P1, &K1_2.P1);
    action_by_translation_z_and_det(&inverses[1], &inverses[5], &K1_4->P2, &K1_2.P2);
    action_by_translation_z_and_det(&inverses[2], &inverses[6], &K2_4->P1, &K2_2.P1);
    action_by_translation_z_and_det(&inverses[3], &inverses[7], &K2_4->P2, &K2_2.P2);

    fp2_batched_inv(inverses, 8);

    action_by_translation_compute_matrix(&Gi[0], &K1_4->P1, &K1_2.P1, &inverses[0], &inverses[4]);
    action_by_translation_compute_matrix(&Gi[1], &K1_4->P2, &K1_2.P2, &inverses[1], &inverses[5]);
    action_by_translation_compute_matrix(&Gi[2], &K2_4->P1, &K2_2.P1, &inverses[2], &inverses[6]);
    action_by_translation_compute_matrix(&Gi[3], &K2_4->P2, &K2_2.P2, &inverses[3], &inverses[7]);

    return 1;
}

// Given the appropriate four torsion, computes the
// change of basis to compute the correct theta null
// point.
// Returns 0 if the order of K1_4 or K2_4 is not 4
static int
gluing_change_of_basis(basis_change_matrix_t *M,
                       const theta_couple_point_t *K1_4,
                       const theta_couple_point_t *K2_4,
                       const theta_couple_curve_t *E12)
{
    // Compute the four 2x2 matrices for the action by translation
    // on the four points:
    translation_matrix_t Gi[4];
    if (!action_by_translation(Gi, K1_4, K2_4, E12))
        return 0;

    // Computation of the 4x4 matrix from Mij
    // t001, t101 (resp t002, t102) first column of M11 * M21 (resp M12 * M22)
    fp2_t t001, t101, t002, t102, tmp;

    fp2_mul(&t001, &Gi[0].g00, &Gi[2].g00);
    fp2_mul(&tmp, &Gi[0].g01, &Gi[2].g10);
    fp2_add(&t001, &t001, &tmp);

    fp2_mul(&t101, &Gi[0].g10, &Gi[2].g00);
    fp2_mul(&tmp, &Gi[0].g11, &Gi[2].g10);
    fp2_add(&t101, &t101, &tmp);

    fp2_mul(&t002, &Gi[1].g00, &Gi[3].g00);
    fp2_mul(&tmp, &Gi[1].g01, &Gi[3].g10);
    fp2_add(&t002, &t002, &tmp);

    fp2_mul(&t102, &Gi[1].g10, &Gi[3].g00);
    fp2_mul(&tmp, &Gi[1].g11, &Gi[3].g10);
    fp2_add(&t102, &t102, &tmp);

    // trace for the first row
    fp2_set_one(&M->m00);
    fp2_mul(&tmp, &t001, &t002);
    fp2_add(&M->m00, &M->m00, &tmp);
    fp2_mul(&tmp, &Gi[2].g00, &Gi[3].g00);
    fp2_add(&M->m00, &M->m00, &tmp);
    fp2_mul(&tmp, &Gi[0].g00, &Gi[1].g00);
    fp2_add(&M->m00, &M->m00, &tmp);

    fp2_mul(&M->m01, &t001, &t102);
    fp2_mul(&tmp, &Gi[2].g00, &Gi[3].g10);
    fp2_add(&M->m01, &M->m01, &tmp);
    fp2_mul(&tmp, &Gi[0].g00, &Gi[1].g10);
    fp2_add(&M->m01, &M->m01, &tmp);

    fp2_mul(&M->m02, &t101, &t002);
    fp2_mul(&tmp, &Gi[2].g10, &Gi[3].g00);
    fp2_add(&M->m02, &M->m02, &tmp);
    fp2_mul(&tmp, &Gi[0].g10, &Gi[1].g00);
    fp2_add(&M->m02, &M->m02, &tmp);

    fp2_mul(&M->m03, &t101, &t102);
    fp2_mul(&tmp, &Gi[2].g10, &Gi[3].g10);
    fp2_add(&M->m03, &M->m03, &tmp);
    fp2_mul(&tmp, &Gi[0].g10, &Gi[1].g10);
    fp2_add(&M->m03, &M->m03, &tmp);

    // Compute the action of (0,out.K2_4.P2) for the second row
    fp2_mul(&tmp, &Gi[3].g01, &M->m01);
    fp2_mul(&M->m10, &Gi[3].g00, &M->m00);
    fp2_add(&M->m10, &M->m10, &tmp);

    fp2_mul(&tmp, &Gi[3].g11, &M->m01);
    fp2_mul(&M->m11, &Gi[3].g10, &M->m00);
    fp2_add(&M->m11, &M->m11, &tmp);

    fp2_mul(&tmp, &Gi[3].g01, &M->m03);
    fp2_mul(&M->m12, &Gi[3].g00, &M->m02);
    fp2_add(&M->m12, &M->m12, &tmp);

    fp2_mul(&tmp, &Gi[3].g11, &M->m03);
    fp2_mul(&M->m13, &Gi[3].g10, &M->m02);
    fp2_add(&M->m13, &M->m13, &tmp);

    // compute the action of (K1_4.P1,0) for the third row
    fp2_mul(&tmp, &Gi[0].g01, &M->m02);
    fp2_mul(&M->m20, &Gi[0].g00, &M->m00);
    fp2_add(&M->m20, &M->m20, &tmp);

    fp2_mul(&tmp, &Gi[0].g01, &M->m03);
    fp2_mul(&M->m21, &Gi[0].g00, &M->m01);
    fp2_add(&M->m21, &M->m21, &tmp);

    fp2_mul(&tmp, &Gi[0].g11, &M->m02);
    fp2_mul(&M->m22, &Gi[0].g10, &M->m00);
    fp2_add(&M->m22, &M->m22, &tmp);

    fp2_mul(&tmp, &Gi[0].g11, &M->m03);
    fp2_mul(&M->m23, &Gi[0].g10, &M->m01);
    fp2_add(&M->m23, &M->m23, &tmp);

    // compute the action of (K1_4.P1,K2_4.P2) for the final row
    fp2_mul(&tmp, &Gi[0].g01, &M->m12);
    fp2_mul(&M->m30, &Gi[0].g00, &M->m10);
    fp2_add(&M->m30, &M->m30, &tmp);

    fp2_mul(&tmp, &Gi[0].g01, &M->m13);
    fp2_mul(&M->m31, &Gi[0].g00, &M->m11);
    fp2_add(&M->m31, &M->m31, &tmp);

    fp2_mul(&tmp, &Gi[0].g11, &M->m12);
    fp2_mul(&M->m32, &Gi[0].g10, &M->m10);
    fp2_add(&M->m32, &M->m32, &tmp);

    fp2_mul(&tmp, &Gi[0].g11, &M->m13);
    fp2_mul(&M->m33, &Gi[0].g10, &M->m11);
    fp2_add(&M->m33, &M->m33, &tmp);

    return 1;
}

/**
 * @brief Compute the gluing isogeny from an elliptic product
 *
 * @param out Output: the theta_gluing
 * @param K1_8 a couple point
 * @param E12 an elliptic curve product
 * @param K2_8 a point in E2[8]
 *
 * out : E1xE2 -> A of kernel [4](K1_8,K2_8)
 * if the kernel supplied has the incorrect order, or gluing seems malformed,
 * returns 0, otherwise returns 1.
 */
static int
gluing_compute(theta_gluing_t *out,
               theta_couple_curve_t *E12,
               const theta_couple_jac_point_t *xyK1_8,
               const theta_couple_jac_point_t *xyK2_8)
{
    // Ensure that we have been given the eight torsion
#ifndef NDEBUG
    {
        int check = test_jac_order_twof(&xyK1_8->P1, &E12->E1, 3);
        if (!check)
            debug_print("xyK1_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK2_8->P1, &E12->E1, 3);
        if (!check)
            debug_print("xyK2_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK1_8->P2, &E12->E2, 3);
        if (!check)
            debug_print("xyK2_8->P1 does not have order 8");
        check = test_jac_order_twof(&xyK2_8->P2, &E12->E2, 3);
        if (!check)
            debug_print("xyK2_8->P2 does not have order 8");
    }
#endif
    // Given points in E[8] x E[8] we need the four torsion below
    double_couple_jac_point(&out->xyK1_4, xyK1_8, E12);
    double_couple_jac_point(&out->xyK2_4, xyK2_8, E12);

    // Convert from (X:Y:Z) coordinates to (X:Z)
    theta_couple_point_t K1_8, K2_8;

    couple_jac_to_xz(&K1_8, xyK1_8);
    couple_jac_to_xz(&K2_8, xyK2_8);
    couple_jac_to_xz(&out->K1_4, &out->xyK1_4);
    couple_jac_to_xz(&out->K2_4, &out->xyK2_4);

    // Set the basis change matrix, if we have not been given a valid K[8] for this computation
    // gluing_change_of_basis will detect this and return 0
    if (!gluing_change_of_basis(&out->M, &out->K1_4, &out->K2_4, E12)) {
        debug_print("gluing failed as kernel does not have correct order");
        return 0;
    }

    // apply the base change to the kernel
    base_change(&out->T1_8, out, &K1_8);
    base_change(&out->T2_8, out, &K2_8);

    // compute the codomain
    theta_point_t TT1, TT2;
    to_squared_theta(&TT1, &out->T1_8);
    to_squared_theta(&TT2, &out->T2_8);

    // If the kernel is well formed then TT1.t and TT2.t are zero
    // if they are not, we exit early as the signature we are validating
    // is probably malformed
    if (!(fp2_is_zero(&TT1.t) || fp2_is_zero(&TT2.t))) {
        debug_print("gluing failed TT1.t or TT2.t is not zero");
        return 0;
    }

    fp2_t t1, t2, t3, t4;
    fp2_copy(&t1, &TT2.z);
    fp2_copy(&t2, &TT1.y);
    fp2_copy(&t3, &TT2.x);
    fp2_copy(&t4, &TT1.x);

    // Invert all values
    fp2_t t_inv[4];
    fp2_copy(&t_inv[0], &t1);
    fp2_copy(&t_inv[1], &t2);
    fp2_copy(&t_inv[2], &t3);
    fp2_copy(&t_inv[3], &t4);
    fp2_batched_inv(t_inv, 4);

    // Compute A, B, C, D
    // (1 : t2 / t4 : t1 / t3 : 0)
    fp2_t tmp;
    fp2_set_one(&out->codomain.x);
    fp2_mul(&tmp, &t2, &t_inv[3]);
    fp2_copy(&out->codomain.y, &tmp);
    fp2_mul(&tmp, &t1, &t_inv[2]);
    fp2_copy(&out->codomain.z, &tmp);
    fp2_set_zero(&out->codomain.t);

    // compute precomp
    // (1 : t4 / t2 : t3 / t1)
    fp2_set_one(&out->codomain.x);
    fp2_mul(&tmp, &t4, &t_inv[1]);
    fp2_copy(&out->precomputation.y, &tmp);
    fp2_mul(&tmp, &t3, &t_inv[0]);
    fp2_copy(&out->precomputation.z, &tmp);
    fp2_set_zero(&out->codomain.t);

    // compute the final codomain
    hadamard(&out->codomain, &out->codomain);

    return 1;
}

// sub routine of the gluing eval
void
gluing_eval_point(theta_point_t *image,
                  const theta_couple_point_t *P,
                  const theta_couple_point_t *Pt,
                  const theta_gluing_t *phi)
{

    theta_point_t T, Tt;

    // apply the basis change
    base_change(&T, phi, P);
    base_change(&Tt, phi, Pt);

    // apply the to_squared_theta transform
    to_squared_theta(&T, &T);
    to_squared_theta(&Tt, &Tt);

    // compute x, y, z
    fp2_copy(&image->x, &T.x);
    fp2_copy(&image->y, &T.y);
    fp2_mul(&image->y, &image->y, &phi->precomputation.y);
    fp2_copy(&image->z, &T.z);
    fp2_mul(&image->z, &image->z, &phi->precomputation.z);

    // Usually here we need to handle the case when z is zero or
    // t is zero, however for our gluing we always have that the
    // zero index is 3 so x,y,z != 0 and t = 0.
#ifndef NDEBUG
    if (!fp2_is_zero(&T.t))
        debug_print("T.t is not zero");
#endif

    fp2_t lam;
    fp2_copy(&lam, &Tt.x);
    fp2_copy(&image->t, &image->y);

    //  normalize using lam
    fp2_mul(&image->x, &image->x, &lam);
    fp2_mul(&image->y, &image->y, &lam);
    fp2_mul(&image->z, &image->z, &lam);

    // recover t
    fp2_mul(&image->t, &image->t, &Tt.z);
    fp2_mul(&image->t, &image->t, &phi->precomputation.z);

    hadamard(image, image);
}

// same as gluing_eval_point but in the very special case where we already know that the point will
// have a zero coordinate at the place where the zero coordinate of the dual_theta_nullpoint would
// have made the computation difficult
void
gluing_eval_point_special_case(theta_point_t *image,
                               const theta_couple_point_t *P,
                               const theta_gluing_t *phi)
{
    theta_point_t T;

    // Apply the basis change
    base_change(&T, phi, P);

    // Apply the to_squared_theta transform
    to_squared_theta(&T, &T);

    // Compute (x, y, z, t)
    fp2_copy(&image->x, &T.x);
    fp2_copy(&image->y, &T.y);
    fp2_mul(&image->y, &image->y, &phi->precomputation.y);
    fp2_copy(&image->z, &T.z);
    fp2_mul(&image->z, &image->z, &phi->precomputation.z);
    fp2_set_zero(&image->t);

    hadamard(image, image);
}

/**
 * @brief Evaluate a gluing isogeny from an elliptic product on an element of a basis
 *
 * @param image Output: the theta_point of the image of the first couple of points
 * @param xyT: A pair of points (X : Y : Z) on E1E2 to glue using phi
 * @param E1E2 : an elliptic product
 * @param phi : a gluing isogeny E1 x E2 -> A
 *
 **/
void
gluing_eval_basis_point(theta_point_t *image,
                        const theta_couple_jac_point_t *xyT,
                        const theta_couple_curve_t *E12,
                        const theta_gluing_t *phi)
{
    theta_couple_jac_point_t T;
    theta_couple_point_t P, Pt;

    // add the point to push with xyK1_4
    add_couple_jac_points(&T, xyT, &phi->xyK1_4, E12);

    // Convert to (X : Z) coordinates
    couple_jac_to_xz(&Pt, &T);
    couple_jac_to_xz(&P, xyT);

    // then we evaluate the gluing
    gluing_eval_point(image, &P, &Pt, phi);
}

/**
 * @brief Evaluate a gluing isogeny from an elliptic product on a basis
 *
 * @param image1 Output: the theta_point of the image of the first couple of points
 * @param image2 Output : the theta point of the image of the second couple of points
 * @param xyT1: A pair of points (X : Y : Z) on E1E2 to glue using phi
 * @param xyT2: A pair of points (X : Y : Z) on E1E2 to glue using phi
 * @param E1E2 : an elliptic product
 * @param phi : a gluing isogeny E1 x E2 -> A
 *
 **/
void
gluing_eval_basis(theta_point_t *image1,
                  theta_point_t *image2,
                  const theta_couple_jac_point_t *xyT1,
                  const theta_couple_jac_point_t *xyT2,
                  const theta_couple_curve_t *E1E2,
                  const theta_gluing_t *phi)
{
    gluing_eval_basis_point(image1, xyT1, E1E2, phi);
    gluing_eval_basis_point(image2, xyT2, E1E2, phi);
}

/**
 * @brief Compute  a (2,2) isogeny in dimension 2 in the theta_model
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_8 a point in A[8]
 * @param T2_8 a point in A[8]
 * @param hadamard_bool_1 a boolean used for the last two steps of the chain
 * @param hadamard_bool_2 a boolean used for the last two steps of the chain
 *
 * out : A -> B of kernel [4](T1_8,T2_8)
 * hadamard_bool_1 controls if the domain is in standard or dual coordinates
 * hadamard_bool_2 controls if the codomain is in standard or dual coordinates
 *
 */
void
theta_isogeny_compute(theta_isogeny_t *out,
                      const theta_structure_t *A,
                      const theta_point_t *T1_8,
                      const theta_point_t *T2_8,
                      bool hadamard_bool_1,
                      bool hadamard_bool_2)
{
    out->hadamard_bool_1 = hadamard_bool_1;
    out->hadamard_bool_2 = hadamard_bool_2;
    out->domain = *A;
    out->T1_8 = *T1_8;
    out->T2_8 = *T2_8;
    out->codomain.precomputation = false;

    theta_point_t TT1, TT2;

    if (hadamard_bool_1) {
        hadamard(&TT1, T1_8);
        to_squared_theta(&TT1, &TT1);
        hadamard(&TT2, T2_8);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT1, T1_8);
        to_squared_theta(&TT2, T2_8);
    }

    fp2_t t1, t2;
    fp2_mul(&t1, &TT1.x, &TT2.y);
    fp2_mul(&t2, &TT1.y, &TT2.x);
    fp2_mul(&out->codomain.null_point.x, &TT2.x, &t1);
    fp2_mul(&out->codomain.null_point.y, &TT2.y, &t2);
    fp2_mul(&out->codomain.null_point.z, &TT2.z, &t1);
    fp2_mul(&out->codomain.null_point.t, &TT2.t, &t2);
    fp2_t t3;
    fp2_mul(&t3, &TT2.z, &TT2.t);
    fp2_mul(&out->precomputation.x, &t3, &TT1.y);
    fp2_mul(&out->precomputation.y, &t3, &TT1.x);
    fp2_copy(&out->precomputation.z, &out->codomain.null_point.t);
    fp2_copy(&out->precomputation.t, &out->codomain.null_point.z);

    if (hadamard_bool_2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
}

/**
 * @brief Compute a (2,2) isogeny when only the 4 torsion above the kernel is known and not the 8
 * torsion
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_4 a point in A[4]
 * @param T2_4 a point in A[4]
 * @param hadamard_bool_1 a boolean
 * @param hadamard_bool_2 a boolean
 *
 * out : A -> B of kernel [2](T1_4,T2_4)
 * hadamard_bool_1 controls if the domain is in standard or dual coordinates
 * hadamard_bool_2 controls if the codomain is in standard or dual coordinates
 *
 */
void
theta_isogeny_compute_4(theta_isogeny_t *out,
                        const theta_structure_t *A,
                        const theta_point_t *T1_4,
                        const theta_point_t *T2_4,
                        bool hadamard_bool_1,
                        bool hadamard_bool_2)
{
    out->hadamard_bool_1 = hadamard_bool_1;
    out->hadamard_bool_2 = hadamard_bool_2;
    out->domain = *A;
    out->T1_8 = *T1_4;
    out->T2_8 = *T2_4;
    out->codomain.precomputation = false;

    theta_point_t TT1, TT2;
    // we will compute:
    // TT1 = (xAB, _ , xCD, _)
    // TT2 = (AA,BB,CC,DD)

    // fp2_t xA_inv,zA_inv,tB_inv;

    if (hadamard_bool_1) {
        hadamard(&TT1, T1_4);
        to_squared_theta(&TT1, &TT1);

        hadamard(&TT2, &A->null_point);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT1, T1_4);
        to_squared_theta(&TT2, &A->null_point);
    }

    fp2_t sqaabb, sqaacc;
    fp2_mul(&sqaabb, &TT2.x, &TT2.y);
    fp2_mul(&sqaacc, &TT2.x, &TT2.z);
    // sqaabb = sqrt(AA*BB)
    fp2_sqrt(&sqaabb);
    // sqaacc = sqrt(AA*CC)
    fp2_sqrt(&sqaacc);

    // we compute out->codomain.null_point = (xAB * sqaacc * AA, xAB *sqaabb *sqaacc, xCD*sqaabb *
    // AA) out->precomputation = (xAB * BB * CC *DD , sqaabb * CC * DD * xAB , sqaacc * BB* DD * xAB
    // , xCD * sqaabb *sqaacc * BB)

    fp2_mul(&out->codomain.null_point.y, &sqaabb, &sqaacc);
    fp2_mul(&out->precomputation.t, &out->codomain.null_point.y, &TT1.z);
    fp2_mul(&out->codomain.null_point.y,
            &out->codomain.null_point.y,
            &TT1.x); // done for out->codomain.null_point.y

    fp2_mul(&out->codomain.null_point.t, &TT1.z, &sqaabb);
    fp2_mul(&out->codomain.null_point.t,
            &out->codomain.null_point.t,
            &TT2.x); // done for out->codomain.null_point.t

    fp2_mul(&out->codomain.null_point.x, &TT1.x, &TT2.x);
    fp2_mul(&out->codomain.null_point.z,
            &out->codomain.null_point.x,
            &TT2.z); // done for out->codomain.null_point.z
    fp2_mul(&out->codomain.null_point.x,
            &out->codomain.null_point.x,
            &sqaacc); // done for out->codomain.null_point.x

    fp2_mul(&out->precomputation.x, &TT1.x, &TT2.t);
    fp2_mul(&out->precomputation.z, &out->precomputation.x, &TT2.y);
    fp2_mul(&out->precomputation.x, &out->precomputation.x, &TT2.z);
    fp2_mul(
        &out->precomputation.y, &out->precomputation.x, &sqaabb); // done for out->precomputation.y
    fp2_mul(
        &out->precomputation.x, &out->precomputation.x, &TT2.y); // done for out->precomputation.x
    fp2_mul(
        &out->precomputation.z, &out->precomputation.z, &sqaacc); // done for out->precomputation.z
    fp2_mul(
        &out->precomputation.t, &out->precomputation.t, &TT2.y); // done for out->precomputation.t

    if (hadamard_bool_2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
}

/**
 * @brief Compute a (2,2) isogeny when only the kernel is known and not the 8 or 4 torsion above
 *
 * @param out Output: the theta_isogeny
 * @param A a theta null point for the domain
 * @param T1_2 a point in A[2]
 * @param T2_2 a point in A[2]
 * @param hadamard_bool_1 a boolean
 * @param boo2 a boolean
 *
 * out : A -> B of kernel (T1_2,T2_2)
 * hadamard_bool_1 controls if the domain is in standard or dual coordinates
 * hadamard_bool_2 controls if the codomain is in standard or dual coordinates
 *
 */
void
theta_isogeny_compute_2(theta_isogeny_t *out,
                        const theta_structure_t *A,
                        const theta_point_t *T1_2,
                        const theta_point_t *T2_2,
                        bool hadamard_bool_1,
                        bool hadamard_bool_2)
{
    out->hadamard_bool_1 = hadamard_bool_1;
    out->hadamard_bool_2 = hadamard_bool_2;
    out->domain = *A;
    out->T1_8 = *T1_2;
    out->T2_8 = *T2_2;
    out->codomain.precomputation = false;

    theta_point_t TT2;
    // we will compute:
    // TT2 = (AA,BB,CC,DD)

    if (hadamard_bool_1) {
        hadamard(&TT2, &A->null_point);
        to_squared_theta(&TT2, &TT2);
    } else {
        to_squared_theta(&TT2, &A->null_point);
    }

    // we compute out->codomain.null_point = (AA,sqaabb, sqaacc, sqaadd)
    // out->precomputation = (  BB * CC *DD , sqaabb * CC * DD , sqaacc * BB* DD , sqaadd * BB * CC)
    fp2_copy(&out->codomain.null_point.x, &TT2.x);
    fp2_mul(&out->codomain.null_point.y, &TT2.x, &TT2.y);
    fp2_mul(&out->codomain.null_point.z, &TT2.x, &TT2.z);
    fp2_mul(&out->codomain.null_point.t, &TT2.x, &TT2.t);
    fp2_sqrt(&out->codomain.null_point.y);
    fp2_sqrt(&out->codomain.null_point.z);
    fp2_sqrt(&out->codomain.null_point.t);

    fp2_mul(&out->precomputation.x, &TT2.z, &TT2.t);
    fp2_mul(&out->precomputation.y,
            &out->precomputation.x,
            &out->codomain.null_point.y); // done for out->precomputation.y
    fp2_mul(
        &out->precomputation.x, &out->precomputation.x, &TT2.y); // done for out->precomputation.x
    fp2_mul(&out->precomputation.z, &TT2.t, &out->codomain.null_point.z);
    fp2_mul(
        &out->precomputation.z, &out->precomputation.z, &TT2.y); // done for out->precomputation.z
    fp2_mul(&out->precomputation.t, &TT2.z, &out->codomain.null_point.t);
    fp2_mul(
        &out->precomputation.t, &out->precomputation.t, &TT2.y); // done for out->precomputation.t

    if (hadamard_bool_2) {
        hadamard(&out->codomain.null_point, &out->codomain.null_point);
    }
}

void
theta_isogeny_eval(theta_point_t *out, const theta_isogeny_t *phi, const theta_point_t *P)
{
    if (phi->hadamard_bool_1) {
        hadamard(out, P);
        to_squared_theta(out, out);
    } else {
        to_squared_theta(out, P);
    }
    fp2_mul(&out->x, &out->x, &phi->precomputation.x);
    fp2_mul(&out->y, &out->y, &phi->precomputation.y);
    fp2_mul(&out->z, &out->z, &phi->precomputation.z);
    fp2_mul(&out->t, &out->t, &phi->precomputation.t);

    if (phi->hadamard_bool_2) {
        hadamard(out, out);
    }
}

int
splitting_compute(theta_splitting_t *out, const theta_structure_t *A)
{
    // init
    uint32_t ctl;
    uint32_t good = 0;
    fp2_t U_cst, t1, t2;

    // enumerate through all indices
    for (int i = 0; i < 10; i++) {
        fp2_set_zero(&U_cst);
        for (int t = 0; t < 4; t++) {
            // Iterate through the null point
            choose_index_theta_point(&t2, t, &A->null_point);
            choose_index_theta_point(&t1, t ^ EVEN_INDEX[i][1], &A->null_point);

            // Compute t1 * t2
            fp2_mul(&t1, &t1, &t2);
            // If CHI_EVAL(i,t) is +1 we want ctl to be 0 and
            // If CHI_EVAL(i,t) is -1 we want ctl to be 0xFF..FF
            ctl = (uint32_t)(CHI_EVAL[EVEN_INDEX[i][0]][t] >> 1);
            assert(ctl == 0 || ctl == 0xffffffff);

            fp2_neg(&t2, &t1);
            fp2_select(&t1, &t1, &t2, ctl);

            // Then we compute U_cst ± (t1 * t2)
            fp2_add(&U_cst, &U_cst, &t1);
        }

        // If U_cst update the good value with 1 + i
        ctl = fp2_is_zero(&U_cst);
        good |= (1 + i) & ctl;
    }

    // now we apply the isomorphism if it was computed
    if (good) {
        // Select the matrix from precomputations
        out->M = SPLITTING_TRANSFORMS[good - 1];
        apply_isomorphism(&out->B.null_point, &out->M, &A->null_point);
    }
    return good;
}

void
theta_product_structure_to_elliptic_product(theta_couple_curve_t *E12, theta_structure_t *A)
{
    fp2_t xx, yy, temp1, temp2;

    ec_curve_init(&(E12->E1));
    ec_curve_init(&(E12->E2));

    // xx = x², yy = y²
    fp2_sqr(&xx, &A->null_point.x);
    fp2_sqr(&yy, &A->null_point.y);

    // A1 = ( (xx² + yy²)² + (xx² - yy²)² )  / (xx² - yy²)
    fp2_add(&temp1, &xx, &yy);
    fp2_sub(&temp2, &xx, &yy);
    fp2_mul(&E12->E1.C, &temp1, &temp2);
    fp2_sqr(&temp1, &temp1);
    fp2_sqr(&temp2, &temp2);
    fp2_add(&temp1, &temp1, &temp2);
    fp2_set_one(&E12->E1.A);
    fp2_mul(&E12->E1.A, &E12->E1.A, &temp1);
    fp2_neg(&E12->E1.A, &E12->E1.A);

    // same with y,t
    fp2_sqr(&xx, &A->null_point.y);
    fp2_sqr(&yy, &A->null_point.t);

    // A2 = ( (xx² + yy²)² + (xx² - yy²)² )  / (xx² - yy²)
    fp2_add(&temp1, &xx, &yy);
    fp2_sub(&temp2, &xx, &yy);
    fp2_mul(&E12->E2.C, &temp1, &temp2);
    fp2_sqr(&temp1, &temp1);
    fp2_sqr(&temp2, &temp2);
    fp2_add(&temp1, &temp1, &temp2);
    fp2_set_one(&E12->E2.A);
    fp2_mul(&E12->E2.A, &E12->E2.A, &temp1);
    fp2_neg(&E12->E2.A, &E12->E2.A);
}

void
theta_point_to_montgomery_point(theta_couple_point_t *P12,
                                const theta_point_t *P,
                                const theta_structure_t *A)
{
    fp2_t temp;

    // P1.X = A.null_point.y * P.x + A.null_point.x * P.y
    // P1.Z = - A.null_point.y * P.x + A.null_point.x * P.y
    fp2_mul(&P12->P1.x, &A->null_point.y, &P->x);
    fp2_mul(&temp, &A->null_point.x, &P->y);
    fp2_sub(&P12->P1.z, &temp, &P12->P1.x);
    fp2_add(&P12->P1.x, &P12->P1.x, &temp);

    // P2.X = A.null_point.t * P.y + A.null_point.y * P.t
    // P2.Z = A.null_point.t * P.y - A.null_point.y * P.t
    fp2_mul(&P12->P2.x, &A->null_point.t, &P->y);
    fp2_mul(&temp, &A->null_point.y, &P->t);
    fp2_sub(&P12->P2.z, &temp, &P12->P2.x);
    fp2_add(&P12->P2.x, &P12->P2.x, &temp);
}

static int
_theta_chain_compute_impl(theta_chain_t *out,
                          unsigned n,
                          theta_couple_curve_t *E12,
                          const theta_couple_point_t *T1,
                          const theta_couple_point_t *T2,
                          const theta_couple_point_t *T1m2,
                          bool eight_above,
                          theta_couple_curve_t *E34,
                          theta_couple_point_t *P12,
                          size_t numP)
{
    theta_structure_t theta;
    theta_point_t R1, R2;

    const unsigned adjusting = eight_above ? 0 : 2;

    const size_t strat_len = sizeof(*strategies) / sizeof(**strategies);
    assert(strat_len - n + adjusting <= sizeof(strategies) / sizeof(*strategies));
    const uint16_t *strategy = strategies[strat_len - n + adjusting];

    // init of the isogeny chain
    if (out) {
        theta_chain_init(out, n - 1);
        out->domain = *E12;
        out->length = n;
        out->T1 = *T1;
        out->T2 = *T2;
    }

    // lift the basis
    theta_couple_jac_point_t xyT1, xyT2;
    ec_basis_t bas1 = { .P = T1->P1, .Q = T2->P1, .PmQ = T1m2->P1 };
    ec_basis_t bas2 = { .P = T1->P2, .Q = T2->P2, .PmQ = T1m2->P2 };
    lift_basis(&xyT1.P1, &xyT2.P1, &bas1, &E12->E1);
    lift_basis(&xyT1.P2, &xyT2.P2, &bas2, &E12->E2);

#ifndef NDEBUG
    if (!test_point_order_twof(&bas2.P, &E12->E2, n + 2 * eight_above))
        debug_print("bas2.P does not have correct order");

    if (!test_jac_order_twof(&xyT2.P2, &E12->E2, n + 2 * eight_above))
        debug_print("xyT2.P2 does not have correct order");
#endif

    theta_couple_jac_point_t jacQ1[n], jacQ2[n];
    unsigned num_dbls[n];

    // the first point
    jacQ1[0] = xyT1;
    jacQ2[0] = xyT2;
    num_dbls[0] = 0;
    size_t count = 1;

    // and then the rest
    for (unsigned tot_dbls = 0; tot_dbls != n - 1 - adjusting; ++count) {
        assert(tot_dbls < n - 1 - adjusting);
        assert(count < n);
        assert(*strategy);
        double_couple_jac_point_iter(&jacQ1[count], *strategy, &jacQ1[count - 1], E12);
        double_couple_jac_point_iter(&jacQ2[count], *strategy, &jacQ2[count - 1], E12);
        tot_dbls += num_dbls[count] = *strategy++;
        assert(tot_dbls <= n - 1 - adjusting);
    }

    // compute the gluing isogeny
    theta_gluing_t first_step;
    if (!gluing_compute(&first_step, E12, &jacQ1[count - 1], &jacQ2[count - 1]))
        return 0;

    if (out)
        out->first_step = first_step;

    // evaluate
    theta_point_t pts[numP ? numP : 1];
    for (size_t j = 0; j < numP; ++j) {
        assert(ec_is_zero(&P12[j].P1) ||
               ec_is_zero(&P12[j].P2)); // only theta_chain_eval_special_case() implemented here
        gluing_eval_point_special_case(&pts[j], &P12[j], &first_step);
    }

    // set-up the theta_structure for the first codomain
    theta.null_point = first_step.codomain;
    theta.precomputation = 0;
    theta_precomputation(&theta);

    assert(count);
    count--;

    theta_point_t thetaQ1[n], thetaQ2[n];

    // push the kernel through the gluing isogeny
    // need to setup the input before
    for (int i = 0; i < count; i++)
        gluing_eval_basis(&thetaQ1[i], &thetaQ2[i], &jacQ1[i], &jacQ2[i], E12, &first_step);

    theta_isogeny_t step;

    // and now we do the remaining steps
    for (unsigned i = 0; i < n - 1 - adjusting; i++) {
        {
            unsigned tot_dbls = 0;
            for (int j = 0; j < count; j++)
                tot_dbls += num_dbls[j];
            for (; tot_dbls != n - i - 2 - adjusting; ++count) {
                assert(tot_dbls < n - 1 - adjusting);
                assert(count < n);
                assert(*strategy);
                double_iter(&thetaQ1[count], &theta, &thetaQ1[count - 1], *strategy);
                double_iter(&thetaQ2[count], &theta, &thetaQ2[count - 1], *strategy);
                tot_dbls += num_dbls[count] = *strategy++;
                assert(tot_dbls <= n - i - 2 - adjusting);
            }
        }

        // computing the next step
        if (eight_above && i == n - 3)
            theta_isogeny_compute(&step, &theta, &thetaQ1[count - 1], &thetaQ2[count - 1], 0, 0);
        else if (eight_above && i == n - 2)
            theta_isogeny_compute(&step, &theta, &thetaQ1[count - 1], &thetaQ2[count - 1], 1, 0);
        else
            theta_isogeny_compute(&step, &theta, &thetaQ1[count - 1], &thetaQ2[count - 1], 0, 1);

        if (out)
            out->steps[i] = step;
        for (size_t j = 0; j < numP; ++j)
            theta_isogeny_eval(&pts[j], &step, &pts[j]);

        // updating the codomain
        theta = step.codomain;
        assert(count);
        count--;

        // pushing the kernel
        if (i < n - 2 - adjusting) {
            for (int j = 0; j < count; j++) {
                theta_isogeny_eval(&thetaQ1[j], &step, &thetaQ1[j]);
                theta_isogeny_eval(&thetaQ2[j], &step, &thetaQ2[j]);
            }
        }
    }
    assert(!*strategy);
    assert(!count);

    if (!eight_above) {
        assert(n - 1 - adjusting - 1 == n - 4); // ensure current value of step is the correct one
        // the last two steps are done here
        theta_isogeny_eval(&thetaQ1[0], &step, &thetaQ1[0]);
        theta_isogeny_eval(&thetaQ2[0], &step, &thetaQ2[0]);
        theta_isogeny_compute_4(&step, &theta, &thetaQ1[0], &thetaQ2[0], 0, 0);
        if (out)
            out->steps[n - 3] = step;
        for (size_t j = 0; j < numP; ++j)
            theta_isogeny_eval(&pts[j], &step, &pts[j]);
        theta = step.codomain;
        theta_isogeny_eval(&thetaQ1[0], &step, &thetaQ1[0]);
        theta_isogeny_eval(&thetaQ2[0], &step, &thetaQ2[0]);
        theta_isogeny_compute_2(&step, &theta, &thetaQ1[0], &thetaQ2[0], 1, 0);
        if (out)
            out->steps[n - 2] = step;
        for (size_t j = 0; j < numP; ++j)
            theta_isogeny_eval(&pts[j], &step, &pts[j]);
        theta = step.codomain;
    }

    // final splitting step
    theta_splitting_t last_step;
    int is_split = splitting_compute(&last_step, &step.codomain);

    if (!is_split) {
        debug_print("kernel did not generate an isogeny between elliptic products");
        return 0;
    }

    if (out) {
        out->last_step = last_step;
        E34 = &out->codomain;
    }
    theta_product_structure_to_elliptic_product(E34, &last_step.B);

    // evaluate
    for (size_t j = 0; j < numP; ++j) {
        apply_isomorphism(&pts[j], &last_step.M, &pts[j]);
        theta_point_to_montgomery_point(&P12[j], &pts[j], &last_step.B);
    }

    return 1;
}

int
theta_chain_compute(theta_chain_t *out,
                    unsigned n,
                    theta_couple_curve_t *E12,
                    const theta_couple_point_t *T1,
                    const theta_couple_point_t *T2,
                    const theta_couple_point_t *T1m2,
                    bool eight_above)
{
    return _theta_chain_compute_impl(out, n, E12, T1, T2, T1m2, eight_above, NULL, NULL, 0);
}

int
theta_chain_compute_and_eval(unsigned n,
                             /*const*/ theta_couple_curve_t *E12,
                             const theta_couple_point_t *T1,
                             const theta_couple_point_t *T2,
                             const theta_couple_point_t *T1m2,
                             bool eight_above,
                             theta_couple_curve_t *E34,
                             theta_couple_point_t *P12,
                             size_t numP)
{
    return _theta_chain_compute_impl(NULL, n, E12, T1, T2, T1m2, eight_above, E34, P12, numP);
}

void
theta_chain_eval_special_case(theta_couple_point_t *out,
                              const theta_chain_t *phi,
                              const theta_couple_point_t *P12)
{
    theta_point_t temp;

    assert(ec_is_zero(&P12->P1) || ec_is_zero(&P12->P2));

    // first, we apply the gluing
    gluing_eval_point_special_case(&temp, P12, &phi->first_step);

    // then, we apply the successive isogenies
    for (int i = 0; i < phi->length - 1; i++) {
        theta_isogeny_eval(&temp, &phi->steps[i], &temp);
    }

    // we send the result to the theta product structure of the codomain
    apply_isomorphism(&temp, &phi->last_step.M, &temp);

    // finally we send the result to the elliptic product
    theta_point_to_montgomery_point(out, &temp, &phi->last_step.B);
}
