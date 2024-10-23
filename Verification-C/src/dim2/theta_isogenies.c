#include <theta_isogenies.h>
#include <stdio.h>
#include <inttypes.h>
#include <assert.h>
#include <tools.h>
#include <dim2_change_coords.h>
#include <ec.h>
#include <fp2_matrix.h>


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

static inline void
table_to_theta_point(theta_point_t *res, const fp2_t *in)
{
    fp2_copy(&res->x,&in[0]);
    fp2_copy(&res->y,&in[1]);
    fp2_copy(&res->z,&in[2]);
    fp2_copy(&res->t,&in[3]);
}

static inline void
theta_point_to_table(fp2_t *res, const theta_point_t *in)
{
    fp2_copy(&in[0],&res->x);
    fp2_copy(&in[1],&res->y);
    fp2_copy(&in[2],&res->z);
    fp2_copy(&in[3],&res->t);
}

static inline void
invert_theta_point_proj(theta_point_t *res, const theta_point_t *in)
{
    fp2_t t[4];
    theta_point_to_table(t,in);
    fp2_proj_batched_inv(t, 4);
    table_to_theta_point(res,t);
}


static void gluing_compute(theta_gluing_t *out,
               theta_couple_curve_t *E12,
               const theta_point_t *prod_null_point,
               const theta_couple_jac_point_t *xyK1_8,
               const theta_couple_jac_point_t *xyK2_8,
               const fp2_t **N)
{
    theta_point_t xAxByCyD, zAtBzYtD;
    int zero_idx;
    fp2_t t1, t2, t3, t4, t_inv[4], ABCD[4], precomp[4];

	// Given points in E[8] x E[8] we need the four torsion below
    double_couple_jac_point(&out->xyK1_4, xyK1_8, E12);
    double_couple_jac_point(&out->xyK2_4, xyK2_8, E12);

    // Convert from (X:Y:Z) coordinates to (X:Z)
    theta_couple_point_t K1_8, K2_8;

    couple_jac_to_xz(&K1_8, xyK1_8);
    couple_jac_to_xz(&K2_8, xyK2_8);
    couple_jac_to_xz(&out->K1_4, &out->xyK1_4);
    couple_jac_to_xz(&out->K2_4, &out->xyK2_4);

    // compute the change of coordinates matrix
    montgomery_to_theta_matrix_dim2(out->M, &out->np_domain, prod_null_point, N);
    
    // apply the change of coordinates to the kernel
    montgomery_to_theta(&out->T1_8, &K1_8, out->M);
    montgomery_to_theta(&out->T2_8, &K2_8, out->M);

    /* Compute the codomain */
    to_squared_theta(&xAxByCyD, &out->T1_8);
    to_squared_theta(&zAtBzYtD, &out->T2_8);

    // Find zero index
    zero_idx = (fp2_is_zero(&xAxByCyD.x)&0)|(fp2_is_zero(&xAxByCyD.y)&1)|(fp2_is_zero(&xAxByCyD.z)&2)|(fp2_is_zero(&xAxByCyD.t)&3);
    out->zero_idx=zero_idx;

    choose_index_theta_point(&t1, 1^zero_idx, &zAtBzYtD);
    choose_index_theta_point(&t2, 2^zero_idx, &xAxByCyD);
    choose_index_theta_point(&t3, 3^zero_idx, &zAtBzYtD);
    choose_index_theta_point(&t4, 3^zero_idx, &xAxByCyD);

    fp2_copy(&t_inv[0], &t1);
    fp2_copy(&t_inv[1], &t2);
    fp2_copy(&t_inv[2], &t3);
    fp2_copy(&t_inv[3], &t4);
    fp2_proj_batched_inv(t_inv, 4);

    // Compute dual codomain theta constants A, B, C, D
    fp2_set_zero(&ABCD[0 ^ zero_idx]);
    fp2_mul(&ABCD[1 ^ zero_idx], &t1, &t_inv[2]);
    fp2_mul(&ABCD[2 ^ zero_idx], &t2, &t_inv[3]);
    fp2_set_one(&ABCD[3 ^ zero_idx]);

    // Compute inverse of dual theta constants 1/A, 1/B, 1/C, 1/D (when defined)
    fp2_set_zero(&precomp[0 ^ zero_idx]);
    fp2_mul(&precomp[1 ^ zero_idx], &t_inv[0], &t3);
    fp2_mul(&precomp[2 ^ zero_idx], &t_inv[1], &t4);
    fp2_set_one(&precomp[3 ^ zero_idx]);

    table_to_theta_point(&out->codomain,ABCD);
    table_to_theta_point(&out->precomputation,precomp);

    // Compute the theta null point a, b, c, d
    hadamard(&out->codomain, &out->codomain);
    copy_theta_point(&out->prod_null_point,prod_null_point);
}

// sub routine of the gluing eval
// P is the point we want to evaluate Pt is its translate by T1_4.
void
gluing_eval_sub(theta_point_t *image,
                  const theta_couple_point_t *P,
                  const theta_couple_point_t *Pt,
                  const theta_gluing_t *phi)
{
    theta_point_t T, Tt;
    fp2_t t1, t2, t3, lam_num, lam_den, xyzt[4];
    uint32_t is_zero_z;

    // apply the basis change
    montgomery_to_theta(&T, P, phi->M);
    montgomery_to_theta(&Tt, Pt, phi->M);

    // apply the to_squared_theta transform
    to_squared_theta(&T, &T);
    to_squared_theta(&Tt, &Tt);

    // Directly compute y,z,t
    choose_index_theta_point(&image->y, 1 ^ phi->zero_idx, &T);
    choose_index_theta_point(&t1, 1 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&image->y, &image->y, &t1);
    choose_index_theta_point(&image->z, 2 ^ phi->zero_idx, &T);
    choose_index_theta_point(&t1, 2 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&image->z, &image->z, &t1);
    choose_index_theta_point(&image->t, 3 ^ phi->zero_idx, &T);

    // Normalization constants
    choose_index_theta_point(&t1, 3 ^ phi->zero_idx, &Tt);
    choose_index_theta_point(&t2, 2 ^ phi->zero_idx, &Tt);
    choose_index_theta_point(&t3, 2 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&t2,&t2,&t3);

    is_zero_z=fp2_is_zero(&z);
    fp2_select(&lam_num,&t,&z,is_zero_z);
    fp2_select(&lam_den,&t2,&t1,is_zero_z);

    choose_index_theta_point(&image->x,1 ^ phi->zero_idx, &Tt);
    choose_index_theta_point(&t3,1 ^ phi->zero_idx, &phi->precomputation);
    fp2_mul(&image->x,&image->x,&t3);

    fp2_mul(&xyzt[0 ^ phi->zero_idx], &image->x, &lam_num);
    fp2_mul(&xyzt[1 ^ phi->zero_idx], &image->y, &lam_den);
    fp2_mul(&xyzt[2 ^ phi->zero_idx], &image->z, &lam_den);
    fp2_mul(&xyzt[3 ^ phi->zero_idx], &image->t, &lam_den);

    // Finally compute image
    table_to_theta_point(image,xyzt);
    hadamard(image,image);
}

/**
 * @brief Evaluate a gluing isogeny from an elliptic product on a couple of points
 *
 * @param image Output: the theta_point of the image of the couple of points
 * @param xyT: A pair of points (X : Y : Z) on E1E2 to glue using phi
 * @param E1E2 : an elliptic product
 * @param phi : a gluing isogeny E1 x E2 -> A
 *
 **/
void
gluing_eval(theta_point_t *image,
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
 * @biref Compute the dual of a gluing isogeny
 * 
 * @param phi_dual Output: the dual isogeny gluing isogeny (splitting isogeny) \tilde{phi}: A->E1*E2.
 * @param phi: the gluing isogeny phi: E1*E2->A.
 * @param splitting_matrix: the precomputed change of coordinates matrix from non product to product
 * theta coordinates on E1*E2.
 **/
void
dual_gluing_compute(theta_splitting_t *phi_dual, const theta_gluing_t *phi, const fp2_t **splitting_matrix)
{
    //copy_theta_point(&phi_dual->domain, &phi->codomain);
    hadamard(&phi_dual->domain,&phi->codomain);
    ec_copy(&phi_dual->codomain.E1,&phi->codomain.E1);
    ec_copy(&phi_dual->codomain.E2,&phi->codomain.E2);

    invert_theta_point_proj(&phi_dual->precomputation,&phi->np_domain);
    copy_theta_point(&phi_dual->prod_null_point,&phi->prod_null_point);

    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            fp2_copy(&M[i][j],&splitting_matrix[i][j]);
        }
    }
}

void
theta_point_to_montgomery_point(theta_couple_point_t *P12,
                                const theta_point_t *P,
                                const theta_point_t *prod_null_point)
{
    fp2_t x1, z1, x2, z2, temp, Pi, Pj, nulli, nullj, one, minus_one;
    uint32_t is_zero, is_m11;
    int iP, jP, inull, jnull;

    // x1 = prod_null_point.z * P.x + prod_null_point.x * P.z
    // z1 = prod_null_point.z * P.x - prod_null_point.x * P.z
    fp2_mul(&x1,&prod_null_point->z,&P->x);
    fp2_mul(&temp,&prod_null_point->x,&P->z);
    fp2_sub(&z1,&x1,&temp);
    fp2_add(&x1,&x1,&temp);

    is_zero=fp2_is_zero(&x1)&fp2_is_zero(&z1);
    
    fp2_add(&temp,&x1,&z1);
    is_m11=(~is_zero)&fp2_is_zero(&temp);

    // Select index in P and null_point to use
    // If (x1:z1)!= (0:0), (-1:1), then
    // x2 = prod_null_point.y * P.x + prod_null_point.x * P.y
    // z2 = prod_null_point.y * P.x - prod_null_point.x * P.y
    // If (x1:z1) == (0:0), then (x2:z2)=(-1:1) and we recompute
    // x1 = prod_null_point.t * P.y + prod_null_point.y * P.t
    // z1 = prod_null_point.t * P.y - prod_null_point.y * P.t
    // If (x1:z1) == (-1:1), then
    // x2 = prod_null_point.t * P.z + prod_null_point.z * P.t
    // z2 = prod_null_point.t * P.z - prod_null_point.z * P.t
    iP=(0&(~is_zero)&(~is_m11))|(1&is_zero)|(2&is_m11);
    jP=(1&(~is_zero)&(~is_m11))|(3&is_zero)|(3&is_m11);
    inull=(1&(~is_zero)&(~is_m11))|(3&is_zero)|(3&is_m11);
    jnull=(0&(~is_zero)&(~is_m11))|(1&is_zero)|(2&is_m11);

    choose_index_theta_point(&Pi,iP,&P);
    choose_index_theta_point(&Pj,jP,&P);
    choose_index_theta_point(&nulli,inull,&prod_null_point);
    choose_index_theta_point(&nullj,jnull,&prod_null_point);
    
    fp2_mul(&x2,&nulli,&Pi);
    fp2_mul(&temp,&nullj,&Pj);
    fp2_sub(&z2,&x2,&temp);
    fp2_add(&x2,&x2,&temp);

    fp2_set_one(&one);
    fp2_neg(&minus_one,&one);

    // If (x1:z1) == (0:0), then (x1:z1) has been computed above as (x2:z2):
    // x1 = prod_null_point.t * P.y + prod_null_point.y * P.t
    // z1 = prod_null_point.t * P.y - prod_null_point.y * P.t
    // and (x2:z2)=(-1:1)
    fp2_select(&P12->P1.x,&x1,&x2,is_zero);
    fp2_select(&P12->P1.z,&z1,&z2,is_zero);
    fp2_select(&P12->P2.x,&x2,&minus_one,is_zero);
    fp2_select(&P12->P2.z,&z2,&one,is_zero);
}

/**
 * 
 **/
void 
splitting_eval(theta_couple_point_t *image, const theta_point_t *P,const theta_splitting_t *phi)
{
    theta_point_t image_tmp;
    fp2_t prod_coords[4], t[4];

    to_squared_theta(&image_tmp, &P);

    fp2_mul(&image_tmp.x,&image_tmp.x,&phi->precomputation.x);
    fp2_mul(&image_tmp.y,&image_tmp.y,&phi->precomputation.y);
    fp2_mul(&image_tmp.z,&image_tmp.z,&phi->precomputation.z);
    fp2_mul(&image_tmp.t,&image_tmp.t,&phi->precomputation.t);

    theta_point_to_table(t,image_tmp);

    mat_vec_prod(prod_coords, phi->M, t, 4, 4);

    table_to_theta_point(image_tmp,prod_coords);

    theta_point_to_montgomery_point(image, image_tmp,phi->prod_null_point);
}