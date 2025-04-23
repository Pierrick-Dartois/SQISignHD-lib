#include <sqisigndim2_heuristic.h>
#include <curve_extras.h>
#include <tools.h>
#include <fips202.h>
#include <stdio.h>
#include <string.h>

#define RESPONSE_LENGTH TORSION_PLUS_EVEN_POWER + 16

const clock_t time_isogenies_odd = 0;
const clock_t time_sample_response = 0;
const clock_t time_change_of_basis_matrix = 0;

void
secret_sig_init(signature_t *sig)
{
    sig->hint_aux = (int *)malloc(2 * sizeof(int));
    ibz_init(&sig->x);
    ibz_init(&sig->b0);
    ibz_init(&sig->d0);
    ibz_init(&sig->b1);
    ibz_init(&sig->d1);
    ibz_init(&sig->e0_adjust);
    ibz_init(&sig->c0_adjust);
}

void
secret_sig_finalize(signature_t *sig)
{
    free(sig->hint_aux);
    ibz_finalize(&sig->x);
    ibz_finalize(&sig->b0);
    ibz_finalize(&sig->d0);
    ibz_finalize(&sig->b1);
    ibz_finalize(&sig->d1);
    ibz_finalize(&sig->e0_adjust);
    ibz_finalize(&sig->c0_adjust);
}

static void
ibz_vec_2_print2(char *name, const ibz_vec_2_t *vec)
{
    printf("%s", name);
    for (int i = 0; i < 2; i++) {
        ibz_printf("%Zd ", &((*vec)[i]));
    }
    ibz_printf("\n");
}

static void
ibz_vec_4_print2(char *name, const ibz_vec_4_t *vec)
{
    printf("%s", name);
    for (int i = 0; i < 4; i++) {
        ibz_printf("%Zd ", &((*vec)[i]));
    }
    ibz_printf("\n");
}

// compute the commitment with fixed degree isogeny
// and apply it to the basis of E0
void
commit(ec_curve_t *E_com, quat_left_ideal_t *lideal_com)
{

    int found = 1;
    ibz_t n, adj;
    ibz_init(&n);
    ibz_init(&adj);
    // generate a random ideal of random norm for the secret ideal
    // TODO make a clean constant for this
    generate_random_prime(&n, 1, ibz_bitsize(&QUATALG_PINFTY.p) / 2);

    theta_chain_t F;
    found = fixed_degree_isogeny(&F, lideal_com, &n, &adj, 1);

    // it's always the second curve
    copy_curve(E_com, &F.codomain.E2);

    assert(found);

    ibz_finalize(&n);
    ibz_finalize(&adj);
}

void
quat_lideal_conjugate_lattice(quat_lattice_t *lat, const quat_left_ideal_t *lideal)
{
    ibz_mat_4x4_copy(&(lat->basis), &(lideal->lattice.basis));
    ibz_copy(&(lat->denom), &(lideal->lattice.denom));

    for (int row = 1; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            ibz_neg(&(lat->basis[row][col]), &(lat->basis[row][col]));
        }
    }

    return;
}

void
norm_from_2_times_gram(ibz_t *norm, ibz_mat_4x4_t *gram, ibz_vec_4_t *vec)
{
    quat_qf_eval(norm, gram, vec);
    assert(ibz_is_even(norm));
    ibz_div_2exp(norm, norm, 1);
}

int
sample_response(quat_alg_elem_t *x,
                const quat_lattice_t *lattice,
                ibz_t const *lattice_content,
                int verbose)
{
    ibz_mat_4x4_t lll;
    ibz_t denom_gram, norm, norm_bound;

    ibz_mat_4x4_init(&lll);
    ibz_init(&denom_gram);
    ibz_init(&norm);
    ibz_init(&norm_bound);

    ibz_pow(&norm_bound, &ibz_const_two, SQIsign2D_response_heuristic_bound);

    int err = quat_lattice_lll(&lll, lattice, &(QUATALG_PINFTY.p));
    assert(!err);
    // The shortest vector found by lll is a candidate

    int found;

    ibz_mat_4x4_t prod, gram;
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_init(&gram);

    ibz_mat_4x4_transpose(&prod, &lll);
    ibz_mat_4x4_mul(&prod, &prod, &(QUATALG_PINFTY.gram));
    ibz_mat_4x4_mul(&gram, &prod, &lll);

    ibz_copy(&denom_gram, &(lattice->denom));
    ibz_mul(&denom_gram, &denom_gram, &(lattice->denom));
    ibz_mul(&denom_gram, &denom_gram, lattice_content);

    assert(ibz_is_even(&denom_gram));
    ibz_div_2exp(&denom_gram, &denom_gram, 1);

    int divides = ibz_mat_4x4_scalar_div(&gram, &denom_gram, &gram);
    assert(divides);

    ibz_copy(&(x->denom), &(lattice->denom));

    ibz_vec_4_t vec;
    ibz_vec_4_init(&vec);

    found = 0;
    int cnt = 0;

    // TODO make these clean constants
    int m = 3;
    while (!found && cnt < 2 * (2 * m + 1) * (2 * m + 1) * (2 * m + 1) * (2 * m + 1)) {
        cnt++;
        for (int i = 0; i < 4; i++) {
            ibz_rand_interval_minm_m(&vec[i], m);
        }
        norm_from_2_times_gram(&norm, &gram, &vec);
        // now we test if the norm is good
        // there are two constraints : the norm must be smaller than some bound
        // and the element must be primitive in O0 (this ensures that there is no backtracking)
        found = ibz_cmp(&norm, &norm_bound) < 0 &&
                (ibz_cmp(&vec[0], &ibz_const_zero) != 0 || ibz_cmp(&vec[1], &ibz_const_zero) != 0 ||
                 ibz_cmp(&vec[2], &ibz_const_zero) != 0 || ibz_cmp(&vec[3], &ibz_const_zero) != 0);
        if (found) {
            ibz_mat_4x4_eval(&(x->coord), &lll, &vec);
            assert(quat_lattice_contains(NULL, lattice, x, &QUATALG_PINFTY));
            if (!quat_alg_is_primitive(x, &MAXORD_O0, &QUATALG_PINFTY)) {
                found = 0;
            }
        }
    }
    assert(quat_lattice_contains(NULL, lattice, x, &QUATALG_PINFTY));

    ibz_finalize(&denom_gram);
    ibz_finalize(&norm);
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&lll);
    ibz_vec_4_finalize(&vec);
    ibz_finalize(&norm_bound);
    return found;
}

// compute the challenge as the hash of the message and the commitment curve and public key
void
hash_to_challenge(ibz_vec_2_t *scalars,
                  const ec_curve_t *com_curve,
                  const unsigned char *message,
                  const public_key_t *pk,
                  size_t length)
{
    unsigned char *buf = malloc(FP2_ENCODED_BYTES + FP2_ENCODED_BYTES + length);
    {
        fp2_t j1, j2;
        ec_j_inv(&j1, com_curve);
        ec_j_inv(&j2, &pk->curve);
        fp2_encode(buf, &j1);
        fp2_encode(buf + FP2_ENCODED_BYTES, &j2); // TODO use defined constant
        memcpy(buf + FP2_ENCODED_BYTES + FP2_ENCODED_BYTES,
               message,
               length); // TODO use defined constant
    }

    // TODO(security) omit some vectors, notably (a,1) with gcd(a,6)!=1 but also things like (2,3)?
    {
        digit_t digits[NWORDS_FIELD];

        // FIXME should use SHAKE128 for smaller parameter sets?
        //  TODO we want to use a bit differently (first we hash first half and then derive the
        //  second half)
        SHAKE256(
            (void *)digits, sizeof(digits), buf, FP2_ENCODED_BYTES + FP2_ENCODED_BYTES + length);
        for (int i = 0; i < SQIsign2D_heuristic_challenge_hash_iteration; i++) {
            SHAKE256((void *)digits, sizeof(digits), (void *)digits, sizeof(digits));
        }

        ibz_set(&(*scalars)[1], 1); // FIXME
        ibz_copy_digit_array(&(*scalars)[1], digits);
    }

    ibz_set(&((*scalars)[0]), 1);
    // ibz_rand_interval(&((*scalars)[1]), &((*scalars)[0]), &TORSION_PLUS_2POWER);

    free(buf);
}

int
protocols_sign(signature_t *sig,
               const public_key_t *pk,
               secret_key_t *sk,
               const unsigned char *m,
               size_t l,
               int verbose)
{
    clock_t t = tic();

    ibz_t lattice_content, pow_chall;
    ec_curve_t E_aux, E_com;
    ec_basis_t Baux0; // basis of 2^TORSION_PLUS_EVEN_POWER
    ec_basis_t B_resp_two;
    ibz_vec_4_t coeffs;
    ibz_vec_2_t vec, vec_chall, vec_resp_two;
    quat_left_ideal_t lideal_tmp;
    quat_left_ideal_t lideal_commit, lideal_chall_two;
    quat_left_ideal_t lideal_chall_secret, lideal_com_resp, lideal_aux, lideal_aux_com,
        lideal_resp_two;
    quat_lattice_t lattice_hom_chall_to_com, lat_commit;
    quat_alg_elem_t resp_quat;
    quat_alg_elem_t elem_tmp;
    ibz_mat_2x2_t mat_Baux0_to_Baux_can, mat_Bchall_can_to_Bchall, mat, sig_mat_pk_can_to_B_pk;
    // ibz_mat_2x2_t mat_sigma_phichall_BA_to_Bcomcan, mat_sigma_phichall_BA0_to_Bcom0;
    ibz_t degree_com_isogeny, tmp, remain;
    ibz_t degree_full_resp, degree_odd_resp;
    ibq_t temp_norm;
    int exp_diadic_val_full_resp;
    int pow_dim2_deg_resp;
    int backtracking;

    int found = 1;

    ec_curve_init(&E_aux);
    ec_curve_init(&E_com);
    ibz_init(&tmp);
    ibz_init(&lattice_content);
    ibz_init(&remain);

    ibz_init(&degree_full_resp);
    ibz_init(&degree_odd_resp);
    ibq_init(&temp_norm);
    ibz_init(&pow_chall);

    ibz_mat_2x2_init(&mat_Bchall_can_to_Bchall);
    ibz_mat_2x2_init(&mat_Baux0_to_Baux_can);
    ibz_mat_2x2_init(&sig_mat_pk_can_to_B_pk);
    ibz_mat_2x2_init(&mat);
    ibz_vec_4_init(&coeffs);

    quat_alg_elem_init(&resp_quat);
    quat_alg_elem_init(&elem_tmp);
    quat_lattice_init(&lattice_hom_chall_to_com);
    quat_lattice_init(&lat_commit);
    quat_left_ideal_init(&lideal_tmp);
    quat_left_ideal_init(&lideal_commit);
    quat_left_ideal_init(&lideal_chall_two);
    quat_left_ideal_init(&lideal_chall_secret);
    quat_left_ideal_init(&lideal_resp_two);
    quat_left_ideal_init(&lideal_com_resp);
    quat_left_ideal_init(&lideal_aux_com);
    quat_left_ideal_init(&lideal_aux);

    ibz_vec_2_init(&vec);
    ibz_vec_2_init(&vec_chall);
    ibz_vec_2_init(&vec_resp_two);

    // computing the commitment
    commit(&E_com, &lideal_commit);

    // TODO make a clean constant for this
    int len_chall = SQIsign2D_heuristic_challenge_length;
    ibz_pow(&pow_chall, &ibz_const_two, len_chall);

    // computing the challenge
    // vec_chall is a pair of coefficients encoding the kernel of the challenge isogeny
    // as vec_chall[0]*B[0] + vec_chall[1]*B[1] where B is the canonical basis of the 2^len_chall
    // torsion of EA
    hash_to_challenge(&vec_chall, &E_com, m, pk, l);

    // now we compute the ideal associated to the challenge
    // for that, we need to find vec such that
    // the kernel of the challenge isogeny is generated by vec[0]*B0[0] + vec[1]*B0[1] where B0 is
    // the image through the secret key isogeny of the canonical basis E0
    ibz_mat_2x2_eval(&vec, &(sk->mat_BAcan_to_BA0_two), &vec_chall);

    // reducing mod 2^len_chall
    ibz_mod(&vec[0], &vec[0], &pow_chall);
    ibz_mod(&vec[1], &vec[1], &pow_chall);

    // lideal_chall_two is the pullback of the ideal challenge through the secret key ideal
    id2iso_kernel_dlogs_to_ideal_two(&lideal_chall_two, &vec, len_chall);
    assert(ibz_cmp(&lideal_chall_two.norm, &pow_chall) == 0);

    // lideal_chall_secret = lideal_secret * lideal_chall_two
    quat_lideal_inter(
        &lideal_chall_secret, &lideal_chall_two, &(sk->secret_ideal), &QUATALG_PINFTY);

    // now we compute lideal_com_to_chall which is dual(Icom)* lideal_chall_secret
    quat_lideal_conjugate_lattice(&lat_commit, &lideal_commit);
    quat_lattice_intersect(&lattice_hom_chall_to_com, &lideal_chall_secret.lattice, &lat_commit);

    // sampling the response
    // TODO : right now the sampling is not done as explained in the paper
    // TODO try with the frobenius conjugate
    ibz_mul(&lattice_content, &(lideal_chall_secret.norm), &(lideal_commit.norm));
    found = sample_response(&resp_quat, &lattice_hom_chall_to_com, &lattice_content, verbose);

    // TODO when it fails, we don't finalize all the ibz
    if (!found) {
        printf("heuristic response sampling failed \n");
        return 0;
    }

    assert(quat_lattice_contains(NULL, &MAXORD_O0, &resp_quat, &QUATALG_PINFTY));
    assert(quat_alg_is_primitive(&resp_quat, &MAXORD_O0, &QUATALG_PINFTY));

    // creating lideal_com * lideal_resp
    // we first compute the norm of lideal_resp
    // norm of the resp_quat
    quat_alg_norm(&temp_norm, &resp_quat, &QUATALG_PINFTY);
    // dividing by n(lideal_com) * n(lideal_secret_chall)
    int is_int = ibq_to_ibz(&degree_full_resp, &temp_norm);
    assert(is_int);
    ibz_div(&degree_full_resp, &remain, &degree_full_resp, &lattice_content);
    assert(ibz_cmp(&remain, &ibz_const_zero) == 0);

    // computing the diadic valuation
    // right now we make the overwhelmingly likely assumption that the diadic valuation of
    // degree_full_resp is smaller than 60
    exp_diadic_val_full_resp = two_adic_valuation(ibz_get(&degree_full_resp));
    assert(exp_diadic_val_full_resp < 60);
    // removing the power of two part
    ibz_pow(&tmp, &ibz_const_two, exp_diadic_val_full_resp);
    ibz_div(&degree_odd_resp, &remain, &degree_full_resp, &tmp);
    assert(ibz_cmp(&remain, &ibz_const_zero) == 0);

    // taking the conjugate of quat_resp so that quat_resp is contained in lideal_com
    quat_alg_conj(&resp_quat, &resp_quat);

    // now we compute the ideal_aux
    // computing the norm
    // TODO make a clean constant for this
    pow_dim2_deg_resp = SQIsign2D_response_heuristic_bound - exp_diadic_val_full_resp;
    ibz_pow(&remain, &ibz_const_two, pow_dim2_deg_resp);
    ibz_sub(&tmp, &remain, &degree_odd_resp);
    assert(ibz_cmp(&tmp, &ibz_const_zero) > 0);

    // sampling the ideal at random
    // TODO replace these two steps with a clean function that samples random ideals from a right
    // order
    sampling_random_ideal_O0(&lideal_aux, &tmp, 0);

    // pushing forward
    quat_lideal_inter(&lideal_aux_com, &lideal_commit, &lideal_aux, &QUATALG_PINFTY);

    // now we evaluate this isogeny on the basis of E0
    // Baux0 = image through aux_com isogeny (odd degree) of canonical basis of E0
    dim2id2iso_arbitrary_isogeny_evaluation(&Baux0, &E_aux, &lideal_aux_com);

#ifndef NDEBUG
    // testing
    assert(test_point_order_twof(&Baux0.P, &E_aux, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Baux0.Q, &E_aux, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Baux0.PmQ, &E_aux, TORSION_PLUS_EVEN_POWER));
#endif

    // now it only remains to format the response for the verification

    // this will be the image of the basis of E0 through  phi_sec quat_resp / deg phi_sec
    // where phi_sec is the isogeny corresponding to sk->lideal_secret
    ec_basis_t bas_sk;

    // we compute the matrix corresponding to resp_quat
    quat_alg_make_primitive(&coeffs, &lattice_content, &resp_quat, &MAXORD_O0, &QUATALG_PINFTY);
    assert(ibz_get(&lattice_content) % 2 == 1);
    ibz_set(&mat[0][0], 0);
    ibz_set(&mat[0][1], 0);
    ibz_set(&mat[1][0], 0);
    ibz_set(&mat[1][1], 0);

    // computing the matrix
    for (unsigned i = 0; i < 2; ++i) {
        ibz_add(&mat[i][i], &mat[i][i], &coeffs[0]);
        for (unsigned j = 0; j < 2; ++j) {
            ibz_mul(&tmp, &ACTION_GEN2[i][j], &coeffs[1]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &ACTION_GEN3[i][j], &coeffs[2]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&tmp, &ACTION_GEN4[i][j], &coeffs[3]);
            ibz_add(&mat[i][j], &mat[i][j], &tmp);
            ibz_mul(&mat[i][j], &mat[i][j], &lattice_content);
            ibz_mod(&mat[i][j], &mat[i][j], &TORSION_PLUS_2POWER);
        }
    }

    // now we inverse the matrix of the secret key to get from the canonical basis of pk to phi(B0)
    ibz_2x2_inv_mod(&sig_mat_pk_can_to_B_pk, &sk->mat_BAcan_to_BA0_two, &TORSION_PLUS_2POWER);

#ifndef NDEBUG
    ec_basis_t bas_test, bas_ref;
    copy_point(&bas_ref.P, &BASIS_EVEN.P);
    copy_point(&bas_ref.Q, &BASIS_EVEN.Q);
    copy_point(&bas_ref.PmQ, &BASIS_EVEN.PmQ);
    copy_point(&bas_test.P, &sk->canonical_basis.P);
    copy_point(&bas_test.Q, &sk->canonical_basis.Q);
    copy_point(&bas_test.PmQ, &sk->canonical_basis.PmQ);
    matrix_application_even_basis(
        &bas_test, &sk->curve, &sig_mat_pk_can_to_B_pk, TORSION_PLUS_EVEN_POWER);
    fp2_t w0, w1, w0_test;
    ec_point_t A24;
    AC_to_A24(&A24, &CURVE_E0);
    ec_normalize_point(&A24);
    weil(&w0, TORSION_PLUS_EVEN_POWER, &bas_ref.P, &bas_ref.Q, &bas_ref.PmQ, &A24);
    weil(&w1, TORSION_PLUS_EVEN_POWER, &bas_test.P, &bas_test.Q, &bas_test.PmQ, &sk->curve.A24);
    digit_t scal[NWORDS_ORDER] = { 0 };
    ibz_to_digit_array(scal, &sk->secret_ideal.norm);
    fp2_pow_vartime(&w0_test, &w0, scal, NWORDS_ORDER);
    assert(fp2_is_equal(&w0_test, &w1));
#endif

    // and we multiply it with the matrix corresponding to resp_quat precomputed above
    // ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk,&mat,&sig_mat_pk_can_to_B_pk,&TORSION_PLUS_2POWER);
    ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk, &sig_mat_pk_can_to_B_pk, &mat, &TORSION_PLUS_2POWER);

    // canonical basis
    ec_basis_t B_aux_can;
    ec_curve_to_basis_2f_to_hint(&B_aux_can, &E_aux, TORSION_PLUS_EVEN_POWER, sig->hint_aux);

    // compute the matrix to go from B_aux0 to B_aux_can
    change_of_basis_matrix_two(
        &mat_Baux0_to_Baux_can, &B_aux_can, &Baux0, &E_aux, TORSION_PLUS_EVEN_POWER);

    // apply the change of basis to the matrix
    // ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk,&mat_Baux0_to_Baux_can,&sig_mat_pk_can_to_B_pk,&TORSION_PLUS_2POWER);
    ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk,
                    &sig_mat_pk_can_to_B_pk,
                    &mat_Baux0_to_Baux_can,
                    &TORSION_PLUS_2POWER);

    // dividing the matrix by the degree of the secret key isogeny
    ibz_copy(&tmp, &sk->secret_ideal.norm);
    assert(ibz_get(&tmp) % 2 == 1);
    ibz_invmod(&tmp, &tmp, &TORSION_PLUS_2POWER);
    ibz_mul(&sig_mat_pk_can_to_B_pk[0][0], &sig_mat_pk_can_to_B_pk[0][0], &tmp);
    ibz_mul(&sig_mat_pk_can_to_B_pk[1][0], &sig_mat_pk_can_to_B_pk[1][0], &tmp);
    ibz_mul(&sig_mat_pk_can_to_B_pk[0][1], &sig_mat_pk_can_to_B_pk[0][1], &tmp);
    ibz_mul(&sig_mat_pk_can_to_B_pk[1][1], &sig_mat_pk_can_to_B_pk[1][1], &tmp);

    // filling the output
    sig->two_resp_length = exp_diadic_val_full_resp;
    // setting sig->E_aux
    fp2_t temp_fp2;
    fp2_copy(&temp_fp2, &E_aux.C);
    fp2_inv(&temp_fp2);
    fp2_mul(&sig->E_aux.A, &temp_fp2, &E_aux.A);
    fp2_set_one(&sig->E_aux.C);
    ec_point_init(&sig->E_aux.A24);
    sig->E_aux.is_A24_computed_and_normalized = 0;

    ibz_pow(&tmp, &ibz_const_two, len_chall + exp_diadic_val_full_resp);
    ibz_pow(&pow_chall,
            &ibz_const_two,
            TORSION_PLUS_EVEN_POWER - (len_chall + exp_diadic_val_full_resp));

    // formatting the challenge info
    if (ibz_get(&sig_mat_pk_can_to_B_pk[0][0]) % 2 != 0) {
        // in that case
        // we will be able to express the kernel as [2^*](P + [x] Q)
        // where P,Q is the canonical basis of Epk
        sig->hint_b = 0;
        ibz_copy(&sig->x, &sig_mat_pk_can_to_B_pk[0][0]);
        ibz_invmod(&sig->x, &sig->x, &tmp);
        ibz_mul(&sig->x, &sig->x, &sig_mat_pk_can_to_B_pk[1][0]);
        ibz_mod(&sig->x, &sig->x, &tmp);
    } else if (ibz_get(&sig_mat_pk_can_to_B_pk[1][0]) % 2 != 0) {
        // in that case
        // we will be able to express the kernel as [2^*](Q + [x] P)
        // where P,Q is the canonical basis of Epk
        sig->hint_b = 1;
        ibz_copy(&sig->x, &sig_mat_pk_can_to_B_pk[1][0]);
        ibz_invmod(&sig->x, &sig->x, &tmp);
        ibz_mul(&sig->x, &sig->x, &sig_mat_pk_can_to_B_pk[0][0]);
        ibz_mod(&sig->x, &sig->x, &tmp);
    } else if (ibz_get(&sig_mat_pk_can_to_B_pk[0][1]) % 2 != 0) {
        // in that case
        // we will be able to express the kernel as [2^*](P + [x] Q)
        // where P,Q is the canonical basis of Epk
        sig->hint_b = 0;
        ibz_copy(&sig->x, &sig_mat_pk_can_to_B_pk[0][1]);
        ibz_invmod(&sig->x, &sig->x, &tmp);
        ibz_mul(&sig->x, &sig->x, &sig_mat_pk_can_to_B_pk[1][1]);
        ibz_mod(&sig->x, &sig->x, &tmp);
    } else {
        assert(ibz_get(&sig_mat_pk_can_to_B_pk[1][1]) % 2 != 0);
        // in that case
        // we will be able to express the kernel as [2^*](Q + [x] P)
        // where P,Q is the canonical basis of Epk
        sig->hint_b = 1;
        ibz_copy(&sig->x, &sig_mat_pk_can_to_B_pk[1][1]);
        ibz_invmod(&sig->x, &sig->x, &tmp);
        ibz_mul(&sig->x, &sig->x, &sig_mat_pk_can_to_B_pk[0][1]);
        ibz_mod(&sig->x, &sig->x, &tmp);
    }

    // now we compute b0,d0
    if (!sig->hint_b) {
        // mat_pk_can_to_B_pk[0][0] = b0 + 2^n b1;
        // mat_pk_can_to_B_pk[1][0] = c0 + 2^n c1;
        // mat_pk_can_to_B_pk[0][1] = d0 + 2^n d1;
        // mat_pk_can_to_B_pk[1][1] = e0 + 2^n e1;
        // with n = TORSION_PLUS_EVEN_POWER - len_chall - exp_diadic_val_full_resp
        // and a = len_chall + exp_diadic_val_full_resp
        // and we set
        // sig->b0 = b0
        // sig->d0 = d0
        // sig->b1 = (c1 - xb1) mod 2^a
        // sig->d1 = (e1 - xd1) mod 2^a
        // when a < n we also compute
        // sig -> c0_adjust = (c0 - (x*b0 mod 2^a))/2^(a)
        // sig -> e0_adjust = (e0 - (x*d0 mod 2^a))/2^(a)

        // computation of sig->b0,sig->d0
        ibz_mod(&sig->b0, &sig_mat_pk_can_to_B_pk[0][0], &pow_chall);
        ibz_mod(&sig->d0, &sig_mat_pk_can_to_B_pk[0][1], &pow_chall);

        // computation of b1
        ibz_sub(&degree_odd_resp, &sig_mat_pk_can_to_B_pk[0][0], &sig->b0);
        ibz_div(&sig->b1, &remain, &degree_odd_resp, &pow_chall);
        assert(ibz_cmp(&remain, &ibz_const_zero) == 0);
        // computation of c0
        ibz_mod(&remain, &sig_mat_pk_can_to_B_pk[1][0], &pow_chall);

        if (len_chall + exp_diadic_val_full_resp <=
            TORSION_PLUS_EVEN_POWER - (len_chall + exp_diadic_val_full_resp)) {
            // computation of sig->c0_adjust
            ibz_mul(&sig->c0_adjust, &sig->x, &sig->b0);
            ibz_mod(&sig->c0_adjust, &sig->c0_adjust, &tmp);
            ibz_sub(&sig->c0_adjust, &remain, &sig->c0_adjust);
            ibz_div(&sig->c0_adjust, &degree_odd_resp, &sig->c0_adjust, &tmp);
            assert(ibz_cmp(&degree_odd_resp, &ibz_const_zero) == 0);
        } else {
            ibz_set(&sig->c0_adjust, 0);
            ibz_set(&sig->e0_adjust, 0);
        }

        // computation of sig->b1 = (c1 - x*b1) mod 2^n
        ibz_sub(&degree_odd_resp, &sig_mat_pk_can_to_B_pk[1][0], &remain);
        ibz_div(&degree_odd_resp, &remain, &degree_odd_resp, &pow_chall); // this is c1
        assert(ibz_cmp(&remain, &ibz_const_zero) == 0);
        ibz_mul(&sig->b1, &sig->b1, &sig->x);
        ibz_sub(&sig->b1, &degree_odd_resp, &sig->b1);
        ibz_mod(&sig->b1, &sig->b1, &tmp); // reducing mod 2^a

        // computation of d1
        ibz_sub(&degree_odd_resp, &sig_mat_pk_can_to_B_pk[0][1], &sig->d0);
        ibz_div(&sig->d1, &remain, &degree_odd_resp, &pow_chall);
        assert(ibz_cmp(&remain, &ibz_const_zero) == 0);
        // computation of e0
        ibz_mod(&remain, &sig_mat_pk_can_to_B_pk[1][1], &pow_chall);

        if (len_chall + exp_diadic_val_full_resp <=
            TORSION_PLUS_EVEN_POWER - (len_chall + exp_diadic_val_full_resp)) {
            // computation of sig->e0_adjust
            ibz_mul(&sig->e0_adjust, &sig->x, &sig->d0);
            ibz_mod(&sig->e0_adjust, &sig->e0_adjust, &tmp);
            ibz_sub(&sig->e0_adjust, &remain, &sig->e0_adjust);
            ibz_div(&sig->e0_adjust, &degree_odd_resp, &sig->e0_adjust, &tmp);
            assert(ibz_cmp(&degree_odd_resp, &ibz_const_zero) == 0);
        }

        // computation of sig->d1 = (e1 - x*d1)
        ibz_sub(&degree_odd_resp, &sig_mat_pk_can_to_B_pk[1][1], &remain);
        ibz_div(&degree_odd_resp, &remain, &degree_odd_resp, &pow_chall); // this is e1
        assert(ibz_cmp(&remain, &ibz_const_zero) == 0);
        ibz_mul(&sig->d1, &sig->d1, &sig->x);
        ibz_sub(&sig->d1, &degree_odd_resp, &sig->d1);
        ibz_mod(&sig->d1, &sig->d1, &tmp); // reducing mod 2^a

#ifndef NDEBUG
        ibz_t I1, I2, I3;
        ibz_init(&I1);
        ibz_init(&I2);
        ibz_init(&I3);
        ibz_mod(&I1, &sig_mat_pk_can_to_B_pk[1][0], &pow_chall);
        ibz_mul(&I2, &sig->b0, &sig->x);
        if (len_chall + exp_diadic_val_full_resp <=
            TORSION_PLUS_EVEN_POWER - (len_chall + exp_diadic_val_full_resp)) {
            ibz_mod(&I2, &I2, &tmp);
            ibz_mul(&I3, &tmp, &sig->c0_adjust);
            ibz_add(&I2, &I2, &I3);

        } else {
            ibz_mod(&I2, &I2, &pow_chall);
        }
        assert(ibz_cmp(&I1, &I2) == 0);
        ibz_mod(&I1, &sig_mat_pk_can_to_B_pk[1][1], &pow_chall);
        ibz_mul(&I2, &sig->d0, &sig->x);
        if (len_chall + exp_diadic_val_full_resp <=
            TORSION_PLUS_EVEN_POWER - (len_chall + exp_diadic_val_full_resp)) {
            ibz_mod(&I2, &I2, &tmp);
            ibz_mul(&I3, &tmp, &sig->e0_adjust);
            ibz_add(&I2, &I2, &I3);
        } else {
            ibz_mod(&I2, &I2, &pow_chall);
        }
        assert(ibz_cmp(&I1, &I2) == 0);
        ibz_finalize(&I1);
        ibz_finalize(&I2);
        ibz_finalize(&I3);
#endif
    } else {
        // TODECIDE if we need to treat this case
        assert(0);
    }

    ibz_finalize(&pow_chall);
    ibz_vec_2_finalize(&vec);
    ibz_vec_2_finalize(&vec_chall);
    ibz_vec_2_finalize(&vec_resp_two);
    ibz_mat_2x2_finalize(&mat_Bchall_can_to_Bchall);
    ibz_mat_2x2_finalize(&mat_Baux0_to_Baux_can);
    ibz_mat_2x2_finalize(&sig_mat_pk_can_to_B_pk);
    ibz_mat_2x2_finalize(&mat);

    quat_alg_elem_finalize(&resp_quat);
    quat_alg_elem_finalize(&elem_tmp);
    quat_lattice_finalize(&lattice_hom_chall_to_com);
    quat_lattice_finalize(&lat_commit);
    quat_left_ideal_finalize(&lideal_commit);
    quat_left_ideal_finalize(&lideal_chall_two);
    quat_left_ideal_finalize(&lideal_chall_secret);
    quat_left_ideal_finalize(&lideal_resp_two);
    quat_left_ideal_finalize(&lideal_tmp);
    quat_left_ideal_finalize(&lideal_com_resp);
    quat_left_ideal_finalize(&lideal_aux_com);
    quat_left_ideal_finalize(&lideal_aux);

    ibz_vec_4_finalize(&coeffs);
    ibz_finalize(&degree_full_resp);
    ibz_finalize(&degree_odd_resp);

    ibq_finalize(&temp_norm);
    ibz_finalize(&tmp);
    ibz_finalize(&lattice_content);
    ibz_finalize(&remain);

    return found;
}

int
protocols_verif(signature_t *sig, const public_key_t *pk, const unsigned char *m, size_t l)
{

    int verif;
    ibz_t tmp, tmp2, remain;
    ibz_vec_2_t vec_chall, check_vec_chall;
    ibz_mat_2x2_t mat;
    ec_isog_even_t phi_chall;
    ec_basis_t bas_EA, B_chall;
    ec_curve_t Epk;

    ec_curve_init(&Epk);
    ibz_mat_2x2_init(&mat);
    ibz_init(&tmp);
    ibz_init(&remain);
    ibz_init(&tmp2);
    ibz_vec_2_init(&vec_chall);
    ibz_vec_2_init(&check_vec_chall);

    // checking that we are given A coefficients and no precomputation
    assert(fp2_is_one(&pk->curve.C) && !pk->curve.is_A24_computed_and_normalized);
    assert(fp2_is_one(&sig->E_aux.C) && !sig->E_aux.is_A24_computed_and_normalized);
    clock_t t = tic();

    // TODO use a constant for this
    int len_chall = SQIsign2D_heuristic_challenge_length;

    copy_curve(&Epk, &pk->curve);
    phi_chall.curve = Epk;
    phi_chall.length = len_chall + sig->two_resp_length;

    // creating a matrix to apply to the basis

    // mat [0][0] = sig->b0
    // mat [0][1] = sig->d0
    ibz_copy(&mat[0][0], &sig->b0);
    ibz_copy(&mat[0][1], &sig->d0);
    if (phi_chall.length <= TORSION_PLUS_EVEN_POWER - phi_chall.length) {
        // mat [1][0] = (sig->b0 * sig->x mod 2^a) + 2^a * sig->cO_adjust + 2^n * sig->b1
        // mat [1][1] = (sig->d0 * sig->x mod 2^a) + 2^a * sig->e0_adjust + 2^n * sig->d1
        ibz_pow(&tmp,
                &ibz_const_two,
                TORSION_PLUS_EVEN_POWER - len_chall - sig->two_resp_length); // 2^n
        ibz_mul(&mat[1][0], &tmp, &sig->b1);
        ibz_mul(&mat[1][1], &tmp, &sig->d1);

        ibz_pow(&tmp, &ibz_const_two, len_chall + sig->two_resp_length); // 2^a
        ibz_mul(&remain, &sig->b0, &sig->x);
        ibz_mod(&remain, &remain, &tmp);
        ibz_mul(&tmp2, &tmp, &sig->c0_adjust);
        ibz_add(&remain, &remain, &tmp2);
        ibz_add(&mat[1][0], &mat[1][0], &remain);
        ibz_mul(&remain, &sig->d0, &sig->x);
        ibz_mod(&remain, &remain, &tmp);
        ibz_mul(&tmp2, &tmp, &sig->e0_adjust);
        ibz_add(&remain, &remain, &tmp2);
        ibz_add(&mat[1][1], &mat[1][1], &remain);
    } else {
        // mat [1][0] = (sig->b0 * sig->x mod 2^n) + 2^n * sig->b1
        // mat [1][1] = (sig->d0 * sig->x mod 2^n) + 2^n * sig->d1
        ibz_pow(&tmp,
                &ibz_const_two,
                TORSION_PLUS_EVEN_POWER - len_chall - sig->two_resp_length); // 2^n
        ibz_mul(&mat[1][0], &tmp, &sig->b1);
        ibz_mul(&mat[1][1], &tmp, &sig->d1);
        ibz_mul(&remain, &sig->b0, &sig->x);
        ibz_mod(&remain, &remain, &tmp);
        ibz_add(&mat[1][0], &mat[1][0], &remain);
        ibz_mul(&remain, &sig->d0, &sig->x);
        ibz_mod(&remain, &remain, &tmp);
        ibz_add(&mat[1][1], &mat[1][1], &remain);
    }

    // computation of the challenge
    // canonical basis
    ec_curve_to_basis_2f_from_hint(
        &bas_EA, &Epk, TORSION_PLUS_EVEN_POWER, pk->hint_pk); // canonical

    assert(test_point_order_twof(&bas_EA.P, &Epk, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&bas_EA.Q, &Epk, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&bas_EA.PmQ, &Epk, TORSION_PLUS_EVEN_POWER));

    t = tic();

    // applying the change matrix on the basis of E_pk
    matrix_application_even_basis(&bas_EA, &Epk, &mat, TORSION_PLUS_EVEN_POWER);

    // TOC_clock(t,"matrix application");

    t = tic();

    // TODO make a clean constant for this
    int pow_dim2_deg_resp = SQIsign2D_response_heuristic_bound - sig->two_resp_length;
    assert(TORSION_PLUS_EVEN_POWER == phi_chall.length + pow_dim2_deg_resp);

    // computing the isogeny
    // computing the small two chain
    ec_curve_t Echall;
    copy_curve(&Echall, &Epk);
    ec_point_t points[3];

    // choosing the right point for the small two_isogenies
    if (ibz_get(&mat[0][0]) % 2 == 0 && ibz_get(&mat[1][0]) % 2 == 0) {
        copy_point(&phi_chall.kernel, &bas_EA.Q);
    } else {
        copy_point(&phi_chall.kernel, &bas_EA.P);
    }

    copy_point(&points[0], &bas_EA.P);
    copy_point(&points[1], &bas_EA.Q);
    copy_point(&points[2], &bas_EA.PmQ);
    ec_dbl_iter(&phi_chall.kernel,
                TORSION_PLUS_EVEN_POWER - (phi_chall.length),
                &Echall,
                &phi_chall.kernel);

    assert(test_point_order_twof(&phi_chall.kernel, &Echall, phi_chall.length));
    ec_eval_even(&Echall, &phi_chall, points, 3);

    assert(ec_is_on_curve(&Echall, &points[0]));
    assert(ec_is_on_curve(&Echall, &points[1]));
    assert(ec_is_on_curve(&Echall, &points[2]));

    assert(test_point_order_twof(&points[0], &Echall, pow_dim2_deg_resp));
    assert(test_point_order_twof(&points[1], &Echall, pow_dim2_deg_resp));
    assert(test_point_order_twof(&points[2], &Echall, pow_dim2_deg_resp));
    copy_point(&B_chall.P, &points[0]);
    copy_point(&B_chall.Q, &points[1]);
    copy_point(&B_chall.PmQ, &points[2]);

    // TOC_clock(t,"challenge isogeny computation");

    t = tic();

    ec_basis_t B_aux_can;

    // recovering the canonical basis
    ec_curve_to_basis_2f_from_hint(&B_aux_can, &sig->E_aux, TORSION_PLUS_EVEN_POWER, sig->hint_aux);

    // setting to the right order
    ec_dbl_iter(
        &B_aux_can.P, TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp, &sig->E_aux, &B_aux_can.P);
    ec_dbl_iter(
        &B_aux_can.Q, TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp, &sig->E_aux, &B_aux_can.Q);
    ec_dbl_iter(
        &B_aux_can.PmQ, TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp, &sig->E_aux, &B_aux_can.PmQ);

#ifndef NDEBUG
    ec_basis_t bas_test, bas_ref;
    copy_point(&bas_ref.P, &B_chall.P);
    copy_point(&bas_ref.Q, &B_chall.Q);
    copy_point(&bas_ref.PmQ, &B_chall.PmQ);
    copy_point(&bas_test.P, &B_aux_can.P);
    copy_point(&bas_test.Q, &B_aux_can.Q);
    copy_point(&bas_test.PmQ, &B_aux_can.PmQ);
    fp2_t w0, w1, w0_test;
    ec_point_t A24;
    AC_to_A24(&A24, &Echall);
    ec_normalize_point(&A24);
    weil(&w0, pow_dim2_deg_resp, &bas_ref.P, &bas_ref.Q, &bas_ref.PmQ, &A24);
    AC_to_A24(&A24, &sig->E_aux);
    ec_normalize_point(&A24);
    weil(&w1, pow_dim2_deg_resp, &bas_test.P, &bas_test.Q, &bas_test.PmQ, &A24);
    fp2_mul(&w0_test, &w0, &w1);
    assert(!fp2_is_one(&w0));
    assert(fp2_is_one(&w0_test));

#endif

    // now compute the dim2 isogeny from E_chall_2 x E_aux -> E_com x E_aux'
    // of kernel B_chall_can x B_aux_can

    // first we set-up the kernel
    theta_couple_curve_t EchallxEaux;
    theta_couple_point_t T1, T2, T1m2;
    theta_chain_t isog;
    copy_curve(&EchallxEaux.E1, &Echall);
    copy_curve(&EchallxEaux.E2, &sig->E_aux);
    copy_point(&T1.P2, &B_aux_can.P);
    copy_point(&T2.P2, &B_aux_can.Q);
    copy_point(&T1m2.P2, &B_aux_can.PmQ);
    copy_point(&T1.P1, &B_chall.P);
    copy_point(&T2.P1, &B_chall.Q);
    copy_point(&T1m2.P1, &B_chall.PmQ);

    // computing the isogeny
    // no points above the twotorsion
    int extra_info = 0;
    theta_chain_comput_strategy_faster_no_eval(
        &isog,
        pow_dim2_deg_resp,
        &EchallxEaux,
        &T1,
        &T2,
        &T1m2,
        strategies[TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp + 2],
        extra_info);

    // TOC_clock(t,"response isogeny");

    // computing the commitment curve
    // apparently its always the second one
    ec_curve_t E_com;
    copy_curve(&E_com, &isog.codomain.E2);

    // recomputing the challenge vector
    hash_to_challenge(&check_vec_chall, &E_com, m, pk, l);

    // comparing
    // recovering the exact vec_chall
    if (sig->hint_b) {
        ibz_copy(&vec_chall[0], &sig->x);
        ibz_set(&vec_chall[1], 1);
    } else {
        ibz_copy(&vec_chall[1], &sig->x);
        ibz_set(&vec_chall[0], 1);
    }
    ibz_pow(&tmp, &ibz_const_two, len_chall);

    // performing the final check
    if (sig->hint_b) {
        ibz_mul(&vec_chall[0], &vec_chall[0], &check_vec_chall[1]);
        ibz_mod(&vec_chall[0], &vec_chall[0], &tmp);
        ibz_mod(&check_vec_chall[0], &check_vec_chall[0], &tmp);
        verif = (ibz_cmp(&vec_chall[0], &check_vec_chall[0]) == 0);
    } else {
        ibz_mul(&vec_chall[1], &vec_chall[1], &check_vec_chall[0]);
        ibz_mod(&vec_chall[1], &vec_chall[1], &tmp);
        ibz_mod(&check_vec_chall[1], &check_vec_chall[1], &tmp);
        verif = (ibz_cmp(&vec_chall[1], &check_vec_chall[1]) == 0);
    }
    // if it didn't work, we start again with the other codomain
    if (!verif) {
        // trying with the other curve
        copy_curve(&E_com, &isog.codomain.E1);
        // recomputing the challenge vector
        hash_to_challenge(&check_vec_chall, &E_com, m, pk, l);

        // performing the final check
        if (sig->hint_b) {
            ibz_mul(&vec_chall[0], &vec_chall[0], &check_vec_chall[1]);
            ibz_mod(&vec_chall[0], &vec_chall[0], &tmp);
            ibz_mod(&check_vec_chall[0], &check_vec_chall[0], &tmp);
            verif = (ibz_cmp(&vec_chall[0], &check_vec_chall[0]) == 0);
        } else {

            ibz_mul(&vec_chall[1], &vec_chall[1], &check_vec_chall[0]);
            ibz_mod(&vec_chall[1], &vec_chall[1], &tmp);
            ibz_mod(&check_vec_chall[1], &check_vec_chall[1], &tmp);
            verif = (ibz_cmp(&vec_chall[1], &check_vec_chall[1]) == 0);
        }
    }

    ibz_finalize(&remain);
    ibz_finalize(&tmp);
    ibz_finalize(&tmp2);
    ibz_vec_2_finalize(&vec_chall);
    ibz_vec_2_finalize(&check_vec_chall);
    ibz_mat_2x2_finalize(&mat);

    return verif;
}
