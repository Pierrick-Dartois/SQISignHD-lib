#include <sqisigndim2.h>
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
    ibz_mat_2x2_init(&(sig->mat_Bchall_can_to_B_chall));
    ibz_init(&sig->chall_coeff);
    sig->hint_aux = (int *)malloc(2 * sizeof(int));
    sig->hint_chall = (int *)malloc(2 * sizeof(int));
}

void
secret_sig_finalize(signature_t *sig)
{
    ibz_mat_2x2_finalize(&(sig->mat_Bchall_can_to_B_chall));
    ibz_finalize(&sig->chall_coeff);
    free(sig->hint_aux);
    free(sig->hint_chall);
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

void
print_signature(const signature_t *sig)
{
    fp2_t j;
    ec_j_inv(&j, &sig->E_aux);
    fp2_print("j_E1 = ", &j);
    // ibz_mat_2x2_print(&sig->mat_sigma_phichall);
    ibz_printf("M_sigma[00] = %Zd, ", &((sig->mat_Bchall_can_to_B_chall)[0][0]));
    ibz_printf("M_sigma[01] = %Zd, ", &((sig->mat_Bchall_can_to_B_chall)[0][1]));
    ibz_printf("M_sigma[10] = %Zd, ", &((sig->mat_Bchall_can_to_B_chall)[1][0]));
    ibz_printf("M_sigma[11] = %Zd", &((sig->mat_Bchall_can_to_B_chall)[1][1]));
}

void
print_public_key(const public_key_t *pk)
{
    fp2_t j;
    ec_j_inv(&j, &pk->curve);
    fp2_print("j_EA = ", &j);
}

// compute the commitment with ideal to isogeny clapotis
// and apply it to the basis of E0 (together with the multiplication by some scalar u)
// the scalar adjusting_factor is a scalar through which the points of the basis are multiplied
void
commit(ec_curve_t *E_com, ec_basis_t *basis_even_com, quat_left_ideal_t *lideal_com)
{

    int found = 1;
    ibz_t n;
    ec_basis_t B_0_two;
    ibz_init(&n);

    // generate a random ideal of random norm for the secret ideal
    generate_random_prime(&n, 1, ibz_bitsize(&QUATALG_PINFTY.p) );
    sampling_random_ideal_O0(lideal_com, &n, 1);

    // ideal to isogeny clapotis
    found = dim2id2iso_arbitrary_isogeny_evaluation(basis_even_com, E_com, lideal_com);
    assert(found);

    ibz_finalize(&n);
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

int
is_good_norm(ibz_t *N)
{
    if ((8 - (ibz_get(N) % 8)) != 5)
        return 0;

    ibz_t pow2, sum_of_squares_candidate;
    ibz_init(&pow2);
    ibz_init(&sum_of_squares_candidate);
    int res = 0;

    ibz_set(&pow2, 1);
    ibz_mul_2exp(&pow2, &pow2, RESPONSE_LENGTH);

    if (ibz_cmp(&pow2, N) < 0) {
        ibz_printf(
            "WARNING: short vectors not short enough...\n2-pow = %Zd\nnorm = %Zd\n", &pow2, &N);
        // assert(0);
        ibz_finalize(&sum_of_squares_candidate);
        return 0;
    }

    ibz_sub(&sum_of_squares_candidate, &pow2, N);

    // unsigned int N_mod_four = ibz_mod_ui (&N, 4);
    assert(ibz_mod_ui(&sum_of_squares_candidate, 8) == 5);

    // if (N_mod_four == 1) {
    res = ibz_probab_prime(&sum_of_squares_candidate, 40);

    ibz_finalize(&pow2);
    ibz_finalize(&sum_of_squares_candidate);
    return res;
}

int
is_good(quat_alg_elem_t *x, ibz_t const *lattice_content)
{
    ibq_t N_q;
    ibz_t N, tmp, pow2;
    ibq_init(&N_q);
    ibz_init(&N);
    ibz_init(&tmp);
    int res = 0;

    quat_alg_norm(&N_q, x, &QUATALG_PINFTY);
    ibq_to_ibz(&N, &N_q);

    // ibz_printf(">>>> %Zd | %Zd\n", lattice_content, &N);

    assert(ibz_divides(&N, lattice_content));

    ibz_div(&N, &tmp, &N, lattice_content);

    res = is_good_norm(&N);

#ifndef NDEBUG
    ibz_init(&pow2);
    ibz_set(&pow2, 1);
    ibz_mul_2exp(&pow2, &pow2, RESPONSE_LENGTH);
    int res2 = 0;
    if (ibz_cmp(&pow2, &N) < 0) {
        ibz_printf(
            "WARNING: short vectors not short enough...\n2-pow = %Zd\nnorm = %Zd\n", &pow2, &N);
    } else {
        ibz_sub(&N, &pow2, &N);

        // unsigned int N_mod_four = ibz_mod_ui (&N, 4);
        unsigned int N_mod_eight = ibz_mod_ui(&N, 8);

        // if (N_mod_four == 1) {
        if (N_mod_eight == 5) {
            res2 = ibz_probab_prime(&N, 40);
        }
    }

    assert(res == res2);
    ibz_finalize(&pow2);
#endif

    ibz_finalize(&N);
    ibq_finalize(&N_q);
    ibz_finalize(&tmp);
    return res;
}

void
norm_from_2_times_gram(ibz_t *norm, ibz_mat_4x4_t *gram, ibz_vec_4_t *vec)
{
    quat_qf_eval(norm, gram, vec);
    assert(ibz_is_even(norm));
    ibz_div_2exp(norm, norm, 1);
}

// TODECIDE : is the current sampling method satisfactory ?
// it seems uniformish but not sure the distribution is exactly as we might want it
void
sample_response(quat_alg_elem_t *x,
                const quat_lattice_t *lattice,
                ibz_t const *lattice_content,
                int verbose)
{
    ibz_mat_4x4_t lll;
    ibz_t denom_gram, norm;
    ibz_t bound;
    ibz_vec_4_t vec;

    ibz_mat_4x4_init(&lll);
    ibz_init(&denom_gram);
    ibz_init(&norm);
    ibz_init(&bound);
    ibz_vec_4_init(&vec);

    int err = quat_lattice_lll(&lll, lattice, &(QUATALG_PINFTY.p));
    assert(!err);

    // The shortest vector found by lll is our response
    ibz_mat_4x4_t prod, gram;
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_init(&gram);
    //
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

    int found = 0;
    int count = 0;

    ibz_vec_4_t b_bound;
    ibz_vec_4_init(&b_bound);

    ibz_pow(&bound, &ibz_const_two, SQIsign2D_response_length);

    // computing the upperbounds for the coefficients of the scalar decomposition
    int first_zero_index = -1;
    for (int j = 0; j < 4; j++) {
        ibz_copy(&b_bound[j], &gram[j][j]);
        ibz_div_2exp(&b_bound[j], &b_bound[j], 1);
        ibz_div(&b_bound[j], &norm, &bound, &b_bound[j]);
        ibz_sqrt_floor(&b_bound[j], &b_bound[j]);
        if (first_zero_index == -1 && ibz_cmp(&b_bound[j], &ibz_const_zero) == 0) {
            first_zero_index = j;
        }
    }
    if (first_zero_index == -1) {
        first_zero_index = 4;
    }

    // TODO make this a proper constant of the scheme
    // loop to find a correct answer
    while (!found && count < 50) {

        for (int i = 0; i < first_zero_index; i++) {
            ibz_rand_interval_minm_m(&vec[i], ibz_get(&b_bound[i]));
        }
        for (int i = first_zero_index; i < 4; i++) {
            ibz_set(&vec[i], 0);
        }

        norm_from_2_times_gram(&norm, &gram, &vec);

        // checking that we got something small enough
        found = ibz_cmp(&norm, &bound) < 0 &&
                (ibz_cmp(&vec[0], &ibz_const_zero) != 0 || ibz_cmp(&vec[1], &ibz_const_zero) != 0 ||
                 ibz_cmp(&vec[2], &ibz_const_zero) != 0 || ibz_cmp(&vec[3], &ibz_const_zero) != 0);
        if (found) {
            // computing the absolute coordinates of the result
            ibz_mat_4x4_eval(&(x->coord), &lll, &vec);
            assert(quat_lattice_contains(NULL, lattice, x, &QUATALG_PINFTY));
        }

        count++;
    }

    // if not found then we just use the smallest vector of the lattice
    if (!found) {
        for (int i = 0; i < 4; i++) {
            ibz_copy(&x->coord[i], &lll[i][0]);
        }
        assert(quat_lattice_contains(NULL, lattice, x, &QUATALG_PINFTY));
        found = 1;
    }

    ibz_vec_4_finalize(&b_bound);

    ibz_finalize(&denom_gram);
    ibz_finalize(&norm);
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&lll);
    ibz_vec_4_finalize(&vec);
    ibz_finalize(&bound);
    return;
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

    {
        digit_t digits[NWORDS_FIELD];

        // FIXME should use SHAKE128 for smaller parameter sets?
        SHAKE256(
            (void *)digits, sizeof(digits), buf, FP2_ENCODED_BYTES + FP2_ENCODED_BYTES + length);

        ibz_set(&(*scalars)[1], 1);
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

    ibz_t lattice_content;
    ec_curve_t E_aux, E_com;
    ec_basis_t Bcom0, Baux0; // basis of 2^TORSION_PLUS_EVEN_POWER-torsion
    ec_basis_t B_resp_two;
    ibz_vec_4_t dummy_coord;
    ibz_vec_2_t vec, vec_chall, vec_resp_two;
    quat_left_ideal_t lideal_tmp;
    quat_left_ideal_t lideal_commit, lideal_chall_two;
    quat_left_ideal_t lideal_chall_secret, lideal_com_resp, lideal_aux, lideal_aux_resp_com,
        lideal_resp_two;
    quat_lattice_t lattice_hom_chall_to_com, lat_commit;
    quat_alg_elem_t resp_quat;
    quat_alg_elem_t elem_tmp;
    ibz_mat_2x2_t mat_Baux2_to_Baux2_can, mat_Bchall_can_to_Bchall;
    // ibz_mat_2x2_t mat_sigma_phichall_BA_to_Bcomcan, mat_sigma_phichall_BA0_to_Bcom0;
    ibz_t degree_com_isogeny, tmp, remain;
    ibz_t degree_full_resp, degree_odd_resp;
    ibq_t temp_norm;
    int exp_diadic_val_full_resp;
    int pow_dim2_deg_resp;
    int backtracking;

    ec_curve_init(&E_aux);
    ec_curve_init(&E_com);

    ibz_init(&tmp);
    ibz_init(&lattice_content);
    ibz_init(&remain);

    ibz_init(&degree_full_resp);
    ibz_init(&degree_odd_resp);
    ibq_init(&temp_norm);

    ibz_mat_2x2_init(&mat_Bchall_can_to_Bchall);
    ibz_mat_2x2_init(&mat_Baux2_to_Baux2_can);
    ibz_vec_4_init(&dummy_coord);
    // ibz_mat_2x2_init(&mat_sigma_phichall_BA_to_Bcomcan);
    // ibz_mat_2x2_init(&mat_sigma_phichall_BA0_to_Bcom0);

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
    quat_left_ideal_init(&lideal_aux_resp_com);
    quat_left_ideal_init(&lideal_aux);

    ibz_vec_2_init(&vec);
    ibz_vec_2_init(&vec_chall);
    ibz_vec_2_init(&vec_resp_two);

    // computing the commitment
    commit(&E_com, &Bcom0, &lideal_commit);


    // computing the challenge
    // vec_chall is a pair of coefficients encoding the kernel of the challenge isogeny
    // as vec_chall[0]*B[0] + vec_chall[1]*B[1] where B is the canonical basis of the
    // 2^TORSION_PLUS_EVEN_TORSION torsion of EA
    hash_to_challenge(&vec_chall, &E_com, m, pk, l);

    // now we compute the ideal associated to the challenge
    // for that, we need to find vec such that
    // the kernel of the challenge isogeny is generated by vec[0]*B0[0] + vec[1]*B0[1] where B0 is
    // the image through the secret key isogeny of the canonical basis E0
    ibz_mat_2x2_eval(&vec, &(sk->mat_BAcan_to_BA0_two), &vec_chall);

    // lideal_chall_two is the pullback of the ideal challenge through the secret key ideal
    id2iso_kernel_dlogs_to_ideal_two(&lideal_chall_two, &vec, TORSION_PLUS_EVEN_POWER);
    assert(ibz_cmp(&lideal_chall_two.norm, &TORSION_PLUS_2POWER) == 0);

    // lideal_chall_secret = lideal_secret * lideal_chall_two
    quat_lideal_inter(
        &lideal_chall_secret, &lideal_chall_two, &(sk->secret_ideal), &QUATALG_PINFTY);

    // now we compute lideal_com_to_chall which is dual(Icom)* lideal_chall_secret
    quat_lideal_conjugate_lattice(&lat_commit, &lideal_commit);
    quat_lattice_intersect(&lattice_hom_chall_to_com, &lideal_chall_secret.lattice, &lat_commit);

    // sampling the smallest response
    ibz_mul(&lattice_content, &(lideal_chall_secret.norm), &(lideal_commit.norm));
    sample_response(&resp_quat, &lattice_hom_chall_to_com, &lattice_content, verbose);

    // computing the amount of backtracking we're making
    // and removing it
    quat_alg_make_primitive(&dummy_coord, &tmp, &resp_quat, &MAXORD_O0, &QUATALG_PINFTY);
    ibz_mul(&resp_quat.denom, &resp_quat.denom, &tmp);
    assert(quat_lattice_contains(NULL, &MAXORD_O0, &resp_quat, &QUATALG_PINFTY));

    // the backtracking is the common part of the response and the challenge
    // the degree of the backtring is the scalar tmp computed above such that quat_resp is in tmp *
    // O0 we assume that the length of the backtracking is smaller than 60;
    backtracking = two_adic_valuation(ibz_get(&tmp));
    assert(backtracking < SQIsign2D_backtracking_bound);
    ibz_pow(&tmp, &ibz_const_two, backtracking);
    ibz_div(&lattice_content, &remain, &lattice_content, &tmp);

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
#ifndef NDEBUG
    ibz_pow(&tmp, &ibz_const_two, SQIsign2D_response_length);
    assert(ibz_cmp(&tmp, &degree_odd_resp) > 0);
#endif

    // creating the ideal
    quat_alg_conj(&resp_quat, &resp_quat);
    // setting the norm
    ibz_mul(&tmp, &lideal_commit.norm, &degree_odd_resp);
    quat_lideal_create_from_primitive(
        &lideal_com_resp, &resp_quat, &tmp, &MAXORD_O0, &QUATALG_PINFTY);

    // now we compute the ideal_aux
    // computing the norm
    pow_dim2_deg_resp = SQIsign2D_response_length - exp_diadic_val_full_resp;
    ibz_pow(&remain, &ibz_const_two, pow_dim2_deg_resp);
    ibz_sub(&tmp, &remain, &degree_odd_resp);

    // multiplying by 4 to account for the fact that we use the 4 torsion above the kernel
    ibz_mul(&remain, &remain, &ibz_const_two);
    ibz_mul(&remain, &remain, &ibz_const_two);

    // sampling the ideal at random
    // TODO replace these two steps with a clean function that samples random ideals from a right
    // order
    sampling_random_ideal_O0(&lideal_aux, &tmp, 0);
    // pushing forward
    quat_lideal_inter(&lideal_aux_resp_com, &lideal_com_resp, &lideal_aux, &QUATALG_PINFTY);


    // now we evaluate this isogeny on the basis of E0
    dim2id2iso_arbitrary_isogeny_evaluation(&Baux0, &E_aux, &lideal_aux_resp_com);

    // notational conventions:
    // B0 = canonical basis of E0
    // Bcom0 = image through commitment isogeny (odd degree) of canonical basis of E0
    // Baux0 = image through aux_resp_com isogeny (odd degree) of canonical basis of E0

#ifndef NDEBUG
    // testing
    assert(test_point_order_twof(&Bcom0.P, &E_com, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Bcom0.Q, &E_com, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Bcom0.PmQ, &E_com, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Baux0.P, &E_aux, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Baux0.Q, &E_aux, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&Baux0.PmQ, &E_aux, TORSION_PLUS_EVEN_POWER));
#endif


    // applying the matrix to compute Baux
    // first, we copy and reduce to the relevant order
    ec_dbl_iter(&Baux0.P,
                TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp - exp_diadic_val_full_resp - 2,
                &E_aux,
                &Baux0.P);
    ec_dbl_iter(&Baux0.Q,
                TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp - exp_diadic_val_full_resp - 2,
                &E_aux,
                &Baux0.Q);
    ec_dbl_iter(&Baux0.PmQ,
                TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp - exp_diadic_val_full_resp - 2,
                &E_aux,
                &Baux0.PmQ);
    ec_dbl_iter(&Bcom0.P,
                TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp - exp_diadic_val_full_resp - 2,
                &E_com,
                &Bcom0.P);
    ec_dbl_iter(&Bcom0.Q,
                TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp - exp_diadic_val_full_resp - 2,
                &E_com,
                &Bcom0.Q);
    ec_dbl_iter(&Bcom0.PmQ,
                TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp - exp_diadic_val_full_resp - 2,
                &E_com,
                &Bcom0.PmQ);

    // now, we compute the isogeny Phi : Ecom x Eaux -> Echl' x Eaux'
    // where Echl' is 2^exp_diadic_val_full_resp isogenous to Echal
    // ker Phi = <(Bcom_can.P,Baux.P),(Bcom_can.Q,Baux.Q)>
    theta_chain_t isog;
    theta_couple_point_t T1, T2, T1m2;
    theta_couple_curve_t EcomXEaux;
    // preparing the domain
    copy_curve(&EcomXEaux.E1, &E_com);
    copy_curve(&EcomXEaux.E2, &E_aux);

    // preparing the kernel
    copy_point(&T1.P1, &Bcom0.P);
    copy_point(&T2.P1, &Bcom0.Q);
    copy_point(&T1m2.P1, &Bcom0.PmQ);

    copy_point(&T1.P2, &Baux0.P);
    copy_point(&T2.P2, &Baux0.Q);
    copy_point(&T1m2.P2, &Baux0.PmQ);

    // multiplying by 1/ deg resp
    ibz_invmod(&tmp, &degree_odd_resp, &remain);
    ec_mul_ibz(&T1.P2, &E_aux, &tmp, &T1.P2);
    ec_mul_ibz(&T2.P2, &E_aux, &tmp, &T2.P2);
    ec_mul_ibz(&T1m2.P2, &E_aux, &tmp, &T1m2.P2);

    // and multiplying by 2^exp_diadic...
    ec_dbl_iter(&T1.P1, exp_diadic_val_full_resp, &E_com, &T1.P1);
    ec_dbl_iter(&T1.P2, exp_diadic_val_full_resp, &E_aux, &T1.P2);
    ec_dbl_iter(&T2.P1, exp_diadic_val_full_resp, &E_com, &T2.P1);
    ec_dbl_iter(&T2.P2, exp_diadic_val_full_resp, &E_aux, &T2.P2);
    ec_dbl_iter(&T1m2.P1, exp_diadic_val_full_resp, &E_com, &T1m2.P1);
    ec_dbl_iter(&T1m2.P2, exp_diadic_val_full_resp, &E_aux, &T1m2.P2);

    int extra_info = 1;

    // computation of the dim2 isogeny
    theta_chain_comput_strategy(&isog,
                                pow_dim2_deg_resp,
                                &EcomXEaux,
                                &T1,
                                &T2,
                                &T1m2,
                                strategies[TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp],
                                extra_info);

    // pushing the points of torsion to recover the kernel of the dual

    theta_couple_point_t Tev1, Tev2, Tev1m2;

    theta_couple_point_t Teval1, Teval2, Teval3;
    copy_point(&Teval1.P1, &Bcom0.P);
    copy_point(&Teval3.P1, &Bcom0.PmQ);
    copy_point(&Teval2.P1, &Bcom0.Q);
    ec_set_zero(&Teval1.P2);
    ec_set_zero(&Teval2.P2);
    ec_set_zero(&Teval3.P2);
    theta_chain_eval_special_case(&Tev1, &isog, &Teval1, &EcomXEaux);
    theta_chain_eval_special_case(&Tev2, &isog, &Teval2, &EcomXEaux);
    theta_chain_eval_special_case(&Tev1m2, &isog, &Teval3, &EcomXEaux);

    assert(test_point_order_twof(
        &Tev1.P1, &isog.codomain.E1, pow_dim2_deg_resp + exp_diadic_val_full_resp + 2));
    assert(test_point_order_twof(
        &Tev1.P2, &isog.codomain.E2, pow_dim2_deg_resp + exp_diadic_val_full_resp + 2));

    ec_curve_t E_chall_2;
    copy_curve(&E_chall_2, &isog.codomain.E2);

    // copying torsion point, it should always be the second curve
    copy_point(&B_resp_two.P, &Tev1.P2);
    copy_point(&B_resp_two.Q, &Tev2.P2);
    copy_point(&B_resp_two.PmQ, &Tev1m2.P2);

    // computation of the remaining small chain of two isogenies when needed
    if (exp_diadic_val_full_resp > 0) {

        // computing the ideal
        ibz_pow(&tmp, &ibz_const_two, exp_diadic_val_full_resp);

        // we compute the generator of the challenge ideal
        quat_lideal_create_from_primitive(
            &lideal_resp_two, &resp_quat, &tmp, &MAXORD_O0, &QUATALG_PINFTY);

        // computing the coefficients of the kernel in terms of the basis of O0
        ec_isog_even_t phi_two;
        id2iso_ideal_to_isogeny_even_dlogs(&phi_two, &vec_resp_two, &lideal_resp_two);

        // dividing by the right power of 2
        ibz_pow(&tmp, &ibz_const_two, TORSION_PLUS_EVEN_POWER - exp_diadic_val_full_resp);
        ibz_div(&vec_resp_two[0], &remain, &vec_resp_two[0], &tmp);
        assert(ibz_cmp(&remain, &ibz_const_zero) == 0);
        ibz_div(&vec_resp_two[1], &remain, &vec_resp_two[1], &tmp);
        assert(ibz_cmp(&remain, &ibz_const_zero) == 0);

        ec_point_t points[3];
        copy_point(&points[0], &B_resp_two.P);
        copy_point(&points[1], &B_resp_two.Q);
        copy_point(&points[2], &B_resp_two.PmQ);

        // getting down to the right order and applying the matrix
        ec_dbl_iter(&B_resp_two.P, pow_dim2_deg_resp + 2, &isog.codomain.E2, &B_resp_two.P);
        ec_dbl_iter(&B_resp_two.Q, pow_dim2_deg_resp + 2, &isog.codomain.E2, &B_resp_two.Q);
        ec_dbl_iter(&B_resp_two.PmQ, pow_dim2_deg_resp + 2, &isog.codomain.E2, &B_resp_two.PmQ);

        assert(test_point_order_twof(&B_resp_two.P, &isog.codomain.E2, exp_diadic_val_full_resp));
        assert(test_point_order_twof(&B_resp_two.Q, &isog.codomain.E2, exp_diadic_val_full_resp));
        assert(test_point_order_twof(&B_resp_two.PmQ, &isog.codomain.E2, exp_diadic_val_full_resp));

        ec_point_t ker;
        // applyling the vector to find the kernel
        ec_biscalar_mul_ibz(&ker,
                            &isog.codomain.E2,
                            &vec_resp_two[0],
                            &vec_resp_two[1],
                            &B_resp_two,
                            exp_diadic_val_full_resp);
#ifndef NDEBUG
        if (ibz_cmp(&vec_resp_two[0], &ibz_const_zero) == 0 &&
            ibz_cmp(&vec_resp_two[1], &ibz_const_one) == 0) {
            assert(ec_is_equal(&ker, &B_resp_two.Q));
        }
#endif
        assert(test_point_order_twof(&ker, &isog.codomain.E2, exp_diadic_val_full_resp));

        // computing the isogeny and pushing the points
        ec_eval_small_chain(&E_chall_2, &ker, exp_diadic_val_full_resp, points, 3);

        // copying the result
        copy_point(&B_resp_two.P, &points[0]);
        copy_point(&B_resp_two.Q, &points[1]);
        copy_point(&B_resp_two.PmQ, &points[2]);

#ifndef NDEBUG
        fp2_t w0;
        ec_point_t AC, A24;
        fp2_copy(&AC.x, &E_chall_2.A);
        fp2_copy(&AC.z, &E_chall_2.C);
        A24_from_AC(&A24, &AC);
        weil(&w0,
             pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp,
             &B_resp_two.P,
             &B_resp_two.Q,
             &B_resp_two.PmQ,
             &A24);
#endif
    }

    // computing the codomain of the challenge
    // we need to settle the exact kernel point we are going to use

    if (ibz_get(&vec_chall[0]) % 2 == 1) {
        sig->chall_b = 0;
        ibz_pow(&tmp, &ibz_const_two, TORSION_PLUS_EVEN_POWER);
        ibz_invmod(&sig->chall_coeff, &vec_chall[0], &tmp);
        ibz_mul(&sig->chall_coeff, &vec_chall[1], &sig->chall_coeff);
        ibz_set(&vec_chall[0], 1);
    } else {
        sig->chall_b = 1;
        ibz_pow(&tmp, &ibz_const_two, TORSION_PLUS_EVEN_POWER);
        ibz_invmod(&sig->chall_coeff, &vec_chall[1], &tmp);
        ibz_mul(&sig->chall_coeff, &vec_chall[0], &sig->chall_coeff);
        ibz_set(&vec_chall[1], 1);
    }

    ec_isog_even_t phi_chall;
    ec_basis_t bas_sk;
    copy_point(&bas_sk.P, &sk->canonical_basis.P);
    copy_point(&bas_sk.Q, &sk->canonical_basis.Q);
    copy_point(&bas_sk.PmQ, &sk->canonical_basis.PmQ);
    phi_chall.curve = sk->curve;
    phi_chall.length = TORSION_PLUS_EVEN_POWER - backtracking;
    assert(test_point_order_twof(&bas_sk.P, &sk->curve, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&bas_sk.Q, &sk->curve, TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&bas_sk.PmQ, &sk->curve, TORSION_PLUS_EVEN_POWER));
    ec_biscalar_mul_ibz(&phi_chall.kernel,
                        &sk->curve,
                        &vec_chall[0],
                        &vec_chall[1],
                        &bas_sk,
                        TORSION_PLUS_EVEN_POWER);
    assert(test_point_order_twof(&phi_chall.kernel, &sk->curve, TORSION_PLUS_EVEN_POWER));
    for (int i = 0; i < backtracking; i++) {
        ec_dbl(&phi_chall.kernel, &sk->curve, &phi_chall.kernel);
    }

    ec_curve_t Echall = sk->curve;
    assert(test_point_order_twof(&phi_chall.kernel, &Echall, phi_chall.length));
    ec_eval_even(&Echall, &phi_chall, &bas_sk.P, 1);

#ifndef NDEBUG
    fp2_t j_chall, j_test1, j_test2;
    ec_j_inv(&j_test1, &isog.codomain.E1);
    ec_j_inv(&j_test2, &E_chall_2);
    ec_j_inv(&j_chall, &Echall);
    // apparently its always the second one
    assert(fp2_is_equal(&j_chall, &j_test2) || fp2_is_equal(&j_chall, &j_test1));

#endif

    // applying the isomorphism from E_chall_2 to Echall
    ec_isom_t isom;
    ec_isomorphism(&isom, &E_chall_2, &Echall);
    ec_iso_eval(&B_resp_two.P, &isom);
    ec_iso_eval(&B_resp_two.Q, &isom);
    ec_iso_eval(&B_resp_two.PmQ, &isom);

#ifndef NDEBUG
    fp2_t w0;
    ec_point_t AC, A24;
    fp2_copy(&AC.x, &Echall.A);
    fp2_copy(&AC.z, &Echall.C);
    A24_from_AC(&A24, &AC);
    weil(&w0,
         pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp,
         &B_resp_two.P,
         &B_resp_two.Q,
         &B_resp_two.PmQ,
         &A24);
#endif

    // now it only remains to format the response for the verification

    // canonical basis
    ec_basis_t B_can_chall, B_aux2, B_aux2_can;
    ec_curve_t E_aux2;
    // it should always be the first curve
    copy_curve(&E_aux2, &isog.codomain.E1);
    copy_point(&B_aux2.P, &Tev1.P1);
    copy_point(&B_aux2.Q, &Tev2.P1);
    copy_point(&B_aux2.PmQ, &Tev1m2.P1);
    ec_curve_to_basis_2f_to_hint(
        &B_can_chall, &Echall, pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp, sig->hint_chall);
    ec_curve_to_basis_2f_to_hint(
        &B_aux2_can, &E_aux2, pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp, sig->hint_aux);

    assert(test_point_order_twof(
        &B_aux2.P, &E_aux2, pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp));
    assert(test_point_order_twof(
        &B_aux2.Q, &E_aux2, pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp));
    assert(test_point_order_twof(
        &B_aux2.PmQ, &E_aux2, pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp));
#ifndef NDEBUG
    fp2_copy(&AC.x, &E_aux2.A);
    fp2_copy(&AC.z, &E_aux2.C);
    A24_from_AC(&A24, &AC);
    weil(&w0,
         pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp,
         &B_aux2.P,
         &B_aux2.Q,
         &B_aux2.PmQ,
         &A24);
    weil(&w0,
         pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp,
         &B_aux2_can.P,
         &B_aux2_can.Q,
         &B_aux2_can.PmQ,
         &A24);
#endif

    // compute the matrix to go from B_aux2 to B_aux2_can
    change_of_basis_matrix_two(&mat_Baux2_to_Baux2_can,
                               &B_aux2_can,
                               &B_aux2,
                               &E_aux2,
                               pow_dim2_deg_resp + 2 + exp_diadic_val_full_resp);

    // apply the change of basis to B_resp_two
    matrix_application_even_basis(&B_resp_two,
                                  &Echall,
                                  &mat_Baux2_to_Baux2_can,
                                  pow_dim2_deg_resp + exp_diadic_val_full_resp + 2);

    // compute the matrix to go from B_chall_can to B_resp_two
    change_of_basis_matrix_two(&mat_Bchall_can_to_Bchall,
                               &B_resp_two,
                               &B_can_chall,
                               &Echall,
                               pow_dim2_deg_resp + exp_diadic_val_full_resp + 2);


    // filling the output
    sig->backtracking = backtracking;
    sig->two_resp_length = exp_diadic_val_full_resp;
    ibz_mat_2x2_copy(&sig->mat_Bchall_can_to_B_chall, &mat_Bchall_can_to_Bchall);
    // setting sig->E_aux
    fp2_t temp_fp2;
    fp2_copy(&temp_fp2, &E_aux2.C);
    fp2_inv(&temp_fp2);
    fp2_mul(&sig->E_aux.A, &temp_fp2, &E_aux2.A);
    fp2_set_one(&sig->E_aux.C);
    ec_point_init(&sig->E_aux.A24);
    sig->E_aux.is_A24_computed_and_normalized = 0;

    ibz_vec_2_finalize(&vec);
    ibz_vec_2_finalize(&vec_chall);
    ibz_vec_2_finalize(&vec_resp_two);
    ibz_mat_2x2_finalize(&mat_Bchall_can_to_Bchall);
    ibz_mat_2x2_finalize(&mat_Baux2_to_Baux2_can);

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
    quat_left_ideal_init(&lideal_aux_resp_com);
    quat_left_ideal_finalize(&lideal_aux);

    ibz_vec_4_finalize(&dummy_coord);
    ibz_finalize(&degree_full_resp);
    ibz_finalize(&degree_odd_resp);

    ibq_finalize(&temp_norm);
    ibz_finalize(&tmp);
    ibz_finalize(&lattice_content);
    ibz_finalize(&remain);
    return 1;
}

int
protocols_verif(signature_t *sig, const public_key_t *pk, const unsigned char *m, size_t l)
{

    int verif;

    ibz_t tmp;
    ibz_vec_2_t vec_chall, check_vec_chall;

    ibz_init(&tmp);
    ibz_vec_2_init(&vec_chall);
    ibz_vec_2_init(&check_vec_chall);

    // checking that we are given A coefficients and no precomputation
    assert(fp2_is_one(&pk->curve.C) && !pk->curve.is_A24_computed_and_normalized);
    assert(fp2_is_one(&sig->E_aux.C) && !sig->E_aux.is_A24_computed_and_normalized);

    clock_t t = tic();

    // computation of the challenge
    ec_isog_even_t phi_chall;
    ec_basis_t bas_EA;
    ec_curve_t Epk;
    copy_curve(&Epk, &pk->curve);
    // ec_curve_normalize_A24(&Epk);
    ec_curve_to_basis_2f_from_hint(
        &bas_EA, &Epk, TORSION_PLUS_EVEN_POWER, pk->hint_pk); // canonical
    phi_chall.curve = Epk;
    phi_chall.length = TORSION_PLUS_EVEN_POWER - sig->backtracking;

    // recovering the exact vec_chall
    if (sig->chall_b) {
        ibz_copy(&vec_chall[0], &sig->chall_coeff);
        ibz_set(&vec_chall[1], 1);
    } else {
        ibz_copy(&vec_chall[1], &sig->chall_coeff);
        ibz_set(&vec_chall[0], 1);
    }
    digit_t scal[NWORDS_ORDER] = { 0 };
    ibz_to_digit_array(scal, &sig->chall_coeff);
    if (sig->chall_b) {
        ec_ladder3pt(&phi_chall.kernel, scal, &bas_EA.Q, &bas_EA.P, &bas_EA.PmQ, &Epk);
    } else {
        ec_ladder3pt(&phi_chall.kernel, scal, &bas_EA.P, &bas_EA.Q, &bas_EA.PmQ, &Epk);
    }

    ec_dbl_iter(&phi_chall.kernel, sig->backtracking, &Epk, &phi_chall.kernel);

    ec_curve_t Echall = Epk;
    ec_eval_even(&Echall, &phi_chall, &bas_EA.P, 1);
    // printf("challenge computation length : %d ",phi_chall.length);
    // TOC_clock(t,"");

    int pow_dim2_deg_resp = SQIsign2D_response_length - sig->two_resp_length;

    ec_basis_t B_chall_can, B_aux_can;
    ec_curve_t E_aux;
    copy_curve(&E_aux, &sig->E_aux);

    // recovering the canonical basis
    ec_curve_to_basis_2f_from_hint(
        &B_chall_can, &Echall, pow_dim2_deg_resp + 2 + sig->two_resp_length, sig->hint_chall);
    ec_curve_to_basis_2f_from_hint(
        &B_aux_can, &E_aux, pow_dim2_deg_resp + 2 + sig->two_resp_length, sig->hint_aux);

    // TOC_clock(t,"challenge and canonical basis");


    // setting to the right order
    ec_dbl_iter(&B_aux_can.P, sig->two_resp_length, &E_aux, &B_aux_can.P);
    ec_dbl_iter(&B_aux_can.Q, sig->two_resp_length, &E_aux, &B_aux_can.Q);
    ec_dbl_iter(&B_aux_can.PmQ, sig->two_resp_length, &E_aux, &B_aux_can.PmQ);

    assert(test_point_order_twof(
        &B_chall_can.P, &Echall, 2 + pow_dim2_deg_resp + sig->two_resp_length));
    assert(test_point_order_twof(
        &B_chall_can.Q, &Echall, 2 + pow_dim2_deg_resp + sig->two_resp_length));
    assert(test_point_order_twof(
        &B_chall_can.PmQ, &Echall, 2 + pow_dim2_deg_resp + sig->two_resp_length));

    // applying the change matrix on the basis of E_chall
    matrix_application_even_basis(&B_chall_can,
                                  &Echall,
                                  &sig->mat_Bchall_can_to_B_chall,
                                  pow_dim2_deg_resp + 2 + sig->two_resp_length);

    // TOC_clock(t,"matrix application");


    ec_curve_t E_chall_2;
    ec_curve_t mem;
    copy_curve(&E_chall_2, &Echall);
    copy_curve(&mem, &Echall);

    if (sig->two_resp_length > 0) {
        // computing the small two chain
        ec_curve_t E_chall_2;
        copy_curve(&E_chall_2, &Echall);
        ec_point_t ker, points[3];

        // choosing the right point for the small two_isogenies
        if (ibz_get(&sig->mat_Bchall_can_to_B_chall[0][0]) % 2 == 0 &&
            ibz_get(&sig->mat_Bchall_can_to_B_chall[1][0]) % 2 == 0) {
            copy_point(&ker, &B_chall_can.Q);
        } else {
            copy_point(&ker, &B_chall_can.P);
        }

        copy_point(&points[0], &B_chall_can.P);
        copy_point(&points[1], &B_chall_can.Q);
        copy_point(&points[2], &B_chall_can.PmQ);
        ec_dbl_iter(&ker, pow_dim2_deg_resp + 2, &Echall, &ker);
        // we extract the kernel by multiplying by 2^(2+pow_dim2_deg_resp)
        // for (int i=0;i<pow_dim2_deg_resp+2;i++) {
        //     ec_dbl(&ker,&Echall,&ker);
        // }
        assert(test_point_order_twof(&ker, &E_chall_2, sig->two_resp_length));
        ec_eval_small_chain(&E_chall_2, &ker, sig->two_resp_length, points, 3);

        assert(test_point_order_twof(&points[0], &E_chall_2, 2 + pow_dim2_deg_resp));
        assert(test_point_order_twof(&points[1], &E_chall_2, 2 + pow_dim2_deg_resp));
        assert(test_point_order_twof(&points[2], &E_chall_2, 2 + pow_dim2_deg_resp));
        copy_point(&B_chall_can.P, &points[0]);
        copy_point(&B_chall_can.Q, &points[1]);
        copy_point(&B_chall_can.PmQ, &points[2]);
        copy_curve(&mem, &E_chall_2);
    }
    // due to a weird error that I couldn't understand
    copy_curve(&E_chall_2, &mem);
    assert(test_point_order_twof(&B_chall_can.P, &E_chall_2, 2 + pow_dim2_deg_resp));
    assert(test_point_order_twof(&B_chall_can.Q, &E_chall_2, 2 + pow_dim2_deg_resp));
    assert(test_point_order_twof(&B_chall_can.PmQ, &E_chall_2, 2 + pow_dim2_deg_resp));
    assert(test_point_order_twof(&B_aux_can.P, &E_aux, 2 + pow_dim2_deg_resp));
    assert(test_point_order_twof(&B_aux_can.Q, &E_aux, 2 + pow_dim2_deg_resp));
    assert(test_point_order_twof(&B_aux_can.PmQ, &E_aux, 2 + pow_dim2_deg_resp));

    // now compute the dim2 isogeny from E_chall_2 x E_aux -> E_com x E_aux'
    // of kernel B_chall_can x B_aux_can

    // first we set-up the kernel
    theta_couple_curve_t EchallxEaux;
    theta_couple_point_t T1, T2, T1m2;
    theta_chain_t isog;
    copy_curve(&EchallxEaux.E1, &E_chall_2);
    copy_curve(&EchallxEaux.E2, &E_aux);
    copy_point(&T1.P2, &B_aux_can.P);
    copy_point(&T2.P2, &B_aux_can.Q);
    copy_point(&T1m2.P2, &B_aux_can.PmQ);
    copy_point(&T1.P1, &B_chall_can.P);
    copy_point(&T2.P1, &B_chall_can.Q);
    copy_point(&T1m2.P1, &B_chall_can.PmQ);

    // computing the isogeny
    int extra_info = 1;
    theta_chain_comput_strategy_faster_no_eval(
        &isog,
        pow_dim2_deg_resp,
        &EchallxEaux,
        &T1,
        &T2,
        &T1m2,
        strategies[TORSION_PLUS_EVEN_POWER - pow_dim2_deg_resp],
        extra_info);

    // TOC_clock(t,"response isogeny");

    // computing the commitment curve
    // apparently its always the second one
    ec_curve_t E_com;
    copy_curve(&E_com, &isog.codomain.E2);

    // recomputing the challenge vector
    hash_to_challenge(&check_vec_chall, &E_com, m, pk, l);

    // performing the final check
    if (sig->chall_b) {
        ibz_mul(&vec_chall[0], &vec_chall[0], &check_vec_chall[1]);
        verif = (ibz_cmp(&vec_chall[0], &check_vec_chall[0]) == 0);
    } else {
        ibz_mul(&vec_chall[1], &vec_chall[1], &check_vec_chall[0]);
        verif = (ibz_cmp(&vec_chall[1], &check_vec_chall[1]) == 0);
    }

    ibz_finalize(&tmp);
    ibz_vec_2_finalize(&vec_chall);
    ibz_vec_2_finalize(&check_vec_chall);

    return verif;
}
