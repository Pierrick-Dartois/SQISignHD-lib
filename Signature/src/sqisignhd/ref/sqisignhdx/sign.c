#include <sqisignhd.h>
#include <curve_extras.h>
#include <tools.h>
#include <fips202.h>
#include <stdio.h>
#include <string.h>
#include <fp2.h>

#define RESPONSE_LENGTH TORSION_PLUS_EVEN_POWER + 16

const clock_t time_isogenies_odd = 0;
const clock_t time_sample_response = 0;
const clock_t time_change_of_basis_matrix = 0;

void
secret_sig_init(signature_t *sig)
{
    sig->hint_com = (int *)malloc(2 * sizeof(int));
    //sig->hint_chal = (int *)malloc(2 * sizeof(int));

    ibz_init(&sig->a);
    ibz_init(&sig->b);
    ibz_init(&sig->c_or_d);
    ibz_init(&sig->q);
    ibz_init(&sig->chal);
    //ibz_init(&sig->x);
    //ibz_init(&sig->b0);
    //ibz_init(&sig->d0);
    //ibz_init(&sig->b1);
    //ibz_init(&sig->d1);
    //ibz_init(&sig->e0_adjust);
    //ibz_init(&sig->c0_adjust);
}

void
secret_sig_finalize(signature_t *sig)
{
    free(sig->hint_com);
    //free(sig->hint_chal);

    ibz_finalize(&sig->a);
    ibz_finalize(&sig->b);
    ibz_finalize(&sig->c_or_d);
    ibz_finalize(&sig->q);
    ibz_finalize(&sig->chal);
    //ibz_finalize(&sig->x);
    //ibz_finalize(&sig->b0);
    //ibz_finalize(&sig->d0);
    //ibz_finalize(&sig->b1);
    //ibz_finalize(&sig->d1);
    //ibz_finalize(&sig->e0_adjust);
    //ibz_finalize(&sig->c0_adjust);
}

void
fprint_signature(FILE *p_file, const signature_t *sig)
{
    fprintf(p_file, "A_com = ");
    fp2_print_to_file(p_file,&(sig->E_com.A));
    ibz_fprintf(p_file,"a = %Zx\n",sig->a);
    ibz_fprintf(p_file,"b = %Zx\n",sig->b);
    ibz_fprintf(p_file,"c_or_d = %Zx\n",sig->c_or_d);
    ibz_fprintf(p_file,"q = %Zx\n",sig->q);
    fprintf(p_file, "hint_com_P = %u\n",(sig->hint_com[0]));
    fprintf(p_file, "hint_com_Q = %u\n",(sig->hint_com[1]));
    //fprintf(p_file, "hint_chal_P = %u\n",(sig->hint_chal[0]));
    //fprintf(p_file, "hint_chal_Q = %u\n",(sig->hint_chal[1]));
    ibz_fprintf(p_file, "chal = %Zx\n",(sig->chal));
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
commit(ec_curve_t *E_com, quat_left_ideal_t *lideal_com, ec_basis_t *B_com)
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

    // pushing the points of the basis
    // first we create the points
    theta_couple_point_t EvP, EvQ, EvPmQ;
    theta_couple_point_t ResP, ResQ, ResPmQ;
    copy_point(&EvP.P1, &BASIS_EVEN.P);
    copy_point(&EvQ.P1, &BASIS_EVEN.Q);
    copy_point(&EvPmQ.P1, &BASIS_EVEN.PmQ);
    ec_set_zero(&EvP.P2);
    ec_set_zero(&EvQ.P2);
    ec_set_zero(&EvPmQ.P2);
    theta_chain_eval_special_case(&ResP, &F, &EvP, &F.domain);
    theta_chain_eval_special_case(&ResQ, &F, &EvQ, &F.domain);
    theta_chain_eval_special_case(&ResPmQ, &F, &EvPmQ, &F.domain);

    assert(found);

    copy_point(&B_com->P, &ResP.P2);
    copy_point(&B_com->Q, &ResQ.P2);
    copy_point(&B_com->PmQ, &ResPmQ.P2);

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
is_good_norm(ibz_t *N)
{
    if ((4 - (ibz_get(N) % 4)) != 1)
        return 0;

    ibz_t sum_of_squares_candidate;
    //ibz_init(&pow2);
    ibz_init(&sum_of_squares_candidate);
    int res = 0;

    //ibz_set(&pow2, 1);
    //ibz_mul_2exp(&pow2, &pow2, ibz_bitsize(N));

    if (ibz_cmp(&DEGREE_RESP_HD, N) < 0) {
        ibz_printf(
            "WARNING: short vectors not short enough...\n2-pow = %Zd\nnorm = %Zd\n", &DEGREE_RESP_HD, N);
        // assert(0);
        ibz_finalize(&sum_of_squares_candidate);
        return 0;
    }

    ibz_sub(&sum_of_squares_candidate, &DEGREE_RESP_HD, N);

    // unsigned int N_mod_four = ibz_mod_ui (&N, 4);
    assert(ibz_mod_ui(&sum_of_squares_candidate, 4) == 1);

    // if (N_mod_four == 1) {
    res = ibz_probab_prime(&sum_of_squares_candidate, 40);

    //ibz_finalize(&pow2);
    ibz_finalize(&sum_of_squares_candidate);
    return res;
}

int
sample_response(quat_alg_elem_t *x,
                const quat_lattice_t *lattice,
                ibz_t const *lattice_content,
                int verbose)
{
    ibz_mat_4x4_t lll;
    ibz_t denom_gram, norm;

    ibz_mat_4x4_init(&lll);
    ibz_init(&denom_gram);
    ibz_init(&norm);

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
    int m = 10;
    while (!found && cnt < 2 * (2 * m + 1) * (2 * m + 1) * (2 * m + 1) * (2 * m + 1)) {
        cnt++;
        for (int i = 0; i < 4; i++) {
            ibz_rand_interval_minm_m(&vec[i], m);
        }
        norm_from_2_times_gram(&norm, &gram, &vec);
        // now we test if the norm is good
        // there are two constraints : the norm n must be such that 2^e - n is a sum of two square
        // where e is the smallest exponent such 2^e > n
        // TODO we could use other values of e (since we can split the chain in two in dim 4, the
        // bound on e is much higher) and the element must be primitive in O0 (this ensures that
        // there is no backtracking -> in SQIsignHD???)
        found = is_good_norm(&norm); //&& Not sure this is useful to require non-zero coordinates 
                //(ibz_cmp(&vec[0], &ibz_const_zero) != 0 || ibz_cmp(&vec[1], &ibz_const_zero) != 0 ||
                 //ibz_cmp(&vec[2], &ibz_const_zero) != 0 || ibz_cmp(&vec[3], &ibz_const_zero) != 0);
        if (found) {
            ibz_mat_4x4_eval(&(x->coord), &lll, &vec);
            assert(quat_lattice_contains(NULL, lattice, x, &QUATALG_PINFTY));
            //if (!quat_alg_is_primitive(x, &MAXORD_O0, &QUATALG_PINFTY)) {
                //found = 0; // Why should this be primitive?
            //}
        }
    }
    assert(quat_lattice_contains(NULL, lattice, x, &QUATALG_PINFTY));

    ibz_finalize(&denom_gram);
    ibz_finalize(&norm);
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&lll);
    ibz_vec_4_finalize(&vec);
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
    //for(int i=0;i<FP2_ENCODED_BYTES + FP2_ENCODED_BYTES + length;i++){
        //printf("%c",buf[i]);
    //}
    //printf("\n");

    // TODO(security) omit some vectors, notably (a,1) with gcd(a,6)!=1 but also things like (2,3)?
    {
        digit_t digits[NWORDS_FIELD];

        // FIXME should use SHAKE128 for smaller parameter sets?
        //  TODO we want to use a bit differently (first we hash first half and then derive the
        //  second half)
        SHAKE256(
            (void *)digits, sizeof(digits), buf, FP2_ENCODED_BYTES + FP2_ENCODED_BYTES + length);
        //for (int i = 0; i < SQIsign2D_heuristic_challenge_hash_iteration; i++) { // Grinding is unnecessary here -> Challenge length is enough
            //SHAKE256((void *)digits, sizeof(digits), (void *)digits, sizeof(digits));
        //}

        ibz_set(&(*scalars)[1], 1); // FIXME
        ibz_copy_digit_array(&(*scalars)[1], digits);
        //ibz_printf("%Zx\n",(*scalars)[1]);
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

    ibz_t lattice_content; //pow_chall;
    ec_curve_t E_com;
    ec_basis_t B_com0; // basis of 2^TORSION_PLUS_EVEN_POWER
    ec_basis_t B_resp_two;
    ibz_vec_4_t coeffs;
    ibz_vec_2_t vec, vec_chall, vec_resp_two;
    quat_left_ideal_t lideal_tmp;
    quat_left_ideal_t lideal_commit, lideal_chall_two;
    quat_left_ideal_t lideal_chall_secret, lideal_com_resp, lideal_resp_two;
    quat_lattice_t lattice_hom_chall_to_com, lat_commit;
    quat_alg_elem_t resp_quat;
    quat_alg_elem_t elem_tmp;
    ibz_mat_2x2_t mat_Bcom0_to_Bcom_can, mat_Bchall_can_to_Bchall, mat, sig_mat_pk_can_to_B_pk, mat_B_chal_can_to_B_chal_pk;
    // ibz_mat_2x2_t mat_sigma_phichall_BA_to_Bcomcan, mat_sigma_phichall_BA0_to_Bcom0;
    ibz_t degree_com_isogeny, tmp, remain, chal;
    ibz_t degree_full_resp; //degree_odd_resp;
    ibq_t temp_norm;
    int exp_diadic_val_full_resp;
    int pow_dim2_deg_resp;
    int backtracking;

    int found = 1;

    ec_curve_init(&E_com);
    ibz_init(&tmp);
    ibz_init(&lattice_content);
    ibz_init(&remain);
    ibz_init(&chal);

    ibz_init(&degree_full_resp);
    //ibz_init(&degree_odd_resp);
    ibq_init(&temp_norm);
    //ibz_init(&pow_chall);

    ibz_mat_2x2_init(&mat_Bchall_can_to_Bchall);
    ibz_mat_2x2_init(&mat_Bcom0_to_Bcom_can);
    ibz_mat_2x2_init(&sig_mat_pk_can_to_B_pk);
    ibz_mat_2x2_init(&mat);
    ibz_mat_2x2_init(&mat_B_chal_can_to_B_chal_pk);
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

    ibz_vec_2_init(&vec);
    ibz_vec_2_init(&vec_chall);
    ibz_vec_2_init(&vec_resp_two);

    // * Computing the commitment
    commit(&E_com, &lideal_commit, &B_com0);

    // * Computing the challenge
    // challenge length
    ///int len_chall = SQIsign2D_heuristic_challenge_length;
    ///ibz_pow(&pow_chall, &ibz_const_two, len_chall);

    // computing the challenge
    // vec_chall is a pair of coefficients encoding the kernel of the challenge isogeny
    // as vec_chall[0]*B[0] + vec_chall[1]*B[1] where B is the canonical basis of the 2^EXPONENT_CHAL_HD
    // torsion of EA
    hash_to_challenge(&vec_chall, &E_com, m, pk, l);

    // computing the challenge isogeny 
    /*
    ec_isog_even_t phi_chall;
    ec_basis_t B_pk_can;
    ec_curve_t Epk;
    copy_curve(&Epk, &pk->curve);
    ec_curve_to_basis_2f_from_hint(
        &B_pk_can, &Epk, TORSION_PLUS_EVEN_POWER, pk->hint_pk); // canonical basis of Epk
    phi_chall.curve = Epk;
    phi_chall.length = EXPONENT_CHAL_HD;

    int chal_b = ibz_get(&vec_chall[0])%2;
    if(chal_b){
        // K_chal = P_pk + [chal]*Q_pk
        ibz_invmod(&tmp,&vec_chall[0],&DEGREE_CHAL_HD);
        ibz_mul(&tmp,&tmp,&vec_chall[1]);
        ibz_mod(&chal,&tmp,&DEGREE_CHAL_HD);
    } 
    else{
        // K_chal = [chal]*P_pk + Q_pk
        ibz_invmod(&tmp,&vec_chall[1],&DEGREE_CHAL_HD);
        ibz_mul(&tmp,&tmp,&vec_chall[0]);
        ibz_mod(&chal,&tmp,&DEGREE_CHAL_HD);
    }

    digit_t scal[NWORDS_ORDER] = { 0 };
    ibz_to_digit_array(scal, &chal);
    if (chal_b) {
        // K_chal = P_pk + [chal]*Q_pk
        ec_ladder3pt(&phi_chall.kernel, scal, &B_pk_can.P, &B_pk_can.Q, &B_pk_can.PmQ, &Epk);
    } else {
        // K_chal = [chal]*P_pk + Q_pk
        ec_ladder3pt(&phi_chall.kernel, scal, &B_pk_can.Q, &B_pk_can.P, &B_pk_can.PmQ, &Epk);
    }

    ec_dbl_iter(&phi_chall.kernel, TORSION_PLUS_EVEN_POWER-EXPONENT_CHAL_HD, &Epk, &phi_chall.kernel);
    */

    //ec_curve_t Echall;
    //ec_point_t points[3]; 
    //copy_point(&points[0],&B_pk_can.P);
    //copy_point(&points[1],&B_pk_can.Q);
    //copy_point(&points[2],&B_pk_can.PmQ);
    //ec_eval_even(&Echall, &phi_chall, points, 1);

    // Computing the matrix of the image B_chal_pk = phi_chal(B_pk_can) in B_chal_can
    //ec_basis_t B_chal_pk;
    //copy_point(&B_chal_pk.P,&points[0]);
    //copy_point(&B_chal_pk.Q,&points[1]);
    //copy_point(&B_chal_pk.PmQ,&points[2]);

    /*
    // Computation of mat_B_chal_can_to_B_chal_pk, the matrix of B_chal_pk=phi_chal(B_pk_can) in B_chal_can
    // mat_B_chal_can_to_B_chal_pk=[[a,c],[b,d]], phi_chal(P_pk)=[a]P_chal+[b]Q_chal, phi_chal(Q_pk)=[c]P_chal+[d]Q_chal.
    ec_basis_t B_chal_can, phi_chall_B_pk;
    ec_curve_t Echall;
    ec_point_t im_phi_chall[3];

    copy_point(&im_phi_chall[0],&B_pk_can.P);
    copy_point(&im_phi_chall[1],&B_pk_can.Q);
    copy_point(&im_phi_chall[2],&B_pk_can.PmQ);

    ec_eval_even(&Echall, &phi_chall, im_phi_chall, 3);

    copy_point(&phi_chall_B_pk.P,&im_phi_chall[0]);
    copy_point(&phi_chall_B_pk.Q,&im_phi_chall[1]);
    copy_point(&phi_chall_B_pk.PmQ,&im_phi_chall[2]);

    ec_curve_to_basis_2f_to_hint(&B_chal_can, &Echall, TORSION_PLUS_EVEN_POWER, sig->hint_chal);

    change_of_basis_matrix_two_robust(&mat_B_chal_can_to_B_chal_pk, &phi_chall_B_pk, &B_chal_can, &Echall, TORSION_PLUS_EVEN_POWER);

    //printf("After challenge basis push\n");*/


    // now we compute the ideal associated to the challenge
    // for that, we need to find vec such that
    // the kernel of the challenge isogeny is generated by vec[0]*B0[0] + vec[1]*B0[1] where B0 is
    // the image through the secret key isogeny of the canonical basis E0
    ibz_mat_2x2_eval(&vec, &(sk->mat_BAcan_to_BA0_two), &vec_chall);

    // reducing mod 2^EXPONENT_CHAL_HD
    ibz_mod(&vec[0], &vec[0], &DEGREE_CHAL_HD);
    ibz_mod(&vec[1], &vec[1], &DEGREE_CHAL_HD);

    // lideal_chall_two is the pullback of the ideal challenge through the secret key ideal
    int e_chal = EXPONENT_CHAL_HD;
    id2iso_kernel_dlogs_to_ideal_two(&lideal_chall_two, &vec, e_chal);
    assert(ibz_cmp(&lideal_chall_two.norm, &DEGREE_CHAL_HD) == 0);

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
        return 0;
    }

    assert(quat_lattice_contains(NULL, &MAXORD_O0, &resp_quat, &QUATALG_PINFTY));
    //assert(quat_alg_is_primitive(&resp_quat, &MAXORD_O0, &QUATALG_PINFTY)); // No need to be primitive

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
    assert(exp_diadic_val_full_resp == 0); // There is no 2-adic valuation
    /// REMOVE
    /// removing the power of two part
    ///ibz_pow(&tmp, &ibz_const_two, exp_diadic_val_full_resp);
    ///ibz_div(&degree_odd_resp, &remain, &degree_full_resp, &tmp);
    ///assert(ibz_cmp(&remain, &ibz_const_zero) == 0);

    // testing the bound of the dim 4 response q<2^f
    ///pow_dim2_deg_resp = ibz_bitsize(&degree_odd_resp);
    ///ibz_pow(&remain, &ibz_const_two, pow_dim2_deg_resp);
    ibz_sub(&tmp, &DEGREE_RESP_HD, &degree_full_resp);
    assert(ibz_cmp(&tmp, &ibz_const_zero) > 0);

    // now it only remains to format the response for the verification
    /*** Goal: evaluate B_pk:=1/(N_sk N_com) phi_chl*phi_sk*bar(quat_resp)*dual(phi_com)(B_com_can)=[2^e_chal]phi_resp(B_com_can) 
     * The evaluation by phi_chl is done during the verification 
     * We express the matrix from B_pk_can to B_pk: 
     * sig_mat_pk_can_to_B_pk = 1/(N_sk N_com) mat_B_chal_can_to_B_chal_pk * mat_BAcan_to_BA0 * mat_B0_to_quat_resp_dual_0 * mat_B0_to_B_com_dual ***/

    // * Commitment: mat_B0_to_B_com_dual = N_com * mat_Bcom0_to_Bcom_can
    // canonical basis
    ec_basis_t B_com_can;
    ec_curve_to_basis_2f_to_hint(&B_com_can, &E_com, TORSION_PLUS_EVEN_POWER, sig->hint_com);

    // compute the matrix to go from B_com0 to B_com_can
    change_of_basis_matrix_two(
        &mat_Bcom0_to_Bcom_can, &B_com_can, &B_com0, &E_com, TORSION_PLUS_EVEN_POWER);

    // * quat_resp: mat = mat_B0_to_quat_resp_dual_0
    // taking the conjugate of quat_resp so that quat_resp is contained in lideal_com
    quat_alg_conj(&resp_quat, &resp_quat);

    /// this will be the image of the basis of E0 through  phi_sec quat_resp / deg phi_sec
    /// where phi_sec is the isogeny corresponding to sk->lideal_secret
    /// ec_basis_t bas_sk;

    // we compute the matrix corresponding to resp_quat (now dual(resp_quat))
    quat_alg_make_primitive(&coeffs, &lattice_content, &resp_quat, &MAXORD_O0, &QUATALG_PINFTY); // `resp_quat` = `lattice_content` · Λ `coeffs`, where Λ is the basis of `order`
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
            ibz_mul(&mat[i][j], &mat[i][j], &lattice_content);// Multiplication by the non-primitive part
            ibz_mod(&mat[i][j], &mat[i][j], &TORSION_PLUS_2POWER);
        }
    }

    // * Matrix multiplications 
    // now we inverse the matrix of the secret key to get from the canonical basis of pk to phi(B0)
    ibz_2x2_inv_mod(&sig_mat_pk_can_to_B_pk, &sk->mat_BAcan_to_BA0_two, &TORSION_PLUS_2POWER);

    // Computation of mat_BAcan_to_BA0 * mat_B0_to_quat_resp_dual_0 (mat)
    ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk, &sig_mat_pk_can_to_B_pk, &mat, &TORSION_PLUS_2POWER);

    // Computation of mat_BAcan_to_BA0 * mat_B0_to_quat_resp_dual_0 * (1/N_com)*mat_B0_to_B_com_dual 
    ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk, &sig_mat_pk_can_to_B_pk, &mat_Bcom0_to_Bcom_can, &TORSION_PLUS_2POWER);

    /*// Computation of mat_B_chal_can_to_B_chal_pk * mat_BAcan_to_BA0 * mat_B0_to_quat_resp_dual_0 * (1/N_com)*mat_B0_to_B_com_dual 
    ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk, &mat_B_chal_can_to_B_chal_pk, &sig_mat_pk_can_to_B_pk, &TORSION_PLUS_2POWER);

    // The matrix should be a multiple of 2^e_chl (to be checked)
    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            ibz_div(&sig_mat_pk_can_to_B_pk[i][j], &remain, &sig_mat_pk_can_to_B_pk[i][j], &DEGREE_CHAL_HD);
            assert(ibz_cmp(&remain, &ibz_const_zero) == 0);
        }
    }
    */

    // Dividing by N_sk
    ibz_copy(&tmp, &sk->secret_ideal.norm);
    //assert(ibz_cmp(&tmp,&FIXED_DEGREE_SK) == 0);// Way to slow to impose this condition
    ibz_invmod(&tmp, &tmp, &TORSION_PLUS_2POWER);

    for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
            ibz_mul(&sig_mat_pk_can_to_B_pk[i][j], &sig_mat_pk_can_to_B_pk[i][j], &tmp);
            ibz_mod(&sig_mat_pk_can_to_B_pk[i][j], &sig_mat_pk_can_to_B_pk[i][j], &TORSION_PLUS_2POWER);
        }
    }

    // * Store into signature
    ibz_div(&tmp, &sig->a, &sig_mat_pk_can_to_B_pk[0][0], &SIGN_PT_ORDER_HD);
    //ibz_copy(&sig->a,&sig_mat_pk_can_to_B_pk[0][0]);

    ibz_mul(&tmp,&sig_mat_pk_can_to_B_pk[0][0],&vec_chall[1]);
    ibz_sub(&tmp,&sig_mat_pk_can_to_B_pk[1][0],&tmp);
    ibz_div(&sig->b, &remain, &tmp, &DEGREE_CHAL_HD);
    ibz_div(&tmp, &sig->b, &sig->b, &SIGN_PT_ORDER_HD);

    ibz_div(&tmp, &remain, &remain, &SIGN_PT_ORDER_HD);
    assert(ibz_cmp(&remain, &ibz_const_zero) == 0);
    
    //ibz_copy(&sig->b,&sig_mat_pk_can_to_B_pk[1][0]);

    if(ibz_get(&sig->a)%2 == 1){
        // a is odd, then c_or_d is c (determined by a)
        ibz_div(&tmp, &sig->c_or_d, &sig_mat_pk_can_to_B_pk[0][1], &SIGN_PT_ORDER_HD);
    }
    else{
        // a is even, then c_or_d is d (determined by b)
        ibz_mul(&tmp,&sig_mat_pk_can_to_B_pk[0][1],&vec_chall[1]);
        ibz_sub(&tmp,&sig_mat_pk_can_to_B_pk[1][1],&tmp);
        ibz_div(&sig->c_or_d, &remain, &tmp, &DEGREE_CHAL_HD);
        ibz_div(&tmp, &sig->c_or_d, &sig->c_or_d, &SIGN_PT_ORDER_HD);

        ibz_div(&tmp, &remain, &remain, &SIGN_PT_ORDER_HD);
        assert(ibz_cmp(&remain, &ibz_const_zero) == 0);
    }

    // Degree q of the response
    ibz_copy(&sig->q,&degree_full_resp);

    // Setting sig->E_com
    fp2_t temp_fp2;
    fp2_copy(&temp_fp2, &E_com.C);
    fp2_inv(&temp_fp2);
    fp2_mul(&sig->E_com.A, &temp_fp2, &E_com.A);
    fp2_set_one(&sig->E_com.C);
    ec_point_init(&sig->E_com.A24);
    sig->E_com.is_A24_computed_and_normalized = 0;

    // Setting challenge
    ibz_copy(&sig->chal,&vec_chall[1]);

    /*
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

    

    // apply the change of basis to the matrix
    // ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk,&mat_Bcom0_to_Bcom_can,&sig_mat_pk_can_to_B_pk,&TORSION_PLUS_2POWER);
    ibz_2x2_mul_mod(&sig_mat_pk_can_to_B_pk,
                    &sig_mat_pk_can_to_B_pk,
                    &mat_Bcom0_to_Bcom_can,
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
    // setting sig->E_com
    fp2_t temp_fp2;
    fp2_copy(&temp_fp2, &E_com.C);
    fp2_inv(&temp_fp2);
    fp2_mul(&sig->E_com.A, &temp_fp2, &E_com.A);
    fp2_set_one(&sig->E_com.C);
    ec_point_init(&sig->E_com.A24);
    sig->E_com.is_A24_computed_and_normalized = 0;

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
    */

    //ibz_finalize(&pow_chall);
    ibz_vec_2_finalize(&vec);
    ibz_vec_2_finalize(&vec_chall);
    ibz_vec_2_finalize(&vec_resp_two);
    ibz_mat_2x2_finalize(&mat_Bchall_can_to_Bchall);
    ibz_mat_2x2_finalize(&mat_Bcom0_to_Bcom_can);
    ibz_mat_2x2_finalize(&sig_mat_pk_can_to_B_pk);
    ibz_mat_2x2_finalize(&mat);
    ibz_mat_2x2_finalize(&mat_B_chal_can_to_B_chal_pk);

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

    ibz_vec_4_finalize(&coeffs);
    ibz_finalize(&degree_full_resp);
    //ibz_finalize(&degree_odd_resp);

    ibq_finalize(&temp_norm);
    ibz_finalize(&tmp);
    ibz_finalize(&lattice_content);
    ibz_finalize(&remain);
    ibz_finalize(&chal);

    return found;
}
