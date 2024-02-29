#include <sqisignhd.h>
#include <curve_extras.h>
#include <toolbox.h>
#include <fips202.h>
#include <stdio.h>
#include <string.h>

#define RESPONSE_LENGTH TORSION_PLUS_EVEN_POWER+16

const clock_t time_isogenies_odd = 0;
const clock_t time_sample_response = 0;
const clock_t time_change_of_basis_matrix = 0;

void secret_sig_init(signature_t *sig) {
    ibz_mat_2x2_init(&(sig->mat_sigma_phichall));
}

void secret_sig_finalize(signature_t *sig) {
    ibz_mat_2x2_finalize(&(sig->mat_sigma_phichall));
}

static void ibz_vec_2_print2(char *name, const ibz_vec_2_t *vec){
    printf("%s", name);
    for(int i = 0; i < 2; i++){
        ibz_printf("%Zd ", &((*vec)[i]));
    }
    ibz_printf("\n");
}

static void ibz_vec_4_print2(char *name, const ibz_vec_4_t *vec){
    printf("%s", name);
    for(int i = 0; i < 4; i++){
        ibz_printf("%Zd ", &((*vec)[i]));
    }
    ibz_printf("\n");
}


static void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016llx", b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016llx", b.im[i]);
    printf(", ");
}

static void point_print(char *name, ec_point_t P){
    fp2_t a;
    if(fp2_is_zero(&P.z)){
        printf("%s = INF, ", name);
    }
    else{
    fp2_copy(&a, &P.z);
    fp2_inv(&a);
    fp2_mul(&a, &a, &P.x);
    fp2_print(name, a);
    }
}

static void curve_print(char *name, ec_curve_t E){
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, a);
}

void print_signature(const signature_t *sig) {
    fp2_t j;
    ec_j_inv(&j, &sig->E_com);
    fp2_print("j_E1 = ", j);
    // ibz_mat_2x2_print(&sig->mat_sigma_phichall);
    ibz_printf("M_sigma[00] = %Zd, ", &((sig->mat_sigma_phichall)[0][0]));
    ibz_printf("M_sigma[01] = %Zd, ", &((sig->mat_sigma_phichall)[0][1]));
    ibz_printf("M_sigma[10] = %Zd, ", &((sig->mat_sigma_phichall)[1][0]));
    ibz_printf("M_sigma[11] = %Zd", &((sig->mat_sigma_phichall)[1][1]));

}

void print_public_key(const public_key_t *pk) {
    fp2_t j;
    ec_j_inv(&j, &pk->curve);
    fp2_print("j_EA = ", j);
}



void commit(ec_curve_t *E_com, ec_basis_t *basis_even_com, quat_left_ideal_t *lideal_commit_three, int verbose) {
    
    quat_alg_elem_t gamma;
    quat_left_ideal_t lideal_even;
    ec_isog_even_t two_isogeny_first_half, two_isogeny_second_half;
    ec_isog_odd_t phi_first_half, phi_second_half;
    ec_point_t list_points[3];

    quat_alg_elem_init(&gamma);
    quat_left_ideal_init(&lideal_even); 
    doublepath(&gamma, &lideal_even, lideal_commit_three, 
    NULL,  // not used ?
    basis_even_com, E_com, verbose); // used only for image of BASIS_EVEN


    // #ifndef NDEBUG 
    //     ec_curve_t E_test;
    //     copy_curve(&E_test, &CURVE_E0);
    //     copy_point(list_points + 0, &BASIS_EVEN.P);
    //     copy_point(list_points + 1, &BASIS_EVEN.Q);
    //     copy_point(list_points + 2, &BASIS_EVEN.PmQ);
    //     TAC("ec_eval_odd in");
    //     ec_eval_odd(&E_test, &phi_first_half, list_points, 3);
    //     ec_eval_odd(&E_test, &phi_second_half, list_points, 3);
    //     TAC("ec_eval_odd out");

    //     assert(ec_is_equal(&(basis_even_com->P), list_points + 0));
    //     assert(ec_is_equal(&(basis_even_com->Q), list_points + 1));
    //     assert(ec_is_equal(&(basis_even_com->PmQ), list_points + 2));


    //     fp2_t j_R,j_L;
    //     ec_j_inv(&j_R, &E_test);
    //     ec_j_inv(&j_L, E_com);
    //     assert(fp2_is_equal(&j_R,&j_L));
    // #endif

    quat_alg_elem_finalize(&gamma);
    quat_left_ideal_finalize(&lideal_even); 
    return;
}



void quat_lideal_conjugate_lattice(quat_lattice_t *lat, const quat_left_ideal_t *lideal) {
    ibz_mat_4x4_copy(&(lat->basis), &(lideal->lattice.basis));
    ibz_copy(&(lat->denom), &(lideal->lattice.denom));
    
    for (int row = 1; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            ibz_neg(&(lat->basis[row][col]),&(lat->basis[row][col]));
        }
    }

    return;
}

int is_good_norm(ibz_t *N) {
    if ((8 - (ibz_get(N) % 8)) != 5) return 0;

    ibz_t pow2, sum_of_squares_candidate;
    ibz_init(&pow2);
    ibz_init(&sum_of_squares_candidate);
    int res = 0;

    ibz_set(&pow2, 1);
    ibz_mul_2exp(&pow2, &pow2, RESPONSE_LENGTH);

    if(ibz_cmp(&pow2, N) < 0) {
        ibz_printf("WARNING: short vectors not short enough...\n2-pow = %Zd\nnorm = %Zd\n", &pow2, &N);
        // assert(0);
        ibz_finalize(&sum_of_squares_candidate);
        return 0;
    }

    ibz_sub(&sum_of_squares_candidate, &pow2, N);

    // unsigned int N_mod_four = ibz_mod_ui (&N, 4);
    assert(ibz_mod_ui (&sum_of_squares_candidate, 8) == 5);

    // if (N_mod_four == 1) {
    res = ibz_probab_prime(&sum_of_squares_candidate, 40);

    ibz_finalize(&pow2);
    ibz_finalize(&sum_of_squares_candidate);
    return res;
}


int is_good(quat_alg_elem_t *x, ibz_t const *lattice_content) {
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
        if(ibz_cmp(&pow2, &N) < 0) {
            ibz_printf("WARNING: short vectors not short enough...\n2-pow = %Zd\nnorm = %Zd\n", &pow2, &N);
        }
        else {
            ibz_sub(&N, &pow2, &N);

            // unsigned int N_mod_four = ibz_mod_ui (&N, 4);
            unsigned int N_mod_eight = ibz_mod_ui (&N, 8);

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

void norm_from_2_times_gram(ibz_t *norm, ibz_mat_4x4_t *gram, ibz_vec_4_t *vec) {
    quat_qf_eval(norm, gram, vec);
    assert(ibz_is_even(norm));
    ibz_div_2exp(norm, norm, 1);
}

// TODO(security): currently just samples smallest vector, instead of random in a ball
void sample_response(quat_alg_elem_t *x, const quat_lattice_t *lattice, ibz_t const *lattice_content, int verbose) {
    ibz_mat_4x4_t lll;
    ibz_t denom_gram, norm;

    ibz_mat_4x4_init(&lll);
    ibz_init(&denom_gram);
    ibz_init(&norm);

    // printf("[");
    // for (int col = 0; col < 4; ++col) {
    //     printf("[");
    //     for (int row = 0; row < 4; ++row) {
    //         ibz_printf("%Zd ", &(lattice->basis[row][col]));
    //     }
    //     printf("]");
    // }
    // printf("]");

    int err = quat_lattice_lll(&lll, lattice, &(QUATALG_PINFTY.p), 1000);
    assert(!err);
    // The shortest vector found by lll is a candidate



    ibz_mat_4x4_t prod, gram;
    ibz_mat_4x4_init(&prod);
    ibz_mat_4x4_init(&gram);

    ibz_mat_4x4_transpose(&prod,&lll);
    ibz_mat_4x4_mul(&prod,&prod,&(QUATALG_PINFTY.gram));
    ibz_mat_4x4_mul(&gram,&prod,&lll);

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

    int k = 0;
    for (int i1 = 0; i1 < 10; i1++){
        for (int i2 = 0; i2 < 10; i2++){
            for (int i3 = 0; i3 < 10; i3++){
                for (int i4 = 0; i4 < 10; i4++){
                    k++;
                    ibz_vec_4_set(&vec, i1,i2,i3,i4);
                    norm_from_2_times_gram(&norm, &gram, &vec);
                    if (is_good_norm(&norm)) {

                        ibz_mat_4x4_eval(&(x->coord), &lll, &vec);

                        assert(quat_lattice_contains(NULL, lattice, x, &QUATALG_PINFTY));
                        assert(is_good(x, lattice_content));

                        #ifndef NDEBUG
                            ibq_t N_q;
                            ibz_t N, tmp, q;
                            ibq_init(&N_q);
                            ibz_init(&N);
                            ibz_init(&tmp);
                            quat_alg_norm(&N_q, x, &QUATALG_PINFTY);
                            ibq_to_ibz(&N, &N_q);
                            assert(ibz_divides(&N, lattice_content));
                            ibz_div(&N, &tmp, &N, lattice_content);

                            assert(ibz_cmp(&N, &norm) == 0);

                            ibz_finalize(&N);
                            ibq_finalize(&N_q);
                            ibz_finalize(&tmp);
                        #endif


                        #ifndef NDEBUG
                            printf("good sample found after %d attempts\n", k);
                        #endif

                        if (verbose) printf("e = %llu, ", RESPONSE_LENGTH);
                        if (verbose) ibz_printf("q = %Zd, ", &norm);
                        // ibz_printf("2^e - q = sum of two squares = prime 1 mod 4\n");
                        i1 = i2 = i3 = i4 = 10;

                    }



                }
            }
        }
    }
    // ibz_copy(&x->coord[i], &lll[i][0]);


    assert(quat_lattice_contains(NULL, lattice, x, &QUATALG_PINFTY));


    ibz_finalize(&denom_gram);
    ibz_finalize(&norm);
    ibz_mat_4x4_finalize(&prod);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&lll);
    ibz_vec_4_finalize(&vec);
    return;
}




// TODO(code): this is also used in verification, move to a common location when verification is implemented
void hash_to_challenge(ibz_vec_2_t *scalars, const ec_curve_t *curve, const unsigned char *message, const public_key_t *pk, size_t length)
{
    unsigned char *buf = malloc(sizeof(fp2_t) + sizeof(fp2_t) + length);
    {
        fp2_t j1, j2;
        ec_j_inv(&j1, curve);
        ec_j_inv(&j2, &pk->curve);
        memcpy(buf, &j1, sizeof(j1));
        memcpy(buf + sizeof(j1), &j2, sizeof(j2));
        memcpy(buf + sizeof(j1) + sizeof(j2), message, length);
    }

    //TODO(security) omit some vectors, notably (a,1) with gcd(a,6)!=1 but also things like (2,3)?
    {
        digit_t digits[NWORDS_FIELD];

        //FIXME should use SHAKE128 for smaller parameter sets?
        SHAKE256((void *) digits, sizeof(digits), buf, sizeof(fp2_t) + sizeof(fp2_t) + length);

        ibz_set(&(*scalars)[1], 1); //FIXME
        ibz_copy_digit_array(&(*scalars)[1], digits);
    }

    // #ifndef NDEBUG
    // {
    // ibz_t gcd;
    // ibz_init(&gcd);
    // ibz_set(&gcd, 6);
    // ibz_gcd(&gcd, &gcd, &(*scalars)[0]);
    // ibz_gcd(&gcd, &gcd, &(*scalars)[1]);
    // assert(ibz_is_one(&gcd));
    // ibz_finalize(&gcd);
    // }
    // #endif

    ibz_set(&((*scalars)[0]), 1);
    // ibz_rand_interval(&((*scalars)[1]), &((*scalars)[0]), &TORSION_PLUS_2POWER);

    free(buf);
}

int protocols_sign(signature_t *sig, const public_key_t *pk, const secret_key_t *sk, const unsigned char* m, size_t l, int verbose) {
    clock_t t = tic();

    ibz_t lattice_content;
    ec_curve_t E_com;
    ec_basis_t Bcom0, Bcom_can; // basis of 2^n-torsion
    ibz_vec_2_t vec, vec_can, vec_zero;
    quat_left_ideal_t lideal_tmp; 
    quat_left_ideal_t lideal_commit_three, lideal_chall_three; 
    quat_left_ideal_t lideal_chall3_secret2, lideal_chall3_secret3;
    quat_lattice_t lattice_hom_chall_to_com, lat_commit;
    quat_alg_elem_t resp_quat;
    quat_alg_elem_t elem_tmp;
    ibz_mat_2x2_t mat_alpha0, mat_Bcom0_to_Bcom;
    ibz_mat_2x2_t mat_sigma_phichall_BA_to_Bcomcan, mat_sigma_phichall_BA0_to_Bcom0;
    ibz_t degree_com_isogeny, tmp;

    ibz_init(&degree_com_isogeny); ibz_init(&tmp); ibz_init(&lattice_content);

    ibz_mat_2x2_init(&mat_alpha0); ibz_mat_2x2_init(&mat_Bcom0_to_Bcom); 
    ibz_mat_2x2_init(&mat_sigma_phichall_BA_to_Bcomcan); 
    ibz_mat_2x2_init(&mat_sigma_phichall_BA0_to_Bcom0);

    quat_alg_elem_init(&resp_quat);
    quat_alg_elem_init(&elem_tmp);
    quat_lattice_init(&lattice_hom_chall_to_com); quat_lattice_init(&lat_commit);
    quat_left_ideal_init(&lideal_tmp);
    quat_left_ideal_init(&lideal_commit_three); quat_left_ideal_init(&lideal_chall_three);
    quat_left_ideal_init(&lideal_chall3_secret2); quat_left_ideal_init(&lideal_chall3_secret3);


    ibz_vec_2_init(&vec); ibz_vec_2_init(&vec_can);

    // t = tic();
    commit(&E_com, &Bcom0, &lideal_commit_three, verbose);

    ibz_copy(&degree_com_isogeny, &(lideal_commit_three.norm));  
    
    // ibz_printf("degree_com_isogeny = %Zd ", &degree_com_isogeny);

    hash_to_challenge(&vec_can, &E_com, m, pk, l);
    // vec_can is a pair or random coefficients
    // the kernel of the challenge isogeny is generated by vec_can[0]*B[0] + vec_can[1]*B[1] where B is the canonical basis of the three^n torsion of EA



    ibz_mat_2x2_eval(&vec, &(sk->mat_BAcan_to_BA0_three), &vec_can);
    // the kernel of the challenge isogeny is generated by vec[0]*B0[0] + vec[1]*B0[1] where B0 is the image through secret isogeny of the canonical basis E0

    // t = tic();
    id2iso_kernel_dlogs_to_ideal_three(&lideal_chall_three, &vec);
    assert(ibz_cmp(&lideal_chall_three.norm, &TORSION_PLUS_3POWER) == 0);
    // TODO(optimization): only 3-torsion is used. Can optimise


    quat_lideal_inter(&lideal_chall3_secret2, &lideal_chall_three, &(sk->secret_ideal_two), &QUATALG_PINFTY);
    quat_lideal_mul(&lideal_chall3_secret3, &lideal_chall3_secret2, &(sk->two_to_three_transporter), &QUATALG_PINFTY, 0); 






    // Careful: want to intersect lideal_chall3_secret3 and dual(lideal_commit_three)
    // Both have norm a power of three, so trouble!! First replace lideal_chall3_secret3 with 
    // an equivalent ideal of norm coprime to 3; compute intersection with that, then sample 
    // in there, and transport the result back to the wanted intersection
    quat_lideal_generator_coprime(&elem_tmp, &lideal_chall3_secret3, &ibz_const_one, &QUATALG_PINFTY, 0);
    quat_alg_conj(&elem_tmp, &elem_tmp);
    ibz_mul(&(elem_tmp.denom), &(elem_tmp.denom) , &(lideal_chall3_secret3.norm));

    quat_lideal_mul(&lideal_tmp, &lideal_chall3_secret3, &elem_tmp, &QUATALG_PINFTY, 0); 
    int test = quat_lideal_isom(&elem_tmp, &lideal_tmp, &lideal_chall3_secret3, &QUATALG_PINFTY);
    assert(test);

    quat_lideal_conjugate_lattice(&lat_commit, &lideal_commit_three);


    quat_lattice_intersect(&lattice_hom_chall_to_com, &lideal_tmp.lattice, &lat_commit);
    // this lattice contains all isogenies that start with chall3_secret3 and end with dual(commit_three)

    // ibz_printf(">>>> lideal_chall3_secret3.norm = %Zd\n", &(lideal_chall3_secret3.norm));
    // ibz_printf(">>>> lideal_commit_three.norm = %Zd\n", &(lideal_commit_three.norm));

    ibz_mul(&lattice_content, &(lideal_tmp.norm), &(lideal_commit_three.norm));
    if (verbose) TOC(t, "sample_response in");

    sample_response(&resp_quat, &lattice_hom_chall_to_com, &lattice_content, verbose);
    assert(is_good(&resp_quat, &lattice_content));
    quat_alg_mul(&resp_quat, &resp_quat, &elem_tmp, &QUATALG_PINFTY); // bring it to intersection of lat_commit and lideal_chall3_secret3
    
    if (verbose) TOC(t, "sample_response out");


    #ifndef NDEBUG
    {
        ibq_t N_q;
        ibz_t N, tmp;

        ibq_init(&N_q);
        ibz_init(&N);
        ibz_init(&tmp);

        assert(quat_lattice_contains(NULL, &(lideal_chall3_secret3.lattice), &resp_quat, &QUATALG_PINFTY));
        assert(!quat_lattice_contains(NULL, &(lideal_commit_three.lattice), &resp_quat, &QUATALG_PINFTY));

        quat_alg_conj(&resp_quat, &resp_quat);
        assert(quat_lattice_contains(NULL, &(lideal_commit_three.lattice), &resp_quat, &QUATALG_PINFTY));
        assert(!quat_lattice_contains(NULL, &(lideal_chall3_secret3.lattice), &resp_quat, &QUATALG_PINFTY));
        quat_alg_conj(&resp_quat, &resp_quat); // repair

        quat_alg_norm(&N_q, &resp_quat, &QUATALG_PINFTY);
        ibq_to_ibz(&N, &N_q);
        // ibz_printf("norm of response: %Zd\n", &N);
    
        ibz_mul(&tmp, &(lideal_chall3_secret3.norm), &(lideal_commit_three.norm)); 
        assert(ibz_divides(&N, &tmp));
        assert(is_good(&resp_quat, &tmp));


        ibq_finalize(&N_q);
        ibz_finalize(&N);
        ibz_finalize(&tmp);
    }
    #endif



    // notational conventions:
    // BA = canonical basis of (the even torsion of) EA
    // B0 = canonical basis of E0
    // BA0 = image through secret isogeny of canonical basis of E0
    // Bcom0 = image through commitment isogeny (odd degree) of canonical basis of E0

    ec_curve_to_basis_2(&Bcom_can, &E_com);


    matrix_of_endomorphism_even(&mat_alpha0, &resp_quat); // matrix of alpha wrt the basis B0

    // printf("mat_alpha0 = ");
    // ibz_mat_2x2_print(&mat_alpha0);

    // M_Bcom_to_Bcom0*Bcom_can = Bcom0

    if (verbose) TOC(t, "change_of_basis_matrix_two in");
    change_of_basis_matrix_two(&mat_Bcom0_to_Bcom, &Bcom0, &Bcom_can, &E_com); // a 2-dimensional DLP
    if (verbose) TOC(t, "change_of_basis_matrix_two out");

    // (M_BA0_to_BA*v).BA = v.BA0, precomputed

    // compute mat_sigma_phichall_BA0_to_Bcom0 = mat_alpha0/degree_com_isogeny 
    // the matrix of sigma_phichall from basis BA0 to basis Bcom0
    ibz_invmod(&tmp, &degree_com_isogeny, &TORSION_PLUS_2POWER);
    for (int row = 0; row < 2; ++row) {
        for (int col = 0; col < 2; ++col) {
            ibz_mul(&(mat_sigma_phichall_BA0_to_Bcom0[row][col]), &(mat_alpha0[row][col]), &tmp);
            ibz_mod(&(mat_sigma_phichall_BA0_to_Bcom0[row][col]), &(mat_sigma_phichall_BA0_to_Bcom0[row][col]), &TORSION_PLUS_2POWER);
        }
    }
    // sigma_phichall_BA_to_Bcom = M_Bcom0_to_Bcom*M_sigma_phichall_BA0_to_Bcom0*M_BA_to_BA0 // the matrix of sigma_phichall from basis BA to basis Bcom

    ibz_2x2_mul_mod(&mat_sigma_phichall_BA_to_Bcomcan, &mat_Bcom0_to_Bcom, &mat_sigma_phichall_BA0_to_Bcom0, &TORSION_PLUS_2POWER);
    ibz_2x2_mul_mod(&mat_sigma_phichall_BA_to_Bcomcan, &mat_sigma_phichall_BA_to_Bcomcan, &(sk->mat_BAcan_to_BA0_two), &TORSION_PLUS_2POWER);

    copy_curve(&(sig->E_com), &E_com);
    ibz_mat_2x2_copy(&(sig->mat_sigma_phichall), &mat_sigma_phichall_BA_to_Bcomcan);




    ec_basis_t B_can_three, B_can_two, B_com_can_two;



    if (verbose) {
        curve_print("A_EA = ", pk->curve);

        ec_curve_to_basis_3(&B_can_three, &(pk->curve));

        point_print("xP3A = ", B_can_three.P);
        point_print("xQ3A = ", B_can_three.Q);
        point_print("xP3AmQ3A = ", B_can_three.PmQ);
        ibz_printf("ker_phi_vect[0] = %Zd, ", &(vec_can[0]));
        ibz_printf("ker_phi_vect[1] = %Zd, ", &(vec_can[1]));


        ec_curve_to_basis_2(&B_can_two, &(pk->curve)); 
        point_print("xPA = ", B_can_two.P);
        point_print("xQA = ", B_can_two.Q);
        point_print("xPAmQA = ", B_can_two.PmQ);


        curve_print("A_E1 = ", sig->E_com);

        ec_curve_to_basis_2(&B_com_can_two, &(sig->E_com)); 
        point_print("xP1 = ", B_com_can_two.P);
        point_print("xQ1 = ", B_com_can_two.Q);
        point_print("xP1mQ1 = ", B_com_can_two.PmQ);


        print_public_key(pk);
        print_signature(sig);
    }








    ibz_vec_2_finalize(&vec);
    ibz_vec_2_finalize(&vec_can);
    ibz_mat_2x2_finalize(&mat_alpha0); ibz_mat_2x2_finalize(&mat_Bcom0_to_Bcom); 
    ibz_mat_2x2_finalize(&mat_sigma_phichall_BA_to_Bcomcan); ibz_mat_2x2_finalize(&mat_sigma_phichall_BA0_to_Bcom0);

    quat_alg_elem_finalize(&resp_quat);
    quat_alg_elem_finalize(&elem_tmp);
    quat_lattice_finalize(&lattice_hom_chall_to_com); quat_lattice_finalize(&lat_commit);
    quat_left_ideal_finalize(&lideal_commit_three); quat_left_ideal_finalize(&lideal_chall_three);
    quat_left_ideal_finalize(&lideal_chall3_secret2); quat_left_ideal_finalize(&lideal_chall3_secret3);
    quat_left_ideal_finalize(&lideal_tmp);

    ibz_finalize(&degree_com_isogeny); ibz_finalize(&tmp); ibz_finalize(&lattice_content);
    return 0;
}






