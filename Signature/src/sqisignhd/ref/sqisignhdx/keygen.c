#include <sqisignhd.h>
#include <curve_extras.h>

#include <inttypes.h>

void
public_key_init(public_key_t *pk)
{
    pk->hint_pk = (int *)malloc(2 * sizeof(int));
    ec_point_init(&pk->curve.A24);
    pk->curve.is_A24_computed_and_normalized = 0;
}

void
public_key_finalize(public_key_t *pk)
{
    free(pk->hint_pk);
}

void
secret_key_init(secret_key_t *sk)
{
    quat_left_ideal_init(&(sk->secret_ideal));
    ibz_mat_2x2_init(&(sk->mat_BAcan_to_BA0_two));
}

void
secret_key_finalize(secret_key_t *sk)
{
    quat_left_ideal_finalize(&(sk->secret_ideal));
    ibz_mat_2x2_finalize(&(sk->mat_BAcan_to_BA0_two));
}

void 
fprint_public_key(FILE *p_file,const public_key_t *pk)
{
    fprintf(p_file, "A_pk = ");
    fp2_print_to_file(p_file,&(pk->curve.A));
    fprintf(p_file, "hint_pk_P = %u\n",pk->hint_pk[0]);
    fprintf(p_file, "hint_pk_Q = %u\n",pk->hint_pk[1]);
}

void
protocols_keygen(public_key_t *pk, secret_key_t *sk)
{

    int found = 1;

    ibz_t n; //replaced by FIXED_DEGREE_SK
    ec_basis_t B_can_two, B_0_two;

    ibz_init(&n); //replaced by FIXED_DEGREE_SK
    ibz_copy(&n,&FIXED_DEGREE_SK);

    // generate a random ideal of random norm for the secret ideal
    // TODO make a clean function for all of that and
    //generate_random_prime(&n, 1, ibz_bitsize(&QUATALG_PINFTY.p) / 2); //replaced by FIXED_DEGREE_SK
    sampling_random_ideal_O0(&sk->secret_ideal, &n, 1);

    // ideal to isogeny clapotis
    found = dim2id2iso_arbitrary_isogeny_evaluation(&B_0_two, &sk->curve, &sk->secret_ideal);
    assert(found);

    assert(test_point_order_twof(&(B_0_two.P), &(sk->curve), TORSION_PLUS_EVEN_POWER));

    ec_curve_to_basis_2f_to_hint(
        &B_can_two, &(sk->curve), TORSION_PLUS_EVEN_POWER, pk->hint_pk); // canonical
    assert(test_point_order_twof(&(B_can_two.P), &(sk->curve), TORSION_PLUS_EVEN_POWER));

    copy_point(&sk->canonical_basis.P, &B_can_two.P);
    copy_point(&sk->canonical_basis.Q, &B_can_two.Q);
    copy_point(&sk->canonical_basis.PmQ, &B_can_two.PmQ);

    // !!! Notation should be mat_BA0_to_BAcan
    change_of_basis_matrix_two(
        &(sk->mat_BAcan_to_BA0_two), &B_can_two, &B_0_two, &(sk->curve), TORSION_PLUS_EVEN_POWER);

    // computing the public key
    fp2_t temp;
    fp2_copy(&temp, &sk->curve.C);
    fp2_inv(&temp);
    fp2_set_one(&pk->curve.C);
    fp2_mul(&pk->curve.A, &temp, &sk->curve.A);

    ibz_finalize(&n); //replaced by FIXED_DEGREE_SK

    return;
}
