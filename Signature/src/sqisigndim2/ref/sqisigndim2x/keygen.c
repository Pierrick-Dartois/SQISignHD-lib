#include <sqisigndim2.h>
#include <curve_extras.h>

#include <inttypes.h>

void
public_key_init(public_key_t *pk)
{
    pk->hint_pk = (int *)malloc(2 * sizeof(int));
    ec_curve_init(&pk->curve);
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
    ec_curve_init(&sk->curve);
}

void
secret_key_finalize(secret_key_t *sk)
{
    quat_left_ideal_finalize(&(sk->secret_ideal));
    ibz_mat_2x2_finalize(&(sk->mat_BAcan_to_BA0_two));
}

void
protocols_keygen(public_key_t *pk, secret_key_t *sk)
{

    int found = 1;

    ibz_t n;
    ec_basis_t B_can_two, B_0_two;

    ibz_init(&n);

    // generate a random ideal of random norm for the secret ideal
    // TODO make a clean function for all of that and
    generate_random_prime(&n, 1, ibz_bitsize(&QUATALG_PINFTY.p));
    sampling_random_ideal_O0(&sk->secret_ideal, &n, 1);

    // ideal to isogeny clapotis
    found = dim2id2iso_arbitrary_isogeny_evaluation(&B_0_two, &sk->curve, &sk->secret_ideal);
    assert(found);

    assert(test_point_order_twof(&(B_0_two.P), &(sk->curve), TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_0_two.Q), &(sk->curve), TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_0_two.PmQ), &(sk->curve), TORSION_PLUS_EVEN_POWER));

    ec_curve_to_basis_2f_to_hint(
        &B_can_two, &(sk->curve), TORSION_PLUS_EVEN_POWER, pk->hint_pk); // canonical
    assert(test_point_order_twof(&(B_can_two.P), &(sk->curve), TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_can_two.Q), &(sk->curve), TORSION_PLUS_EVEN_POWER));
    assert(test_point_order_twof(&(B_can_two.PmQ), &(sk->curve), TORSION_PLUS_EVEN_POWER));

    copy_point(&sk->canonical_basis.P, &B_can_two.P);
    copy_point(&sk->canonical_basis.Q, &B_can_two.Q);
    copy_point(&sk->canonical_basis.PmQ, &B_can_two.PmQ);

    change_of_basis_matrix_two(
        &(sk->mat_BAcan_to_BA0_two), &B_can_two, &B_0_two, &(sk->curve), TORSION_PLUS_EVEN_POWER);

    // computing the public key
    fp2_t temp;
    fp2_copy(&temp, &sk->curve.C);
    fp2_inv(&temp);
    fp2_set_one(&pk->curve.C);
    fp2_mul(&pk->curve.A, &temp, &sk->curve.A);

    ibz_finalize(&n);

    return;
}
