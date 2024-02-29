#include "id2iso_tests.h"

// TODO deduplicate with id2ker_odd.c


static void random_scalar(ibz_t *k)
{
    ibz_pow(k, &ibz_const_two, 123);
    ibz_rand_interval(k, &ibz_const_zero, k);
}


// defined in src/klpt/ref/generic/tools.c
int represent_integer(quat_alg_elem_t *gamma, ibz_t *n_gamma, const quat_alg_t *Bpoo);

// both ideals will have norm dividing T
static void random_pair_of_equivalent_ideals(quat_left_ideal_t *I, quat_left_ideal_t *J)
{
    ibz_t normIJ;
    ibz_init(&normIJ);
    ibz_mul(&normIJ, &TORSION_ODD, &TORSION_ODD);

// arbitrarily divide out one of the prime powers
{
ibz_t rem;
ibz_init(&rem);
ibz_div(&normIJ, &rem, &normIJ, &TORSION_ODD_PRIMEPOWERS[rand() % sizeof(TORSION_ODD_PRIMEPOWERS)/sizeof(*TORSION_ODD_PRIMEPOWERS)]);
assert(ibz_is_zero(&rem));
ibz_finalize(&rem);
}

    quat_alg_elem_t genI;
    quat_alg_elem_init(&genI);

    {
        int r = represent_integer(&genI, &normIJ, &QUATALG_PINFTY);
        assert(r);
    }

//   quat_alg_elem_print(&genI);

    quat_alg_elem_t genJ;
    quat_alg_elem_init(&genJ);
    quat_alg_conj(&genJ, &genI);

    ibz_t normI, normJ;
    ibz_init(&normI);
    ibz_init(&normJ);
    ibz_gcd(&normI, &normIJ, &TORSION_ODD);
    {
        ibz_t rem;
        ibz_init(&rem);
        ibz_div(&normJ, &rem, &normIJ, &normI);
        assert(ibz_is_zero(&rem));
        ibz_finalize(&rem);
    }

    quat_lideal_create_from_primitive(I, &genI, &normI, &MAXORD_O0, &QUATALG_PINFTY);
    quat_lideal_create_from_primitive(J, &genJ, &normJ, &MAXORD_O0, &QUATALG_PINFTY);

    ibz_finalize(&normJ);
    ibz_finalize(&normI);
    quat_alg_elem_finalize(&genJ);
    quat_alg_elem_finalize(&genI);
    ibz_finalize(&normIJ);
}


int _id2iso_test_id2kerdlogs_odd()
{
    int res = 1;

    quat_left_ideal_t ideal, ignored;
    quat_left_ideal_init(&ideal);
    {
        quat_left_ideal_init(&ignored);
        random_pair_of_equivalent_ideals(&ignored, &ideal);
        quat_left_ideal_finalize(&ignored);
    }
//    quat_left_ideal_print(&ideal);

    ibz_vec_2_t vec1, vec2;
    ibz_vec_2_init(&vec1);
    ibz_vec_2_init(&vec2);

    // call it twice and make sure we get equivalent results
    void id2iso_ideal_to_kernel_dlogs_odd(ibz_vec_2_t *vec, ec_degree_odd_t *deg, const quat_left_ideal_t *lideal);
    ec_degree_odd_t deg1, deg2;
    id2iso_ideal_to_kernel_dlogs_odd(&vec1, &deg1, &ideal);
    id2iso_ideal_to_kernel_dlogs_odd(&vec2, &deg2, &ideal);

//    ibz_printf("vec1 = (%Zd, %Zd)\n", &vec1[0], &vec1[1]);
//    ibz_printf("vec2 = (%Zd, %Zd)\n", &vec2[0], &vec2[1]);
//    printf("deg1 = "); for (size_t i = 0; i < sizeof(deg1)/sizeof(*deg1); ++i) printf(",%d"+(!i), (int) deg1[i]); printf("\n");
//    printf("deg2 = "); for (size_t i = 0; i < sizeof(deg2)/sizeof(*deg2); ++i) printf(",%d"+(!i), (int) deg2[i]); printf("\n");

    ibz_t lhs, rhs;
    ibz_init(&lhs);
    ibz_init(&rhs);
    ibz_mul(&lhs, &vec1[0], &vec2[1]);
    ibz_mul(&rhs, &vec2[0], &vec1[1]);
    ibz_mod(&lhs, &lhs, &ideal.norm);
    ibz_mod(&rhs, &rhs, &ideal.norm);
//   ibz_printf("lhs = %Zd\n", &lhs);
//   ibz_printf("rhs = %Zd\n", &rhs);
    res &= !ibz_cmp(&lhs, &rhs);
    ibz_finalize(&lhs);
    ibz_finalize(&rhs);

    quat_left_ideal_finalize(&ideal);

    ibz_vec_2_finalize(&vec1);
    ibz_vec_2_finalize(&vec2);

    return res;
}

#include <endomorphism_action.h>

int _id2iso_test_id2ker_odd()
{
    int res = 1;

    ec_isog_odd_t isog1, isog2;
    {
        quat_left_ideal_t ideal1, ideal2;
        quat_left_ideal_init(&ideal1);
        quat_left_ideal_init(&ideal2);

        random_pair_of_equivalent_ideals(&ideal1, &ideal2);

//        quat_left_ideal_print(&ideal1);
//        quat_left_ideal_print(&ideal2);

        id2iso_ideal_to_isogeny_odd(&isog1, &CURVE_E0, &BASIS_ODD_PLUS, &BASIS_ODD_MINUS, &ideal1);
        id2iso_ideal_to_isogeny_odd(&isog2, &CURVE_E0, &BASIS_ODD_PLUS, &BASIS_ODD_MINUS, &ideal2);

        // check orders
        {
            if (!ibz_is_one(&ideal1.norm) && fp2_is_zero(&isog1.ker_plus.z) && fp2_is_zero(&isog1.ker_minus.z)) {
                printf("norm of ideal1 is >1 but kernel is trivial; something is wrong!\n");
                res &= 0;
            }
            if (!ibz_is_one(&ideal2.norm) && fp2_is_zero(&isog2.ker_plus.z) && fp2_is_zero(&isog2.ker_minus.z)) {
                printf("norm of ideal2 is >1 but kernel is trivial; something is wrong!\n");
                res &= 0;
            }
        }

        quat_left_ideal_finalize(&ideal1);
        quat_left_ideal_finalize(&ideal2);
    }

    {
        ec_curve_t curve1, curve2;
        ec_eval_odd(&curve1, &isog1, NULL, 0);
        ec_eval_odd(&curve2, &isog2, NULL, 0);

        if (fp2_is_zero(&curve1.C) || fp2_is_zero(&curve2.C)) {
            printf("broken curve constant after isogeny; something is wrong!\n");
            res &= 0;
        }

        fp2_t j1, j2;
        ec_j_inv(&j1, &curve1);
        ec_j_inv(&j2, &curve2);

        res &= fp2_is_equal(&j1, &j2);
    }

    return res;
}


int id2iso_test_id2ker_odd() {
    int res = 1;
    printf("\n \nRunning id2iso tests for ideal_to_kernel_odd() \n \n");

    for (int i = 0; i < 5; i++) {
        res &= _id2iso_test_id2kerdlogs_odd();
        res &= _id2iso_test_id2ker_odd();
    }

    if (!res) {
        printf("ID2ISO unit test ideal_to_kernel_odd() failed\n");
    }

    return res;
}
