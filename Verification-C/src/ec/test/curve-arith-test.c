#include <bench.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>

#include "test_extras.h"
#include <ec.h>

/******************************
Util functions
******************************/

static int BENCH_LOOPS = 500; // Number of iterations per bench
static int TEST_LOOPS = 500; // Number of iterations per bench

/******************************
Test functions
******************************/

int
test_xDBL_xADD(ec_curve_t *curve, unsigned int Ntest)
{
    int i;

    ec_point_t P, Q, PQ, R1, R2;

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        ec_random_test(&Q, curve);
        projective_difference_point(&PQ, &P, &Q, curve);

        // 2(P + Q) = 2P + 2Q
        xADD(&R1, &P, &Q, &PQ);
        ec_dbl(&R1, &R1, curve);
        ec_dbl(&P, &P, curve);
        ec_dbl(&Q, &Q, curve);
        ec_dbl(&PQ, &PQ, curve);
        xADD(&R2, &P, &Q, &PQ);
        if(!ec_is_equal(&R1, &R2)){
            printf("Failed 2(P + Q) = 2P + 2Q\n");
            return 1;
        }

        // (P+Q) + (P-Q) = 2P
        xADD(&R1, &P, &Q, &PQ);
        ec_dbl(&Q, &Q, curve);
        xADD(&R1, &R1, &PQ, &Q);
        ec_dbl(&P, &P, curve);
        ec_dbl(&PQ, &PQ, curve);
        if(!ec_is_equal(&R1, &P)){
            printf("Failed (P+Q) + (P-Q) = 2P\n");
            return 1;
        }
    }

    return 0;
}

int
test_xDBLADD(ec_curve_t *curve, unsigned int Ntest)
{
    int i;

    ec_point_t P, Q, PQ, R1, R2;

    ec_point_t A24;
    AC_to_A24(&A24, curve);

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        ec_random_test(&Q, curve);
        projective_difference_point(&PQ, &P, &Q, curve);
        
        xDBLADD(&R1, &R2, &P, &Q, &PQ, &A24);
        xADD(&PQ, &P, &Q, &PQ);
        if(!ec_is_equal(&R2, &PQ)){
            printf("Failed addition in xDBLADD\n");
            return 1;
        }
        ec_dbl(&P, &P, curve);
        if(!ec_is_equal(&R1, &P)){
            printf("Failed doubling in xDBLADD\n");
            return 1;
        }
    }
    return 0;
}

int
test_xDBL_variants(ec_curve_t *curve, unsigned int Ntest)
{
    int i;

    ec_point_t P, R1, R2, R3;
    fp2_t z;

    ec_point_t A24, A24norm;
    AC_to_A24(&A24, curve);
    copy_point(&A24norm, &A24);
    ec_normalize_point(&A24norm);

    // Randomize projective representation
    fp2_random_test(&z);
    fp2_mul(&(A24.x), &(A24.x), &z);
    fp2_mul(&(A24.z), &(A24.z), &z);


    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);
        xDBL(&R1, &P, (const ec_point_t *)curve);
        xDBL_A24(&R2, &P, &A24);
        xDBL_A24_normalized(&R3, &P, &A24norm);
        if(!ec_is_equal(&R1, &R2)){
            printf("xDBL and xDBL_A24 dont match\n");
            return 1;
        }
        if(!ec_is_equal(&R1, &R3)){
            printf("xDBL and xDBL_A24_normalized dont match\n");
            return 1;
        }
    }
    return 0;
}

int
test_zero_identities(ec_curve_t *curve, unsigned int Ntest)
{
    int i;

    ec_point_t P, Q, R, ec_zero;

    fp2_set_one(&(P.x));
    fp2_set_zero(&(P.z));

    assert(ec_is_zero(&P));

    for (i = 0; i < Ntest; i++) {
        ec_random_test(&P, curve);


        
        xADD(&R, &ec_zero, &ec_zero, &ec_zero);
        if(!ec_is_zero(&R)){
            printf("Failed 0 + 0 = 0\n");
            return 1;
        }
        
        ec_dbl(&R, &P, curve);
        xADD(&R, &P, &P, &R);
        if(!ec_is_zero(&R)){
            printf("Failed P - P = 0\n");
            return 1;
        }
    
        // TODO: current formulas can't handle 2*0 nor P+0
        // ec_dbl(&R, &ec_zero, curve);
        // if(!ec_is_zero(&R)){
        //     printf("Failed 2*0 = 0\n");
        //     return 1;
        // }

        // xADD(&R, &P, &ec_zero, &P);
        // if(!ec_is_equal(&R, &P)){
        //     printf("Failed P + 0 = P\n");
        //     return 1;
        // }
        // xADD(&R, &ec_zero, &P, &P);
        // if(!ec_is_equal(&R, &P)){
        //     printf("Failed P + 0 = P\n");
        //     return 1;
        // }

        // xDBLADD(&R, &Q, &P, &ec_zero, &P, &(curve->A24));
        // if(!ec_is_equal(&Q, &P)){
        //     printf("Failed P + 0 = P in xDBLADD\n");
        //     return 1;
        // }
        // xDBLADD(&R, &Q, &ec_zero, &P, &P, &(curve->A24));
        // if(!ec_is_equal(&Q, &P)){
        //     printf("Failed P + 0 = P in xDBLADD\n");
        //     return 1;
        // }
        // if(!ec_is_zero(&R)){
        //     printf("Failed 2*0 = 0 in xDBLADD\n");
        //     return 1;
        // }
    }
    return 0;
}


int
main()
{
    bool fail = 0;
    
    // Curve A=6
    ec_curve_t curve;
    ec_curve_init(&curve);
    fp2_set_small(&(curve.A), 6);
    fp2_random_test(&(curve.C));
    fp2_mul(&(curve.A), &(curve.A), &(curve.C));
    ec_curve_normalize_A24(&curve);

    fail |= test_xDBL_xADD(&curve, TEST_LOOPS);
    fail |= test_xDBLADD(&curve, TEST_LOOPS);
    fail |= test_xDBL_variants(&curve, TEST_LOOPS);
    fail |= test_zero_identities(&curve, TEST_LOOPS);

    if (fail) {
        printf("Tests failed!\n");
    } else {
        printf("All ec arithmetic tests passed.\n");
    }
    return fail;
}