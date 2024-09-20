#include <bench.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include <ec.h>
#include <gf_constants.h>

/******************************
Util functions
******************************/

static int BENCH_LOOPS = 500; // Number of iterations per bench

int
cmp_u64(const void *v1, const void *v2)
{
    uint64_t x1 = *(const uint64_t *)v1;
    uint64_t x2 = *(const uint64_t *)v2;
    if (x1 < x2) {
        return -1;
    } else if (x1 == x2) {
        return 0;
    } else {
        return 1;
    }
}

/******************************
Test functions
******************************/

int
inner_test_generated_basis(ec_basis_t *basis, ec_curve_t *curve, unsigned int n)
{
    int i;
    int PASSED = 1;

    ec_point_t P, Q;
    copy_point(&P, &basis->P);
    copy_point(&Q, &basis->Q);

    // Double points to get point of order 2
    for (i = 0; i < n - 1; i++) {
        xDBL_A24(&P, &P, &curve->A24);
        xDBL_A24(&Q, &Q, &curve->A24);
    }
    if (ec_is_zero(&P)) {
        printf("Point P generated does not have full order\n");
        PASSED = 0;
    }
    if (ec_is_zero(&Q)) {
        printf("Point Q generated does not have full order\n");
        PASSED = 0;
    }
    if (ec_is_equal(&P, &Q)) {
        printf("Points P, Q are linearly dependent\n");
        PASSED = 0;
    }

    if (!fp2_is_zero(&Q.x)) {
        printf("Points Q is not above the Montgomery point\n");
        PASSED = 0;
    }

    // This should give the identity
    xDBL_A24(&P, &P, &curve->A24);
    xDBL_A24(&Q, &Q, &curve->A24);
    if (!ec_is_zero(&P)) {
        printf("Point P generated does not have order exactly 2^n\n");
        PASSED = 0;
    }
    
    if (!ec_is_zero(&Q)) {
        printf("Point Q generated does not have order exactly 2^n\n");
        PASSED = 0;
    }

    if (PASSED == 0) {
        printf("Test failed with n = %d\n", n);
    }

    return PASSED;
}

int
inner_test_hint_basis(ec_basis_t *basis, ec_basis_t *basis_hint)
{
    int PASSED = 1;

    if (!ec_is_equal(&basis->P, &basis_hint->P)) {
        printf("The points P do not match using the hint\n");
        PASSED = 0;
    }

    if (!ec_is_equal(&basis->Q, &basis_hint->Q)) {
        printf("The points Q do not match using the hint\n");
        PASSED = 0;
    }

    if (!ec_is_equal(&basis->PmQ, &basis_hint->PmQ)) {
        printf("The points PmQ do not match using the hint\n");
        PASSED = 0;
    }

    if (PASSED == 0) {
        printf("Test failed\n");
    }

    return PASSED;
}

int
inner_test_basis_3(ec_basis_t *basis, ec_curve_t *curve, unsigned int n)
{
    int i;
    int PASSED = 1;

    ec_point_t P, Q, A3;
    copy_point(&P, &basis->P);
    copy_point(&Q, &basis->Q);

    fp2_sub(&A3.z, &(curve->A24.x), &(curve->A24.z));
    fp2_copy(&A3.x, &(curve->A24.x));

    // Triple points to get point of order 3
    for (i = 0; i < n - 1; i++) {
        xTPL(&P, &P, &A3);
        xTPL(&Q, &Q, &A3);
    }
    if (ec_is_zero(&P)) {
        printf("Point P generated does not have full order\n");
        PASSED = 0;
    }
    if (ec_is_zero(&Q)) {
        printf("Point Q generated does not have full order\n");
        PASSED = 0;
    }
    if (ec_is_equal(&P, &Q)) {
        printf("Points P, Q are linearly dependent\n");
        PASSED = 0;
    }

    // This should give the identity
    xTPL(&P, &P, &A3);
    xTPL(&Q, &Q, &A3);
    if (!ec_is_zero(&P)) {
        printf("Point P generated does not have order exactly 3^n\n");
        PASSED = 0;
    }
    
    if (!ec_is_zero(&Q)) {
        printf("Point Q generated does not have order exactly 3^n\n");
        PASSED = 0;
    }

    if (PASSED == 0) {
        printf("Test failed with n = %d\n", n);
    }

    return PASSED;
}

/******************************
Test wrapper functions
******************************/

int
test_basis_generation(unsigned int n)
{
    ec_basis_t basis;
    ec_curve_t curve;

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis
    ec_curve_to_basis_2f(&basis, &curve, n);

    // Test result
    return inner_test_generated_basis(&basis, &curve, n);
}

int
test_basis_generation_with_hints(unsigned int n)
{
    ec_basis_t basis, basis_hint;
    ec_curve_t curve;
    ec_curve_init(&curve);

    uint8_t hint[2];

    int check_1, check_2;

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis with hints
    ec_curve_to_basis_2f_to_hint(&basis, hint, &curve, n);

    // Ensure the basis from the hint is good
    check_1 = inner_test_generated_basis(&basis, &curve, n);

    // Generate a basis using hints
    ec_curve_to_basis_2f_from_hint(&basis_hint, &curve, n, hint);

    // These two bases should be the same
    check_2 = inner_test_hint_basis(&basis, &basis_hint);

    return check_1 && check_2;
}

int
test_basis_generation_3()
{
    ec_basis_t basis;
    ec_curve_t curve;

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis
    ec_curve_to_basis_3(&basis, &curve);

    // Test result
    return inner_test_basis_3(&basis, &curve, POWER_OF_3);
}

void
bench_basis_generation(unsigned int n)
{
    int i, j;
    uint64_t cycles1, cycles2;
    uint64_t cycle_runs[20];

    ec_basis_t basis;
    ec_curve_t curve;
    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Full even torsion generation without hints
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (j = 0; j < BENCH_LOOPS; j++) {
            ec_curve_to_basis_2f(&basis, &curve, n);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  2^%d torsion generation takes .................................... %" PRIu64
           " cycles\n",
           n,
           cycle_runs[4] / (BENCH_LOOPS));
}

void
bench_basis_generation_from_hint(unsigned int n)
{
    int i, j;
    uint64_t cycles1, cycles2;
    uint64_t cycle_runs[20];

    ec_basis_t basis;
    ec_curve_t curve;
    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    uint8_t hint[2];
    ec_curve_to_basis_2f_to_hint(&basis, hint, &curve, n);

    // Full even torsion generation without hints
    for (i = 0; i < 20; i++) {
        cycles1 = cpucycles();
        for (j = 0; j < BENCH_LOOPS; j++) {
            ec_curve_to_basis_2f_from_hint(&basis, &curve, n, hint);
        }
        cycles2 = cpucycles();
        cycle_runs[i] = cycles2 - cycles1;
    }
    qsort(cycle_runs + 10, 10, sizeof cycle_runs[0], cmp_u64);
    printf("  2^%d torsion generation takes .................................... %" PRIu64
           " cycles\n",
           n,
           cycle_runs[4] / (BENCH_LOOPS));
}

int
test_basis(void)
{
    int passed;

    // Test full order 2-torsion
    passed = test_basis_generation(POWER_OF_2);
    passed &= test_basis_generation_with_hints(POWER_OF_2);

    // Test partial order 2-torsion
    passed &= test_basis_generation(POWER_OF_2>>1);
    passed &= test_basis_generation_with_hints(POWER_OF_2>>1);

    // Test 3-torsion (full order)
    printf("Starting 3-torsion tests.\n");
    passed &= test_basis_generation_3();

    return passed;
}

void
bench_basis(void)
{
    printf("\n-------------------------------------------------------------------------------------"
           "-------------------\n\n");
    printf("Benchmarking E[2^n] basis generation: \n\n");
    bench_basis_generation(POWER_OF_2);
    bench_basis_generation(POWER_OF_2>>1);

    printf("\nBenchmarking E[2^n] basis generation with hint: \n\n");
    bench_basis_generation_from_hint(POWER_OF_2);
    bench_basis_generation_from_hint(POWER_OF_2>>1);
}

int
main(void)
{
    bool ok;
    ok = test_basis();
    /*
        bench_basis();
    */
    if (!ok) {
        printf("Tests failed!\n");
    } else {
        printf("All basis generation tests passed.\n");
    }
    return !ok;
}
