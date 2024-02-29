
#include <ec.h>
#include <quaternion.h>
#include <sqisignhd.h>
#include <inttypes.h>
#include "test_sqisignhd.h"


static void fp2_print(char *name, fp2_t const a){
    fp2_t b;
    fp2_set(&b, 1);
    fp2_mul(&b, &b, &a);
    printf("%s = 0x", name);
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.re[i]);
    printf(" + i*0x");
    for(int i = NWORDS_FIELD - 1; i >=0; i--)
        printf("%016" PRIx64, b.im[i]);
    printf("\n");
}

static void curve_print(char *name, ec_curve_t E){
    fp2_t a;
    fp2_copy(&a, &E.C);
    fp2_inv(&a);
    fp2_mul(&a, &a, &E.A);
    fp2_print(name, a);
}

bool curve_is_canonical(ec_curve_t const *E);   // test_sqisign.c

int test_doublepath()
{
    int res = 1;

    quat_left_ideal_t ideal;
    quat_left_ideal_init(&ideal);

    ec_curve_t E1;
    ec_basis_t E1basis = BASIS_CHALLENGE;

    res &= 1;

    quat_left_ideal_finalize(&ideal);

    return res;
}

