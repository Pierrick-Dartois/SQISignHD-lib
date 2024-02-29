#include <rng.h>
#include <stdio.h>
#include <ec.h>
#include <inttypes.h>

#include "test_sqisignhd.h"
#include "toolbox.h"




//XXX FIXME stolen from src/ec/opt/generic/test/isog-test.c
void fp2_print(char *name, fp2_t const a){
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

static void point_print(char *name, ec_point_t P){
    fp2_t a;
    if(fp2_is_zero(&P.z)){
        printf("%s = INF\n", name);
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
//XXX


bool curve_is_canonical(ec_curve_t const *E)
{
    ec_curve_t EE;
    ec_isom_t isom;
    ec_curve_normalize(&EE, &isom, E);

    fp2_t lhs, rhs;
    fp2_mul(&lhs, &E->A, &EE.C);
    fp2_mul(&rhs, &E->C, &EE.A);
    return fp2_is_equal(&lhs, &rhs);
}

int test_sqisign(int repeat)
{
    int res = 1;

    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    unsigned char msg[32] = { 0 };

    secret_key_init(&sk);
    secret_sig_init(&sig);

    clock_t t = tic();
    protocols_keygen(&pk, &sk);
    TOC(t, "protocols_keygen");


    printf("Printing details of first signature\n");
    int val = protocols_sign(&sig, &pk, &sk, msg, 32, 1);


    printf("\n\nTesting more signatures\n");
    t = tic();
    for (int i = 0; i < repeat; ++i)
    {
        int val = protocols_sign(&sig, &pk, &sk, msg, 32, 0);
    }
    // TOC(t, "protocols_sign");

    float ms = (1000. * (float) (clock() - t) / CLOCKS_PER_SEC);
    printf("average signing time [%.2f ms]\n", (float) (ms/repeat));


    secret_key_finalize(&sk);

    return res;
}

// run all tests in module
int main(){
    int res = 1;

    randombytes_init((unsigned char *) "some", (unsigned char *) "string", 128);

    // printf("\nRunning encoding tests\n");
    // res &= test_encode();

    printf("\nRunning sqisignhd tests\n \n");

    printf("Format of printed data:\n\n");
    printf("A_EA = Montgomery coefficient of public key,\n");
    printf("xP3A, xQ3A, and xP3AmQ3A = cannonical differential basis of the power-of-3 torsion of the public key, \n");
    printf("ker_phi_vect[0] and ker_phi_vect[1] = coefficients defining the challenge, \n");
    printf("xPA, xQA, xPAmQA = cannonical differential basis of the power-of-2 torsion of the public key, \n");
    printf("A_E1 = Montgomery coefficient of commitment, \n");
    printf("xP1, xQ1, xP1mQ1 = cannonical differential basis of the power-of-2 torsion of the commitment, \n");
    printf("j_EA, j_E1 = j-invariants of the public key and commitment, \n");
    printf("M_sigma[ij] = matrix defining the response isogeny with respect to the cannonical bases of the power-of-2 torsion\n\n");

    res &= test_sqisign(10);

    if(!res){
        printf("\nSome tests failed!\n");
    } 
    else {
        printf("All tests passed!\n");
    }
    return(!res);
}
