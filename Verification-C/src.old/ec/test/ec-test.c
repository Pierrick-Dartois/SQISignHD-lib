//#include "ec-tests.h"
//Replaced by
#include "test_arith_data.h"
#include "../../gf/test/test_extras.h"
#include <stdio.h>

bool ec_test(){
    // Tests for ecc arithmetic
    bool OK = true;

    ec_point_t R = {0}, S = {0};

    fp2_tomont(&AC.z, &AC.z);
        
    fp2_tomont(&R.x, &P.x);
    fp2_tomont(&R.z, &P.z);
    xDBL(&S, &R, &AC); // S=2P=twoP
    //fp2_copy(&SS.x, &S.x);    // Copy of S = SS <- 2P 
    //fp2_copy(&SS.z, &S.z);
    fp2_inv(&S.z);
    fp2_mul(&S.x, &S.x, &S.z);
    fp2_frommont(&S.x, &S.x);

    if (compare_words((digit_t*)&S.x, (digit_t*)&twoP.x, NWORDS_FIELD*2)!=0) { passed=0; goto out0; }

out0:
    if (passed==1) printf("  ECC arithmetic tests ............................................ PASSED");
    else { printf("  ECC arithmetic tests... FAILED"); printf("\n"); return false; }
    printf("\n");
 
    return OK;
}

int main(int argc, char* argv[])
{
    if (argc < 3) {
        printf("Please enter an argument: 'test' or 'bench' and <reps>\n");
        exit(1);
    }
    if (!strcmp(argv[1], "test")) {
        TEST_LOOPS = atoi(argv[2]);
        return ec_test();
    } else if (!strcmp(argv[1], "bench")) {
        BENCH_LOOPS = atoi(argv[2]);
        exit(1);
        //return !(ec_run() & dlog_run());
    } else {
        exit(1);
    }
}