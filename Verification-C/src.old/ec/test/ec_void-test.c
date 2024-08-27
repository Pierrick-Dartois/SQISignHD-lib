//#include "ec-tests.h"
//Replaced by
#include "../include/ec.h"
#include "test_arith_data.h"
//#include "../include/curve_extras.h"
//#include "../../gf/test/test_extras.h"
#include <stdio.h>

bool ec_mont_test(){
    ec_point_t R = {0}, S={0}, T={0}, U={0}, V={0}, W={0}, A24={0};
    ec_curve_t E0;
    digit_t m[NWORDS_FIELD]={1};
    fp2_t j_inv={0}, j0={1728}, t;

    // R=2P
    xDBL(&R, &P, &AC);

    // S=P+Q
    xADD(&S, &P, &Q, &PmQ);

    point_to_curve(&E0,&AC);

    // T=k*P
    xMUL(&T, &P, k, &E0);

    curve_to_A24(&A24, &E0);

    // U=k*P
    xMULv2(&U, &P, k, NWORDS_FIELD*RADIX, &A24);

    // V=P+k*Q
    xDBLMUL(&V, &P, m, &Q, k, &PmQ, &E0);

    // W=P+k*Q
    ec_ladder3pt(&W, k, &P, &Q, &PmQ, &E0);

    // j(E0)
    ec_j_inv(&j_inv, &E0);
    fp2_sub(&t, &j0, &j_inv);
    
    return is_point_equal(&R,&twoP)&&is_point_equal(&S,&PpQ)&&is_point_equal(&T,&kP)&&is_point_equal(&U,&kP)&&is_point_equal(&V,&PpkQ)&&is_point_equal(&W,&PpkQ)&&fp2_is_zero(&t);
}

bool ec_jac_test(){
    jac_point_t R = {0}, S={0};
    ec_curve_t E0;

    point_to_curve(&E0,&AC);

    //fp2_tomont(&E0.A,&E0.A);
    //fp2_tomont(&E0.C,&E0.C);
    fp2_tomont(&jac_P.x,&jac_P.x);
    fp2_tomont(&jac_P.y,&jac_P.y);
    fp2_tomont(&jac_P.z,&jac_P.z);
    fp2_tomont(&R.x,&R.x);
    fp2_tomont(&R.y,&R.y);
    fp2_tomont(&R.z,&R.z);

    DBL(&R, &jac_P, &E0);

    fp2_frommont(&jac_P.x,&jac_P.x);
    fp2_frommont(&jac_P.y,&jac_P.y);
    fp2_frommont(&jac_P.z,&jac_P.z);
    fp2_frommont(&R.x,&R.x);
    fp2_frommont(&R.y,&R.y);
    fp2_frommont(&R.z,&R.z);
    
    return is_jac_equal(&R,&jac_twoP);
}

int main(int argc, char* argv[])
{
    bool b=ec_mont_test();
    printf("Are all ec tests correct ? %d\n",b);
    bool c=ec_jac_test();
    printf("Are all jac tests correct ? %d\n",c);
}