//#include "ec-tests.h"
//Replaced by
#include "../include/ec.h"
#include "test_arith_data.h"
#include "../include/curve_extras.h"
//#include "../../gf/test/test_extras.h"
#include <stdio.h>

bool ec_test(){
    ec_point_t R = {0}, S={0}, T={0}, U={0}, A24={0};

    xDBL(&R, &P, &AC);
    xADD(&S, &P, &Q, &PmQ);

    ec_curve_t E0;
    point_to_curve(&E0,&AC);

    xMUL(&T, &P, k, &E0);

    curve_to_A24(&A24, &E0);
    
    xMULv2(&U, &P, k, NWORDS_FIELD*RADIX, &A24);

    return is_point_equal(&R,&twoP)&&is_point_equal(&S,&PpQ)&&is_point_equal(&T,&kP)&&is_point_equal(&U,&kP);
}

int main(int argc, char* argv[])
{
    bool b=ec_test();
    printf("Are all tests correct ? %d\n",b);
}