#include <bench.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include <ec.h>
#include <biextension.h>
#include <gf_constants.h>

int test_weil_pairing_2f(unsigned int n){
	int PASSED = 1;

	fp2_t wPQ, wQP, w2PQ, wP2Q, wP2P, wPQ_sqr, wPQwQP;
	ec_point_t P, Q, PQ, dblP, dblQ, dblPQ, PdblQ;
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

    copy_point(&P,&basis.P);
    copy_point(&Q,&basis.Q);
    xADD(&PQ,&basis.P,&basis.Q,&basis.PmQ);

    xADD(&dblPQ, &PQ, &P, &Q);
    xADD(&PdblQ, &PQ, &Q, &P);
    xDBL_A24_normalized(&dblP, &P, &curve.A24);
    xDBL_A24_normalized(&dblQ, &Q, &curve.A24);

    weil_2e(&wPQ, n, &P, &Q, &PQ, &curve);
    weil_2e(&wQP, n, &Q, &P, &PQ, &curve);
    weil_2e(&w2PQ, n, &dblP, &Q, &dblPQ, &curve);
    weil_2e(&wP2Q, n, &P, &dblQ, &PdblQ, &curve);
    weil_2e(&wP2P, n, &P, &dblP, &P, &curve);

    fp2_sqr(&wPQ_sqr,&wPQ);
    fp2_mul(&wPQwQP,&wPQ,&wQP);

    if(!fp2_is_one(&wPQwQP)){
    	printf("Pairing is not skew-symmetric.\n");
    	PASSED = 0;
    }

    if(!fp2_is_one(&wP2P)){
    	printf("Self-pairing (up to scalar) is not trivial.\n");
    	PASSED = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&w2PQ)){
    	printf("Pairing is not left-linear.\n");
    	PASSED = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&wP2Q)){
    	printf("Pairing is not right-linear.\n");
    	PASSED = 0;
    }

    return PASSED;
}

int test_weil_pairing_3g(){
	int PASSED = 1;

	fp2_t wPQ, wQP, w2PQ, wP2Q, wP2P, wPQ_sqr, wPQwQP;
	ec_point_t P, Q, PQ, dblP, dblQ, dblPQ, PdblQ, S;
	ec_basis_t basis;
    ec_curve_t curve;
    const unsigned int nwords=NWORDS_FIELD/2+1;
    uint64_t n[nwords],m[nwords],one[nwords];

    mp_set_small(one,1,nwords);

    mp_set_small(n,3,nwords);
    for(int i=0;i<POWER_OF_3-1;i++){
    	mp_add(m,n,n,nwords);
    	mp_add(n,m,n,nwords);
    }

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis
    ec_curve_to_basis_3(&basis, &curve);

    copy_point(&P,&basis.P);
    copy_point(&Q,&basis.Q);
    xADD(&PQ,&basis.P,&basis.Q,&basis.PmQ);

    xADD(&dblPQ, &PQ, &P, &Q);
    xADD(&PdblQ, &PQ, &Q, &P);
    xDBL_A24_normalized(&dblP, &P, &curve.A24);
    xDBL_A24_normalized(&dblQ, &Q, &curve.A24);

    int nbits=mp_nbits(n, nwords);

    xDBLMUL(&S,&P,one,&Q,n,&PQ,nbits,&curve);
    printf("%u\n",ec_is_equal(&S,&P));

    printf("wPQ\n");
    weil(&wPQ, n, nwords, &P, &Q, &PQ, &curve);
    printf("wQP\n");
    weil(&wQP, n, nwords, &Q, &P, &PQ, &curve);
    printf("w2PQ\n");
    weil(&w2PQ, n, nwords, &dblP, &Q, &dblPQ, &curve);
    printf("wP2Q\n");
    weil(&wP2Q, n, nwords, &P, &dblQ, &PdblQ, &curve);
    printf("wP2P\n");
    weil(&wP2P, n, nwords, &P, &dblP, &P, &curve);

    fp2_sqr(&wPQ_sqr,&wPQ);
    fp2_mul(&wPQwQP,&wPQ,&wQP);

    if(!fp2_is_one(&wPQwQP)){
    	printf("Pairing is not skew-symmetric.\n");
    	PASSED = 0;
    }

    if(!fp2_is_one(&wP2P)){
    	printf("Self-pairing (up to scalar) is not trivial.\n");
    	PASSED = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&w2PQ)){
    	printf("Pairing is not left-linear.\n");
    	PASSED = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&wP2Q)){
    	printf("Pairing is not right-linear.\n");
    	PASSED = 0;
    }

    return PASSED;
}

int
main(void)
{
    bool ok;
    ok = test_weil_pairing_2f(POWER_OF_2);
    ok = ok&test_weil_pairing_3g();

    if (!ok) {
        printf("Tests failed!\n");
    } else {
        printf("All pairing tests passed.\n");
    }
    return !ok;
}
