#include <bench.h>
#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include <ec.h>
#include <biextension.h>
#include <gf_constants.h>
#include <test_utils.h>

int test_weil_pairing_2f(unsigned int n){
	int passed = 1;

	fp2_t wPQ, wQP, w2PQ, wP2Q, wP2P, wPQ_sqr, wPQwQP, wPQ_2f;
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

    fp2_copy(&wPQ_2f,&wPQ_sqr);
    for(int i=1; i<POWER_OF_2-1;i++){
        fp2_sqr(&wPQ_2f,&wPQ_2f);
    }

    if(!fp2_is_one(&wPQwQP)){
    	printf("Pairing is not skew-symmetric.\n");
    	passed = 0;
    }

    if(!fp2_is_one(&wP2P)){
    	printf("Self-pairing (up to scalar) is not trivial.\n");
    	passed = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&w2PQ)){
    	printf("Pairing is not left-linear.\n");
    	passed = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&wP2Q)){
    	printf("Pairing is not right-linear.\n");
    	passed = 0;
    }

    if(fp2_is_one(&wPQ_2f)){
        printf("Pairing does not have full order.\n");
        passed = 0;
    }

    fp2_sqr(&wPQ_2f,&wPQ_2f);
    if(!fp2_is_one(&wPQ_2f)){
        printf("Pairing is not a 2^f-th root of unity.\n");
        passed = 0;
    }

    return passed;
}

int test_weil_pairing_3g(){
	int passed = 1;

	fp2_t wPQ, wQP, w2PQ, wP2Q, wP2P, wPQ_sqr, wPQwQP, wPQ_3g, temp;
	ec_point_t P, Q, PQ, dblP, dblQ, tplP, dblPQ, PdblQ;
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
    xADD(&tplP, &dblP, &P, &P);

    int nbits=mp_nbits(n, nwords);

    weil(&wPQ, n, nwords, &P, &Q, &PQ, &curve);
    weil(&wQP, n, nwords, &Q, &P, &PQ, &curve);
    weil(&w2PQ, n, nwords, &dblP, &Q, &dblPQ, &curve);
    weil(&wP2Q, n, nwords, &P, &dblQ, &PdblQ, &curve);
    weil(&wP2P, n, nwords, &P, &dblP, &tplP, &curve);

    fp2_sqr(&wPQ_sqr,&wPQ);
    fp2_mul(&wPQwQP,&wPQ,&wQP);

    fp2_mul(&wPQ_3g,&wPQ_sqr,&wPQ);
    for(int i=1; i<POWER_OF_3-1;i++){
        fp2_copy(&temp,&wPQ_3g);
        fp2_sqr(&wPQ_3g,&wPQ_3g);
        fp2_mul(&wPQ_3g,&wPQ_3g,&temp);
    }

    if(!fp2_is_one(&wPQwQP)){
    	printf("Pairing is not skew-symmetric.\n");
    	passed = 0;
    }

    if(!fp2_is_one(&wP2P)){
    	printf("Self-pairing (up to scalar) is not trivial.\n");
    	passed = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&w2PQ)){
    	printf("Pairing is not left-linear.\n");
    	passed = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&wP2Q)){
    	printf("Pairing is not right-linear.\n");
    	passed = 0;
    }

    if(fp2_is_one(&wPQ_3g)){
        printf("Pairing does not have full order.\n");
        passed = 0;
    }

    fp2_copy(&temp,&wPQ_3g);
    fp2_sqr(&wPQ_3g,&wPQ_3g);
    fp2_mul(&wPQ_3g,&wPQ_3g,&temp);
    if(!fp2_is_one(&wPQ_3g)){
        printf("Pairing is not a 3^g-th root of unity.\n");
        passed = 0;
    }

    return passed;
}

int test_weil_pairing_2f3g(){
    int passed = 1;

    fp2_t wPQ, wQP, w2PQ, wP2Q, wP2P, wPQ_sqr, wPQwQP, wPQ_2f3g, temp;
    ec_point_t P, Q, PQ, dblP, dblQ, tplP, dblPQ, PdblQ;
    ec_basis_t basis2, basis3;
    jac_point_t xyP2, xyP3, xyQ2, xyQ3, xyP2P3, xyQ2Q3, xyPQ;
    ec_curve_t curve;
    const unsigned int nwords=NWORDS_FIELD;
    uint64_t n[nwords],m[nwords],one[nwords];

    mp_set_small(one,1,nwords);

    mp_set_small(n,3,nwords);
    for(int i=0;i<POWER_OF_3-1;i++){
        mp_add(m,n,n,nwords);
        mp_add(n,m,n,nwords);
    }
    multiple_mp_shiftl(n,POWER_OF_2,nwords);

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis
    ec_curve_to_basis_2f(&basis2, &curve, POWER_OF_2);
    ec_curve_to_basis_3(&basis3, &curve);

    lift_basis(&xyP2, &xyQ2, &basis2, &curve);
    lift_basis(&xyP3, &xyQ3, &basis3, &curve);
    
    ADD(&xyP2P3, &xyP2, &xyP3, &curve);
    jac_to_xz(&P, &xyP2P3);

    ADD(&xyQ2Q3, &xyQ2, &xyQ3, &curve);
    jac_to_xz(&Q, &xyQ2Q3);

    ADD(&xyPQ, &xyP2P3, &xyQ2Q3, &curve);
    jac_to_xz(&PQ, &xyPQ);

    xADD(&dblPQ, &PQ, &P, &Q);
    xADD(&PdblQ, &PQ, &Q, &P);
    xDBL_A24_normalized(&dblP, &P, &curve.A24);
    xDBL_A24_normalized(&dblQ, &Q, &curve.A24);
    xADD(&tplP, &dblP, &P, &P);

    int nbits=mp_nbits(n, nwords);

    weil(&wPQ, n, nwords, &P, &Q, &PQ, &curve);
    weil(&wQP, n, nwords, &Q, &P, &PQ, &curve);
    weil(&w2PQ, n, nwords, &dblP, &Q, &dblPQ, &curve);
    weil(&wP2Q, n, nwords, &P, &dblQ, &PdblQ, &curve);
    weil(&wP2P, n, nwords, &P, &dblP, &tplP, &curve);

    fp2_sqr(&wPQ_sqr,&wPQ);
    fp2_mul(&wPQwQP,&wPQ,&wQP);

    fp2_mul(&wPQ_2f3g,&wPQ_sqr,&wPQ);
    for(int i=1; i<POWER_OF_3;i++){
        fp2_copy(&temp,&wPQ_2f3g);
        fp2_sqr(&wPQ_2f3g,&wPQ_2f3g);
        fp2_mul(&wPQ_2f3g,&wPQ_2f3g,&temp);
    }
    for(int i=0; i<POWER_OF_2-1;i++){
        fp2_sqr(&wPQ_2f3g,&wPQ_2f3g);
    }

    if(!fp2_is_one(&wPQwQP)){
        printf("Pairing is not skew-symmetric.\n");
        passed = 0;
    }

    if(!fp2_is_one(&wP2P)){
        printf("Self-pairing (up to scalar) is not trivial.\n");
        passed = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&w2PQ)){
        printf("Pairing is not left-linear.\n");
        passed = 0;
    }

    if(!fp2_is_equal(&wPQ_sqr,&wP2Q)){
        printf("Pairing is not right-linear.\n");
        passed = 0;
    }

    if(fp2_is_one(&wPQ_2f3g)){
        printf("Pairing does not have full order.\n");
        passed = 0;
    }

    fp2_sqr(&wPQ_2f3g,&wPQ_2f3g);
    if(!fp2_is_one(&wPQ_2f3g)){
        printf("Pairing is not a 2^f3^g-th root of unity.\n");
        passed = 0;
    }

    return passed;
}

int test_multiple_weil_dlog_2f(unsigned int n){
    int passed = 1;

    ec_point_t R, S, RmS, test;
    ec_basis_t PQ, RS;
    ec_curve_t curve;
    const unsigned int nwords=n/RADIX+!!(n%RADIX);
    uint64_t r1[nwords], r2[nwords], s1[nwords], s2[nwords], cr1[nwords], cr2[nwords], cs1[nwords], cs2[nwords], r1ms1[nwords], r2ms2[nwords];

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis
    ec_curve_to_basis_2f(&PQ, &curve, n);

    // Compute r1, r2, s1, s2, r1ps1, r2ps2
    for(int i=0;i<nwords;i++){
        r1[i]=rand();
        r2[i]=rand();
        s1[i]=rand();
        s2[i]=rand();
    }
    mp_sub(r1ms1,r1,s1,nwords);
    mp_sub(r2ms2,r2,s2,nwords);

    // Compute R, S
    xDBLMUL(&R,&PQ.P,r1,&PQ.Q,r2,&PQ.PmQ,nwords*RADIX,&curve);
    xDBLMUL(&S,&PQ.P,s1,&PQ.Q,s2,&PQ.PmQ,nwords*RADIX,&curve);
    xDBLMUL(&RmS,&PQ.P,r1ms1,&PQ.Q,r2ms2,&PQ.PmQ,nwords*RADIX,&curve);

    copy_point(&RS.P,&R);
    copy_point(&RS.Q,&S);
    copy_point(&RS.PmQ,&RmS);

    ec_dlog_2_weil(cr1,cr2,cs1,cs2,&PQ,&RS,&curve,n,nwords);

    xDBLMUL(&test,&PQ.P,cr1,&PQ.Q,cr2,&PQ.PmQ,nwords*RADIX,&curve);

    if(ec_is_equal(&test,&R)==0){
        printf("Wrong value for R.\n");
        passed = 0;
    }

    ec_biscalar_mul(&test, cs1, cs2, n, &PQ, &curve);

    if(ec_is_equal(&test,&S)==0){
        printf("Wrong value for S.\n");
        passed = 0;
    }

    return passed;

}

int test_signle_weil_dlog_3g(){
    int passed = 1;

    fp2_t wPQ, wRP, wRQ;
    ec_point_t P, Q, PQ, R, test;
    ec_basis_t basis;
    ec_curve_t curve;
    const unsigned int nwords=NWORDS_FIELD/2+1;
    uint64_t n[nwords], m[nwords], one[nwords], r1[nwords], r2[nwords], cr1[nwords], cr2[nwords];

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

    for(int i=0;i<nwords;i++){
        r1[i]=rand();
        r2[i]=rand();
    }

    xDBLMUL(&R,&P,r1,&Q,r2,&basis.PmQ,nwords*RADIX,&curve);

    ec_single_dlog_le_weil(cr1,cr2,&basis,&R,&curve,3,POWER_OF_3,nwords);

    xDBLMUL(&test,&basis.P,cr1,&basis.Q,cr2,&basis.PmQ,nwords*RADIX,&curve);

    if(ec_is_equal(&test,&R)==0){
        printf("Wrong value for R.\n");
        passed = 0;
    }

    return passed;
}

int
main(void)
{
    bool ok;
    printf("Testing full 2^f-torsion basis.\n");
    ok = test_weil_pairing_2f(POWER_OF_2);
    printf("Testing full 3^g-torsion basis.\n");
    ok = ok&test_weil_pairing_3g();
    printf("Testing full 2^f*3^g-torsion basis.\n");
    ok = ok&test_weil_pairing_2f3g();
    printf("Testing full 2^f-torsion multiple discrete log.\n");
    ok = ok&test_multiple_weil_dlog_2f(POWER_OF_2);
    printf("Testing full 3^g-torsion multiple discrete log.\n");
    ok = ok&test_signle_weil_dlog_3g();

    if (!ok) {
        printf("Tests failed!\n");
    } else {
        printf("All pairing tests passed.\n");
    }
    return !ok;
}
