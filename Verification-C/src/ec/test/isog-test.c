#include <isog.h>
#include <strategies.h>
#include <ec.h>
#include <mp.h>
#include <stdio.h>
#include <inttypes.h>
#include <gf_constants.h>
#include <assert.h>


int ec_2_isog_chain_test(digit_t *lambda, unsigned int n){
	int passed = 1;

	ec_basis_t basis;
	ec_point_t R, phiP, phiQ, phiR, psiP, psiQ, lambphiQ;
    ec_curve_t curve;
    unsigned int nwords=NWORDS_FIELD/2;
    //digit_t lambda[nwords], one[nwords];
    digit_t one[nwords];
    ec_2_isog_chain_t phi, psi;
    unsigned int strategy[n-1];
    float M=1.0, S=0.8; 

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    // Generate a basis (P,Q) with [2^(n-1)]Q=(0,0)
    ec_curve_to_basis_2f(&basis, &curve, n);

    mp_set_small(one,1,nwords);
    //for(int i=0;i<nwords;i++){
    	//lambda[i]=rand();
    //}

    // Compute R=P+[lambda]Q, never above (0,0).
    ec_biscalar_mul(&R,one,lambda,nwords*RADIX,&basis,&curve);

    // Strategy
    optimised_strategy(strategy,n, M+2*S, 4*M);

    // Compute isogeny with kernel <R>.
    for(int i=0;i<NWORDS_FIELD;i++){
        printf("R.x.re[%i]=%llu\n",i,R.x.re[i]);
    }
    ec_2_isog_chain(&phi,&R,&curve,n,strategy);

    // Evaluate phi(P), phi(Q), phi(R) 
    ec_eval_2_isog_chain(&phiP,&basis.P,&phi);
    ec_eval_2_isog_chain(&phiQ,&basis.Q,&phi);
    ec_eval_2_isog_chain(&phiR,&R,&phi);
    
    uint16_t n_bits=mp_nbits(lambda,nwords);

    ec_mul(&lambphiQ,lambda,n_bits,&phiQ,&phi.codomain);
    //xMUL(&lambphiQ,
                    //&phiQ,
                    //lambda,
                    //n_bits,
                    //&phi.codomain);

    printf("%llu\n",phiP.x.re[0]);
    printf("%llu\n",phiQ.x.re[0]);
    printf("%llu\n",lambda[0]);
    printf("%llu\n",lambda[1]);
    printf("%i\n",n_bits);
    printf("%llu\n",phi.codomain.A.re[0]);
    printf("%llu\n",phi.codomain.C.re[0]);

    
    // Compute isogeny with kernel <Q>.
    ec_2_isog_chain(&psi,&basis.Q,&curve,n,strategy);

    // Evaluate psi(P), psi(Q) 
    ec_eval_2_isog_chain(&psiP,&basis.P,&psi);
    ec_eval_2_isog_chain(&psiQ,&basis.Q,&psi);

    if(!ec_is_on_curve(&phiQ,&phi.codomain)){
    	printf("phi(Q) not on image curve.\n");
    	passed = 0;
    }

    if(!ec_is_zero(&phiR)){
    	printf("phi(R)!=0.\n");
    	passed = 0;
    }

    if(!ec_is_equal(&lambphiQ,&phiP)){
        printf("[lambda]phi(Q)!=phi(P).\n");
        passed = 0;
    }


    if(!ec_is_on_curve(&psiP,&psi.codomain)){
    	printf("psi(P) not on image curve.\n");
    	passed = 0;
    }

    if(!ec_is_zero(&psiQ)){
    	printf("psi(Q)!=0.\n");
    	passed = 0;
    }

    del_2_isog_chain(&phi);
    del_2_isog_chain(&psi);

    return passed;

}

int ec_2_isog_chain_test_2(digit_t *lambda, unsigned int n){
    int passed = 1;

    ec_basis_t basis2;
    ec_point_t S, phi2P, phi2Q, phi2S, psi2P, psi2Q, lambphi2Q;
    ec_curve_t curve2;
    unsigned int nwords=NWORDS_FIELD/2;
    //digit_t lambda[nwords], one[nwords];
    digit_t one[nwords];
    ec_2_isog_chain_t phi2, psi2;
    unsigned int strategy2[n-1];
    float M=1.0, T=0.8; 

    ec_curve_init(&curve2);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve2.A), 6);
    fp2_set_one(&(curve2.C));
    ec_curve_normalize_A24(&curve2);

    // Generate a basis (P,Q) with [2^(n-1)]Q=(0,0)
    ec_curve_to_basis_2f(&basis2, &curve2, n);

    mp_set_small(one,1,nwords);
    //for(int i=0;i<nwords;i++){
        //lambda[i]=rand();
    //}

    // Compute R=P+[lambda]Q, never above (0,0).
    ec_biscalar_mul(&S,one,lambda,nwords*RADIX,&basis2,&curve2);


    // Strategy
    optimised_strategy(strategy2,n, M+2*T, 4*M);

    // Compute isogeny with kernel <R>.
    for(int i=0;i<NWORDS_FIELD;i++){
        printf("S.x.re[%i]=%llu\n",i,S.x.re[i]);
    }
    ec_2_isog_chain(&phi2,&S,&curve2,n,strategy2);

    // Evaluate phi(P), phi(Q), phi(R) 
    ec_eval_2_isog_chain(&phi2P,&basis2.P,&phi2);
    ec_eval_2_isog_chain(&phi2Q,&basis2.Q,&phi2);
    ec_eval_2_isog_chain(&phi2S,&S,&phi2);
    
    uint16_t n_bits=mp_nbits(lambda,nwords);

    ec_mul(&lambphi2Q,lambda,n_bits,&phi2Q,&phi2.codomain);
    //xMUL(&lambphi2Q,
                    //&phi2Q,
                    //lambda,
                    //n_bits,
                    //&phi2.codomain);
    
    // Compute isogeny with kernel <Q>.
    //ec_2_isog_chain(&psi2,&basis2.Q,&curve2,n,strategy2);

    // Evaluate psi(P), psi(Q) 
    //ec_eval_2_isog_chain(&psi2P,&basis2.P,&psi2);
    //ec_eval_2_isog_chain(&psi2Q,&basis2.Q,&psi2);

    printf("%llu\n",phi2P.x.re[0]);
    printf("%llu\n",phi2Q.x.re[0]);
    printf("%llu\n",lambda[0]);
    printf("%llu\n",lambda[1]);
    printf("%i\n",n_bits);
    printf("%llu\n",phi2.codomain.A.re[0]);
    printf("%llu\n",phi2.codomain.C.re[0]);

    if(!ec_is_on_curve(&phi2Q,&phi2.codomain)){
        printf("phi(Q) not on image curve.\n");
        passed = 0;
    }

    if(!ec_is_zero(&phi2S)){
        printf("phi(R)!=0.\n");
        passed = 0;
    }

    if(!ec_is_equal(&lambphi2Q,&phi2P)){
        printf("[lambda]phi(Q)!=phi(P).\n");
        passed = 0;
    }

    //printf("%u\n",ec_is_equal(&basis2.P,&basis2.P));


    //if(!ec_is_on_curve(&psi2P,&psi2.codomain)){
        //printf("psi(P) not on image curve.\n");
        //passed = 0;
    //}

    //if(!ec_is_zero(&psi2Q)){
        //printf("psi(Q)!=0.\n");
        //passed = 0;
    //}

    del_2_isog_chain(&phi2);
    //del_2_isog_chain(&psi2);

    return passed;

}

int ec_odd_isog_chain_test(unsigned int n){
    int passed = 1;

    ec_basis_t basis;
    ec_point_t A3, R, phiP, phiQ, phiR, lambphiQ;
    ec_curve_t curve, E1;
    unsigned int nwords=NWORDS_FIELD/2;
    digit_t lambda[nwords], one[nwords];
    ec_odd_isog_chain_t phi, psi;
    unsigned int strategy[n-1];
    float M=1.0, S=0.8;
    fp2_t t0; 

    ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    fp2_set_small(&(curve.A), 6);
    fp2_set_one(&(curve.C));
    ec_curve_normalize_A24(&curve);

    fp2_add(&t0,&curve.C,&curve.C); // 2C
    fp2_add(&A3.x,&curve.A,&t0); // A+2C
    fp2_sub(&A3.z,&curve.A,&t0); // A-2C

    // Generate a basis of 3^g-torsion
    ec_curve_to_basis_3(&basis, &curve);

    for(int i=0;i<POWER_OF_3-n;i++){
        xTPL(&basis.P,&basis.P,&A3);
        xTPL(&basis.Q,&basis.Q,&A3);
        xTPL(&basis.PmQ,&basis.PmQ,&A3);
    }

    mp_set_small(one,1,nwords);
    for(int i=0;i<nwords;i++){
        lambda[i]=rand();
    }

    // Compute R=P+[lambda]Q.
    int n_bits=mp_nbits(lambda,nwords);

    ec_biscalar_mul(&R,one,lambda,n_bits,&basis,&curve);


    // Strategy
    optimised_strategy(strategy,n, 4*M+2*S, 7*M+5*S);

    // Compute isogeny with kernel <R>.
    ec_odd_isog_chain(&phi,&R,&curve,3,n,strategy);

    // Evaluate phi(Q), phi(R) 
    ec_eval_odd_isog_chain(&phiP,&basis.P,&phi);
    ec_eval_odd_isog_chain(&phiQ,&basis.Q,&phi);
    ec_eval_odd_isog_chain(&phiR,&R,&phi);

    ec_mul(&lambphiQ,lambda,nwords*RADIX,&phiQ,&phi.codomain);

    // Compute isogeny with kernel <Q>.
    //ec_odd_isog_chain(&phi,&R,&curve,3,n,strategy,false);

    // Evaluate psi(P), psi(Q) 
    //ec_eval_odd_isog_chain(&phiQ,&basis.Q,&phi,false);
    //ec_eval_odd_isog_chain(&phiR,&R,&phi,false);
    printf("%llu\n",phi.codomain.A.re[0]);

    if(!ec_is_on_curve(&phiQ,&phi.codomain)){
        printf("phi(Q) not on image curve.\n");
        passed = 0;
    }

    if(!ec_is_zero(&phiR)){
        printf("phi(R)!=0.\n");
        passed = 0;
    }

    if(!ec_is_equal(&lambphiQ,&phiP)){
        printf("[lambda]phi(Q)!=phi(P).\n");
        passed = 0;
    }

    //if(!ec_is_on_curve(&psiP,&psi.codomain)){
        //printf("psi(P) not on image curve.\n");
        //passed = 0;
    //}

    //if(!ec_is_zero(&psiQ)){
        //printf("psi(Q)!=0.\n");
        //passed = 0;
    //}

    del_odd_isog_chain(&phi);

    return passed;

}

int main(){
	int ok=1;

    ec_curve_t curve;
    ec_basis_t basis;
    ec_point_t R, S;
    unsigned int nwords=NWORDS_FIELD/2;
    digit_t lambda[nwords], mu[nwords], one[nwords];

    mp_set_small(one,1,nwords);
    for(int i=0;i<nwords;i++){
        lambda[i]=rand();
        mu[i]=rand();
    }

    mu[0]=282475249;
    mu[1]=984943658;
    lambda[0]=16807;
    lambda[1]=1622650073;


    //ec_curve_init(&curve);

    // Set a supersingular elliptic curve
    // E : y^2 = x^3 + 6*x^2 + x
    //fp2_set_small(&(curve.A), 6);
    //fp2_set_one(&(curve.C));
    //ec_curve_normalize_A24(&curve);


    //ec_curve_to_basis_2f(&basis, &curve, POWER_OF_2);
    //ec_biscalar_mul(&R,one,lambda,nwords*RADIX,&basis,&curve);

    //ec_biscalar_mul(&R,one,lambda,nwords*RADIX,&basis,&curve);
    //ec_biscalar_mul(&S,one,mu,nwords*RADIX,&basis,&curve);

    //for(int i=0;i<NWORDS_FIELD;i++){
        //printf("R.x.re[%i]=%llu\n",i,R.x.re[i]);
    //}

    //for(int i=0;i<NWORDS_FIELD;i++){
        //printf("S.x.re[%i]=%llu\n",i,S.x.re[i]);
    //}

    for(int i=0;i<nwords;i++){
        printf("lambda[%i]=%llu\n",i,lambda[i]);
    }
    for(int i=0;i<nwords;i++){
        printf("mu[%i]=%llu\n",i,mu[i]);
    }


	printf("Testing 2-isogeny chain with full torsion.\n");
	ok=ok&ec_2_isog_chain_test(lambda,POWER_OF_2);
	printf("Testing 2-isogeny chain with half torsion.\n");
	ok=ok&ec_2_isog_chain_test_2(mu,POWER_OF_2);
    printf("Testing 3-isogeny chain with full torsion.\n");
    ok=ok&ec_odd_isog_chain_test(POWER_OF_3);


	if(ok){
		printf("All tests passed.\n");
	}
	else{
		printf("Tests failed.\n");
	}

}

