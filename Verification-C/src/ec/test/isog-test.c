#include <isog.h>
#include <strategies.h>
#include <ec.h>
#include <mp.h>
#include <stdio.h>
#include <inttypes.h>
#include <gf_constants.h>
#include <assert.h>


int ec_2_isog_chain_test(unsigned int n){
	int passed = 1;

	ec_basis_t basis;
	ec_point_t R, phiP, phiQ, phiR, psiP, psiQ, lambphiQ;
    ec_curve_t curve;
    unsigned int nwords=NWORDS_FIELD/2;
    digit_t lambda[nwords], one[nwords];
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
    for(int i=0;i<nwords;i++){
    	lambda[i]=rand();
    }

    // Compute R=P+[lambda]Q, never above (0,0).
    ec_biscalar_mul(&R,one,lambda,nwords*RADIX,&basis,&curve);

    // Strategy
    optimised_strategy(strategy,n, M+2*S, 4*M);

    // Compute isogeny with kernel <R>.
    ec_2_isog_chain(&phi,&R,&curve,n,strategy);

    // Evaluate phi(P), phi(Q), phi(R) 
    ec_eval_2_isog_chain(&phiP,&basis.P,&phi);
    ec_eval_2_isog_chain(&phiQ,&basis.Q,&phi);
    ec_eval_2_isog_chain(&phiR,&R,&phi);
    
    uint16_t n_bits=mp_nbits(lambda,nwords);

    ec_mul(&lambphiQ,lambda,n_bits,&phiQ,&phi.codomain);
    
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

int ec_odd_isog_chain_test(unsigned int n){
    int passed = 1;

    ec_basis_t basis;
    ec_point_t A3, R, phiP, phiQ, phiR, lambphiQ, ker, tmp, test, A24;
    ec_kps_t kps;
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

    // Evaluate phi(P), phi(Q), phi(R) 
    ec_eval_odd_isog_chain(&phiP,&basis.P,&phi);
    ec_eval_odd_isog_chain(&phiQ,&basis.Q,&phi);
    ec_eval_odd_isog_chain(&phiR,&R,&phi);

    ec_mul(&lambphiQ,lambda,nwords*RADIX,&phiQ,&phi.codomain);

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

	printf("Testing 2-isogeny chain with full torsion.\n");
	ok=ok&ec_2_isog_chain_test(POWER_OF_2);
	printf("Testing 2-isogeny chain with half torsion.\n");
	ok=ok&ec_2_isog_chain_test(POWER_OF_2>>1);
    printf("Testing 3-isogeny chain with full torsion.\n");
    ok=ok&ec_odd_isog_chain_test(POWER_OF_3);


	if(ok){
		printf("All tests passed.\n");
	}
	else{
		printf("Tests failed.\n");
	}

}

