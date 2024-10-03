#include <isog.h>
#include <strategies.h>
#include <ec.h>
#include <mp.h>
#include <stdio.h>
#include <inttypes.h>
#include <gf_constants.h>

int ec_2_isog_chain_test(unsigned int n){
	int passed = 1;

	ec_basis_t basis;
	ec_point_t R, phiQ, phiR, psiP, psiQ, dblsP, dblsQ, A24, B24, rhodblsP, rhodblsQ;
    ec_curve_t curve, E1;
    unsigned int nwords=NWORDS_FIELD/2;
    digit_t lambda[nwords], one[nwords];
    ec_2_isog_chain_t phi, psi;
    ec_kps2_t rho;
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

    AC_to_A24(&A24,&curve);

    copy_point(&dblsP,&basis.P);
    copy_point(&dblsQ,&basis.Q);
    for(int i=0;i<n-1;i++){
        xDBL_A24(&dblsP,&dblsP,&A24);
        xDBL_A24(&dblsQ,&dblsQ,&A24);
    }

    printf("%u\n",ec_is_zero(&dblsP));
    printf("%u\n",ec_is_zero(&dblsQ));

    xisog_2(&rho, &B24, dblsP);
    A24_to_AC(&E1,&B24);
    xeval_2(&rhodblsP, &dblsP, 1, &rho);
    xeval_2(&rhodblsQ, &dblsQ, 1, &rho);
    printf("rho(P)=0 ? %u\n",ec_is_zero(&rhodblsP));
    printf("rho(Q) is on curve ? %u\n",ec_is_on_curve(&rhodblsQ,&E1));
    printf("%u\n",ec_is_zero(&rhodblsQ));


    // Strategy
    optimised_strategy(strategy,n, M+2*S, 4*M);

    // Compute isogeny with kernel <R>.
    ec_2_isog_chain(&phi,&R,&curve,n,strategy);

    // Evaluate phi(Q), phi(R) 
    ec_eval_2_isog_chain(&phiQ,&basis.Q,&phi);
    ec_eval_2_isog_chain(&phiR,&R,&phi);

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

    if(!ec_is_on_curve(&psiP,&psi.codomain)){
    	printf("psi(P) not on image curve.\n");
    	passed = 0;
    }

    if(!ec_is_zero(&psiQ)){
    	printf("psi(Q)!=0.\n");
    	passed = 0;
    }

    return passed;

}

int main(){
	int ok=1;

	printf("Testing 2-isogeny chain with full torsion.\n");
	ok=ok&ec_2_isog_chain_test(POWER_OF_2);
	printf("Testing 2-isogeny chain with half torsion.\n");
	ok=ok&ec_2_isog_chain_test(POWER_OF_2>>1);

	if(ok){
		printf("All tests passed.\n");
	}
	else{
		printf("Tests failed.\n");
	}

}

