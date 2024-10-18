#include <dim2_change_basis.h>
#include <fp2.h>
#include <fp2_matrix.h>
#include <mp.h>

static inline void
choose_non_vanishing_index(int ir0, int ir1, fp2_t *v, const digit_t **M, const fp2_t *zeta, unsigned int nwords, bool cst_time){
	fp2_t zeta_powers[4];

	fp2_set_one(&zeta_powers[0]);
	fp2_copy(&zeta_powers[1],zeta);
	fp2_sqr(&zeta_powers[2],zeta);
	fp2_mul(&zeta_powers[3],&zeta_powers[2],zeta);

	if(cst_time){
		fp2_t w[4];
		int k0, k1, l0, l1, iw;
		bool not_all_zero;
		for(int i0=0;i0<2;i0++){
			for(int i1=0;i1<2;i1++){
				w[0]=0;
				w[1]=0;
				w[2]=0;
				w[3]=0;
				for(int j0=0;j0<2;j0++){
					for(int j1=0;j1<2;j1++){
						k0 = M[2][0] * j0 + C[2][1] * j1;
                		k1 = C[3][0] * j0 + C[3][1] * j1;

                		l0 = D[2][0] * j0 + D[2][1] * j1;
                		l1 = D[3][0] * j0 + D[3][1] * j1;

                		e = (-(k0 + 2 * i0) * l0 - (k1 + 2 * i1) * l1) % 4;
                		iw = (k0 + i0) % 2 + 2 * ((k1 + i1) % 2)

                		fp2_add(&w[iw],&w[iw],&zeta_powers[e]);
					}
				}
				not_all_zero=0;
				for(int i=0;i<4;i++){
					not_all_zero=not_all_zero|fp2_is_zero(&w[i]);
				}
				ir0=((not_all_zero-1)&i0)|(~(not_all_zero-1)&ir0);
				ir1=((not_all_zero-1)&i1)|(~(not_all_zero-1)&ir1);
				for(int i=0;i<4;i++){
					fp2_cswap(&v[i], &w[i], not_all_zero);
				}
			}
		}
	}
} 


void 
basis_to_theta_change_dim2(fp2_t **N, const digit_t **M, const fp2_t *zeta, unsigned int nwords){

}