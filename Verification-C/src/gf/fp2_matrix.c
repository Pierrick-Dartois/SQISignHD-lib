#include <fp2_matrix.h>

void 
mat_vec_prod(fp2_t *res, const fp2_t **M,const fp2_t *v, unsigned int n, unsigned int m){
	fp2_t t;
	for(int i=0;i<n;i++){
		fp2_set_zero(&res[i]);
		for(int j=0;j<m;j++){
			fp2_mul(&t,&M[i][j],&v[j]);
			fp2_add(&res[i],&res[i],&t);
		}
	}
}