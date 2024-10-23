#include <dim2_change_coords.h>

static inline void
choose_non_vanishing_index(int ir0, int ir1, fp2_t *v, fp2_t *zeta_powers, const digit_t **M, const fp2_t *zeta){
	/* Auxiliary function meant to chose a non vanishing index to define the coordinate theta'_0. 
	This is used in basis_to_theta_change_dim2. */
	fp2_t w[4];
	int k0, k1, l0, l1, e, iw;
	uint32_t not_all_zero;

	fp2_set_one(&zeta_powers[0]);
	fp2_copy(&zeta_powers[1],zeta);
	fp2_sqr(&zeta_powers[2],zeta);
	fp2_mul(&zeta_powers[3],&zeta_powers[2],zeta);

	for(int i0=0;i0<2;i0++){
		for(int i1=0;i1<2;i1++){
			w[0]=0;
			w[1]=0;
			w[2]=0;
			w[3]=0;
			for(int j0=0;j0<2;j0++){
				for(int j1=0;j1<2;j1++){
					k0 = M[0][2] * j0 + M[0][3] * j1;
               		k1 = M[1][2] * j0 + M[1][3] * j1;

               		l0 = M[2][2] * j0 + M[2][3] * j1;
               		l1 = M[3][2] * j0 + M[3][3] * j1;

               		e = (-(k0 + 2 * i0) * l0 - (k1 + 2 * i1) * l1) % 4;
               		iw = (k0 + i0) % 2 + 2 * ((k1 + i1) % 2);

               		fp2_add(&w[iw],&w[iw],&zeta_powers[e]);
				}
			}
			not_all_zero=0;
			for(int i=0;i<4;i++){
				not_all_zero=not_all_zero|(~fp2_is_zero(&w[i]));
			}

			ir0=((~not_all_zero)&i0)|(not_all_zero&ir0);
			ir1=((~not_all_zero)&i1)|(not_all_zero&ir1);
			for(int i=0;i<4;i++){
				fp2_select(&v[i], &v[i], &w[i], not_all_zero);
			}
		}
	}
}


void 
basis_to_theta_change_dim2(fp2_t **N, const digit_t **M, const fp2_t *zeta){
	/** @brief Computes the change of theta coordinates matrix N (in level 2),
	 * given a symplectic change of basis of the 4-torsion.
	 * 
	 * @param N Output: the change of coordinates matrix.
	 * @param M: a symplectic change of basis of the 4-torsion.
	 * @param zeta: the weil-paring of the symplectic basis.
	*/
	int ir0, ir1, k0, k1, l0, l1, e, jN;
	fp2_t v[4], zeta_powers[4];

	choose_non_vanishing_index(ir0, ir1, v, zeta_powers, M, zeta, nwords);

	for(int i=0;i<4;i++){
		if(i==0){
			for(int j=0;j<4;j++){
				fp2_copy(&N[i][j],&v[j]);
			}
		}
		else{
			for(int j=0;j<4;j++){
				fp2_set_zero(&N[i][j]);
			}
		}
	}

	for(int i0=0;i0<2;i0++){
		for(int i1=0;i1<2;i1++){
			if((i0!=0)||(i1!=0)){
				for(int j0=0;j0<2;j0++){
					for(int j1=0;j1<2;j1++){
						k0 = M[0][0] * i0 + M[0][1] * i1 + M[0][2] * j0 + M[0][3] * j1;
        				k1 = M[1][0] * i0 + M[1][1] * i1 + M[1][2] * j0 + M[1][3] * j1;

        				l0 = M[2][0] * i0 + M[2][1] * i1 + M[2][2] * j0 + M[2][3] * j1;
        				l1 = M[3][0] * i0 + M[3][1] * i1 + M[3][2] * j0 + M[3][3] * j1;

        				e = (i0 * j0 + i1 * j1 - (k0 + 2 * ir0) * l0 - (k1 + 2 * ir1) * l1)%4;
        				jN = (k0 + ir0) % 2 + 2 * ((k1 + ir1) % 2);
        				fp2_add(&N[i0 + 2 * i1][jN],&N[i0 + 2 * i1][jN],&zeta_powers[e]);
					}
				}
			}	
		}
	}
}


void montgomery_to_theta_matrix_dim2(fp2_t **out_matrix, theta_point_t *out_null_point, 
	const theta_point_t *prod_null_point, const fp2_t **N)
{
	fp2_t tmp_matrix[4][4];

	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			fp2_set_zero(&out_matrix[i][j]);
			fp2_mul(&tmp_matrix[i][j],&N[i][j],&prod_null_point[j]);
		}
	}

	for(int i=0;i<4;i++){
		// out_matrix[i][0]=tmp_matrix[i][0]+tmp_matrix[i][1]+tmp_matrix[i][2]+tmp_matrix[i][3]
		fp2_copy(&out_matrix[i][0],&tmp_matrix[i][0]);
		fp2_add(&out_matrix[i][0],&out_matrix[i][0],&tmp_matrix[i][1]);
		fp2_add(&out_matrix[i][0],&out_matrix[i][0],&tmp_matrix[i][2]);
		fp2_add(&out_matrix[i][0],&out_matrix[i][0],&tmp_matrix[i][3]);

		// out_matrix[i][1]=-tmp_matrix[i][0]-tmp_matrix[i][1]+tmp_matrix[i][2]+tmp_matrix[i][3]
		fp2_neg(&out_matrix[i][1],&tmp_matrix[i][0]);
		fp2_sub(&out_matrix[i][1],&out_matrix[i][1],&tmp_matrix[i][1]);
		fp2_add(&out_matrix[i][1],&out_matrix[i][1],&tmp_matrix[i][2]);
		fp2_add(&out_matrix[i][1],&out_matrix[i][1],&tmp_matrix[i][3]);

		// out_matrix[i][2]=-tmp_matrix[i][0]+tmp_matrix[i][1]-tmp_matrix[i][2]+tmp_matrix[i][3]
		fp2_neg(&out_matrix[i][2],&tmp_matrix[i][0]);
		fp2_add(&out_matrix[i][2],&out_matrix[i][2],&tmp_matrix[i][1]);
		fp2_sub(&out_matrix[i][2],&out_matrix[i][2],&tmp_matrix[i][2]);
		fp2_add(&out_matrix[i][2],&out_matrix[i][2],&tmp_matrix[i][3]);

		// out_matrix[i][3]=tmp_matrix[i][0]-tmp_matrix[i][1]-tmp_matrix[i][2]+tmp_matrix[i][3]
		fp2_copy(&out_matrix[i][3],&tmp_matrix[i][0]);
		fp2_sub(&out_matrix[i][3],&out_matrix[i][3],&tmp_matrix[i][1]);
		fp2_sub(&out_matrix[i][3],&out_matrix[i][3],&tmp_matrix[i][2]);
		fp2_add(&out_matrix[i][3],&out_matrix[i][3],&tmp_matrix[i][3]);
	}

	fp2_copy(&out_null_point->x,&out_matrix[0][0]);
	fp2_copy(&out_null_point->y,&out_matrix[1][0]);
	fp2_copy(&out_null_point->z,&out_matrix[2][0]);
	fp2_copy(&out_null_point->t,&out_matrix[3][0]);
}

void montgomery_to_theta(theta_point_t *out, const theta_couple_point_t *in, const fp2_t **N){
	fp2_t t[4], u[4];

	fp2_mul(&t[0],&in->P1.x,&in->P2.x);
	fp2_mul(&t[1],&in->P1.x,&in->P2.z);
	fp2_mul(&t[2],&in->P1.z,&in->P2.x);
	fp2_mul(&t[3],&in->P1.z,&in->P2.z);

	mat_vec_prod(u,N,t,4, 4);

	fp2_copy(&out->x,&u[0]);
	fp2_copy(&out->y,&u[1]);
	fp2_copy(&out->z,&u[2]);
	fp2_copy(&out->t,&u[3]);
}