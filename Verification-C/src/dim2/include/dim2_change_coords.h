#include <hd.h>
#include <fp2.h>
#include <fp2_matrix.h>
#include <mp.h>

/** @brief Computes the change of theta coordinates matrix N (in level 2),
	 * given a symplectic change of basis of the 4-torsion.
	 * 
	 * @param N Output: the change of coordinates matrix.
	 * @param M: a symplectic change of basis of the 4-torsion.
	 * @param zeta: the weil-paring of the symplectic basis.
	*/
void basis_to_theta_change_dim2(fp2_t **N, const digit_t **M, const fp2_t *zeta);

/** @brief Computes the matrix that transforms Montgomery coordinates on E1*E2 into theta coordinates 
 * with respect to a certain theta structure given by a base change matrix N.
 * 
 * @param out_matrix Output: a matrix that maps (x1*x2,x1*z2,z1*x2,z1*z2) to the theta coordinates
 * of P1=(x1:z1) and P2=(x2:z2) after change of theta coordinates from product theta coordinates via N.
 * @param out_theta_struct Output: the theta structure obtained after change of theta coordinates 
 * from product theta coordinates via N.
 * @param prod_null_point: thet theta null point for the product theta structure.
 * @param N: the change of theta coordinates matrix.  
 **/
void montgomery_to_theta_matrix_dim2(fp2_t **out_matrix, theta_point_t *out_null_point, 
	const theta_point_t *prod_null_point, const fp2_t **N);