from sage.all import *

import itertools

def complete_symplectic_matrix_dim2(C, D, n=4):
	Zn = Integers(n)

	# Compute first row
	F = D.stack(-C).transpose()
	xi = F.solve_right(vector(Zn, [1, 0]))

	# Rearrange TODO: why?
	v = vector(Zn, [xi[2], xi[3], -xi[0], -xi[1]])

	# Compute second row
	F2 = F.stack(v)
	yi = F2.solve_right(vector(Zn, [0, 1, 0]))
	M = Matrix(Zn,
		[
			[xi[0], yi[0], C[0, 0], C[0, 1]],
			[xi[1], yi[1], C[1, 0], C[1, 1]],
			[xi[2], yi[2], D[0, 0], D[0, 1]],
			[xi[3], yi[3], D[1, 0], D[1, 1]],
		],
	)

	return M

def block_decomposition(M):
	"""
	Given a 4x4 matrix, return the four 2x2 matrices as blocks
	"""
	A = M.matrix_from_rows_and_columns([0, 1], [0, 1])
	B = M.matrix_from_rows_and_columns([2, 3], [0, 1])
	C = M.matrix_from_rows_and_columns([0, 1], [2, 3])
	D = M.matrix_from_rows_and_columns([2, 3], [2, 3])
	return A, B, C, D

def is_symplectic_matrix_dim2(M):
	A, B, C, D = block_decomposition(M)
	if B.transpose() * A != A.transpose() * B:
		return False
	if C.transpose() * D != D.transpose() * C:
		return False
	if A.transpose() * D - B.transpose() * C != identity_matrix(2):
		return False
	return True

def base_change_theta_dim2(M, zeta):
    r"""
    Computes the matrix N (in row convention) of the new level 2
    theta coordinates after a symplectic base change of A[4] given
    by M\in Sp_4(Z/4Z). We have:
    (theta'_i)=N*(theta_i)
    where the theta'_i are the new theta coordinates and theta_i,
    the theta coordinates induced by the product theta structure.
    N depends on the fourth root of unity zeta=e_4(Si,Ti) (where
    (S0,S1,T0,T1) is a symplectic basis if A[4]).

    Inputs:
    - M: symplectic base change matrix.
    - zeta: a primitive 4-th root of unity induced by the Weil-pairings
    of the symplectic basis of A[4].

    Output: Matrix N of base change of theta-coordinates.
    """
    # Split 4x4 matrix into 2x2 blocks
    A, B, C, D = block_decomposition(M)

    # Initialise N to 4x4 zero matrix
    N = [[0 for _ in range(4)] for _ in range(4)]

    def choose_non_vanishing_index(C, D, zeta):
        """
        Choice of reference non-vanishing index (ir0,ir1)
        """
        for ir0, ir1 in itertools.product([0, 1], repeat=2):
            L = [0, 0, 0, 0]
            for j0, j1 in itertools.product([0, 1], repeat=2):
                k0 = C[0, 0] * j0 + C[0, 1] * j1
                k1 = C[1, 0] * j0 + C[1, 1] * j1

                l0 = D[0, 0] * j0 + D[0, 1] * j1
                l1 = D[1, 0] * j0 + D[1, 1] * j1

                e = -(k0 + 2 * ir0) * l0 - (k1 + 2 * ir1) * l1
                L[ZZ(k0 + ir0) % 2 + 2 * (ZZ(k1 + ir1) % 2)] += zeta ** (ZZ(e))

            # Search if any L value in L is not zero
            if any([x != 0 for x in L]):
                return ir0, ir1

    ir0, ir1 = choose_non_vanishing_index(C, D, zeta)

    for i0, i1, j0, j1 in itertools.product([0, 1], repeat=4):
        k0 = A[0, 0] * i0 + A[0, 1] * i1 + C[0, 0] * j0 + C[0, 1] * j1
        k1 = A[1, 0] * i0 + A[1, 1] * i1 + C[1, 0] * j0 + C[1, 1] * j1

        l0 = B[0, 0] * i0 + B[0, 1] * i1 + D[0, 0] * j0 + D[0, 1] * j1
        l1 = B[1, 0] * i0 + B[1, 1] * i1 + D[1, 0] * j0 + D[1, 1] * j1

        e = i0 * j0 + i1 * j1 - (k0 + 2 * ir0) * l0 - (k1 + 2 * ir1) * l1
        N[i0 + 2 * i1][ZZ(k0 + ir0) % 2 + 2 * (ZZ(k1 + ir1) % 2)] += zeta ** (ZZ(e))

    return Matrix(N)

def montgomery_to_theta_matrix_dim2(zero12,N=identity_matrix(4),return_null_point=False):
	r"""
	Computes the matrix that transforms Montgomery coordinates on E1*E2 into theta coordinates 
	with respect to a certain theta structure given by a base change matrix N.

	Input: 
	- zero12: theta-null point for the product theta structure on E1*E2[2]:
	zero12=[zero1[0]*zero2[0],zero1[1]*zero2[0],zero1[0]*zero2[1],zero1[1]*zero2[1]],
	where the zeroi are the theta-null points of Ei for i=1,2.
	- N: base change matrix from the product theta-structure.
	- return_null_point: True if the theta-null point obtained after applying N to zero12 is returned.

	Output:
	- Matrix for the change of coordinates:
	(X1*X2,Z1*X2,X1*Z2,Z1*Z2)-->(theta_00,theta_10,theta_01,theta_11).
	- If return_null_point, the theta-null point obtained after applying N to zero12 is returned.
	"""

	Fp2=zero12[0].parent()

	M=zero_matrix(Fp2,4,4)

	for i in range(4):
		for j in range(4):
			M[i,j]=N[i,j]*zero12[j]

	M2=[]
	for i in range(4):
		M2.append([M[i,0]+M[i,1]+M[i,2]+M[i,3],-M[i,0]-M[i,1]+M[i,2]+M[i,3],-M[i,0]+M[i,1]-M[i,2]+M[i,3],M[i,0]-M[i,1]-M[i,2]+M[i,3]])

	if return_null_point:
		null_point=[M[i,0]+M[i,1]+M[i,2]+M[i,3] for i in range(4)]
		return Matrix(M2), null_point
	else:
		return Matrix(M2)

def apply_base_change_theta_dim2(N,P):
	Q=[]
	for i in range(4):
		Q.append(0)
		for j in range(4):
			Q[i]+=N[i,j]*P[j]
	return Q

