from sage.all import *
import itertools

from theta_structures.theta_helpers_dim4 import multindex_to_index

def bloc_decomposition(M):
	I1=[0,1,2,3]
	I2=[4,5,6,7]

	A=M[I1,I1]
	B=M[I2,I1]
	C=M[I1,I2]
	D=M[I2,I2]

	return A,B,C,D

def mat_prod_vect(A,I):
	J=[]
	for i in range(4):
		J.append(0)
		for k in range(4):
			J[i]+=A[i,k]*I[k]
	return tuple(J)

def add_tuple(I,J):
	K=[]
	for k in range(4):
		K.append(I[k]+J[k])
	return tuple(K)

def scal_prod_tuple(I,J):
	s=0
	for k in range(4):
		s+=I[k]*J[k]
	return s

def red_mod_2(I):
	J=[]
	for x in I:
		J.append(ZZ(x)%2)
	return tuple(J)

def choose_non_vanishing_index(C,D,zeta):
	for I0 in itertools.product([0,1],repeat=4):
		L=[0 for k in range(16)]
		for J in itertools.product([0,1],repeat=4):
			CJ=mat_prod_vect(C,J)
			DJ=mat_prod_vect(D,J)
			e=-scal_prod_tuple(CJ,DJ)-2*scal_prod_tuple(I0,DJ)

			I0pDJ=add_tuple(I0,CJ)

			L[multindex_to_index(red_mod_2(I0pDJ))]+=zeta**(ZZ(e))
		for k in range(16):
			if L[k]!=0:
				return I0,L


def base_change_theta_dim4(M,zeta):
	
	Z4=Integers(4)
	A,B,C,D=bloc_decomposition(M.change_ring(Z4))

	I0,L0=choose_non_vanishing_index(C,D,zeta)


	N=[L0]+[[0 for j in range(16)] for i in range(15)]
	for I in itertools.product([0,1],repeat=4):
		if I!=(0,0,0,0):
			AI=mat_prod_vect(A,I)
			BI=mat_prod_vect(B,I)
			for J in itertools.product([0,1],repeat=4):
				CJ=mat_prod_vect(C,J)
				DJ=mat_prod_vect(D,J)

				AIpCJ=add_tuple(AI,CJ)
				BIpDJ=add_tuple(BI,DJ)

				e=scal_prod_tuple(I,J)-scal_prod_tuple(AIpCJ,BIpDJ)-2*scal_prod_tuple(I0,BIpDJ)
				N[multindex_to_index(I)][multindex_to_index(red_mod_2(add_tuple(AIpCJ,I0)))]+=zeta**(ZZ(e))

	Fp2=zeta.parent()
	return matrix(Fp2,N)

def apply_base_change_theta_dim4(N,P):
	Q=[]
	for i in range(16):
		Q.append(0)
		for j in range(16):
			Q[i]+=N[i,j]*P[j]
	return Q

def complete_symplectic_matrix_dim4(C,D,n=4):
	Zn=Integers(n)

	Col_I4=[matrix(Zn,[[1],[0],[0],[0]]),matrix(Zn,[[0],[1],[0],[0],[0]]),
	matrix(Zn,[[0],[0],[1],[0],[0],[0]]),matrix(Zn,[[0],[0],[0],[1],[0],[0],[0]])]

	L_DC=block_matrix([[D.transpose(),-C.transpose()]])
	Col_AB_i=L_DC.solve_right(Col_I4[0])
	A_t=Col_AB_i[[0,1,2,3],0].transpose()
	B_t=Col_AB_i[[4,5,6,7],0].transpose()

	for i in range(1,4):
		F=block_matrix(2,1,[L_DC,block_matrix(1,2,[B_t,-A_t])])
		Col_AB_i=F.solve_right(Col_I4[i])
		A_t=block_matrix(2,1,[A_t,Col_AB_i[[0,1,2,3],0].transpose()])
		B_t=block_matrix(2,1,[B_t,Col_AB_i[[4,5,6,7],0].transpose()])

	A=A_t.transpose()
	B=B_t.transpose()

	M=block_matrix([[A,C],[B,D]])

	return M

def random_symmetric_matrix(n,r):
	Zn=Integers(n)
	M=zero_matrix(Zn,r,r)
	for i in range(r):
		M[i,i]=randint(0,n-1)
	for i in range(r):
		for j in range(i+1,r):
			M[i,j]=randint(0,n-1)
			M[j,i]=M[i,j]
	return M


def random_symplectic_matrix(n):
	Zn=Integers(n)

	C=random_symmetric_matrix(n,4)
	Cp=random_symmetric_matrix(n,4)

	A=random_matrix(Zn,4,4)
	while not A.is_invertible():
		A=random_matrix(Zn,4,4)

	tA_inv=A.inverse().transpose()
	ACp=A*Cp
	return block_matrix([[ACp,A],[C*ACp-tA_inv,C*A]])

def is_symplectic_matrix_dim4(M):
	A,B,C,D=bloc_decomposition(M)
	if B.transpose()*A!=A.transpose()*B:
		return False
	if C.transpose()*D!=D.transpose()*C:
		return False
	if A.transpose()*D-B.transpose()*C!=identity_matrix(4):
		return False
	return True

