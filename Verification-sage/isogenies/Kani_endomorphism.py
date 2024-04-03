from sage.all import *

from basis_change.kani_base_change import (base_change_canonical_dim2, 
	gluing_base_change_matrix_dim2, 
	gluing_base_change_matrix_dim2_dim4, splitting_base_change_matrix_dim4,
	gluing_base_change_matrix_dim2_F1, gluing_base_change_matrix_dim2_F2,
	starting_two_symplectic_matrices, gluing_base_change_matrix_dim2_dim4_F1, 
	gluing_base_change_matrix_dim2_dim4_F2, kernel_basis)
from basis_change.canonical_basis_dim1 import make_canonical
from basis_change.base_change_dim2 import base_change_theta_dim2, montgomery_to_theta_matrix_dim2
from basis_change.base_change_dim4 import base_change_theta_dim4
from theta_structures.Theta_dim1 import ThetaStructureDim1, ThetaPointDim1
from theta_structures.Theta_dim2 import ThetaStructureDim2, ProductThetaStructureDim2, ThetaPointDim2
from theta_structures.Tuple_point import TuplePoint
from theta_structures.Theta_dim4 import ProductThetaStructureDim1To4, ProductThetaStructureDim2To4
from utilities.discrete_log import weil_pairing_pari
from isogenies_dim2.isogeny_dim2 import ThetaIsogenyDim2
from isogenies_dim2.gluing_isogeny_dim2 import GluingThetaIsogenyDim2
from isogenies.gluing_isogeny_dim4 import GluingIsogenyDim4
from isogenies.isogeny_dim4 import IsogenyDim4
from isogenies.Kani_gluing_isogeny_chain_dim4 import KaniGluingIsogenyChainDim4, KaniGluingIsogenyChainDim4Half
from isogenies.isogeny_chain_dim4 import IsogenyChainDim4

def form_kernel(x,y,P1,Q1,R2,S2,T1,U1):
	r""" Given the data outputted by generate_kernel, computes a basis B_Kpp of a maximal isotropic
	subgroup of E_1^2\times E_2^2[2**(e_A+2)] lying above the kernel of: 
	F:=[[\alpha_1,Diag(\hat{\sigma},\hat{\sigma})],
	    [-Diag(\sigma,\sigma),\tilde{\alpha_2}]]\in\End(E_1^2\times E_2^2)
	given by:
	K:=\{(\tilde{\alpha_1}(P),\Diag(\sigma,\sigma)(P)), P\in E_1^2[2**e_A]\}.
	(see generate_kernel description).

	Input:
	- x,y: integers defining
	\alpha_i:=[[x, y],
			   [-y, x]]\in\End(E_i^2).
	- P1,Q1: a basis of E1[2**(e_A+2)] inducing a canonical basis of E1[4].
	- R2, S2: the image of (P1, Q1) by a 3**e_B isogeny \sigma: E1-->E2, forming a basis of E2[2**(e_A+2)].
	- T1, U1: the canonical basis of E1[4] induced by (P1,Q1) (T1=2**e_A*P1 and U1=2**e_A*Q1).

	Output:
	- B_Kpp: basis B_Kpp of a maximal isotropic subgroup of E_1^2\times E_2^2[2**(e_A+2)] 
	lying above the kernel of F (4*B_Kpp=B_K).
	"""
	xP1=x*P1
	yP1=y*P1
	xQ1=x*Q1
	yQ1=y*Q1
	
	E2=R2.curve()
	OE2=E2(0)

	if x%2==1:
		a=-inverse_mod(x,4)
		B_Kpp=[TuplePoint(xP1+a*T1,yP1,R2,OE2),TuplePoint(xQ1,yQ1,S2,OE2),TuplePoint(-yP1,xP1,OE2,R2),TuplePoint(-yQ1,xQ1+a*U1,OE2,S2)]
	else:
		b=-inverse_mod(y,4)
		B_Kpp=[TuplePoint(xP1,yP1+b*T1,R2,OE2),TuplePoint(xQ1,yQ1,S2,OE2),TuplePoint(-yP1,xP1,OE2,R2),TuplePoint(-yQ1-b*U1,xQ1,OE2,S2)]

	return B_Kpp

class KaniEndo(IsogenyChainDim4):
	r"""Class representing a full isogeny chain E1^2*E2^2 --> E1^2*E2^2 derived from Kani's lemma 
	when given torsion point information and degrees.

	INPUT:
	- P1, Q1: Basis of E1[2**f].
	- R2, S2: image of (P1,Q1) via a q-isogeny \sigma: E1-->E2.
	- q, a1, a2, e: integers such that a1**2+a2**2+q=2**e.
	- f: integer such that f>=e+2 specifying the available 2**f-torsion in E1 and E2. 
	By default, f=e+2.

	OUTPUT: A chain of 2-isogenies in dimension 4 representing the endomorphism:
	F:=[[\alpha_1,Diag(\hat{\sigma},\hat{\sigma)],
	    [-Diag(\sigma,\sigma),\tilde{\alpha_2}]]\in\End(E_1^2\times E_2^2)
	with:
	\alpha_i:=[[a1, a2],
			   [-a2, a1]]\in\End(E_i^2).
	"""
	def __init__(self,P1,Q1,R2,S2,q,a1,a2,e,f=None,strategy=None):
		if q+a1**2+a2**2!=2**e:
			raise ValueError("Wrong parameters: q+a1^2+a2^2!=2^e")

		# Make sure that a1 and a2 are in the right order
		self.swap=False
		if a1%2==0:
			a1,a2=a2,a1
			self.swap=True

		if f==None:
			f=e+2
		elif f!=e+2:
			N=2**(f-e-2)
			P1=N*P1
			Q1=N*Q1
			R2=N*R2
			S2=N*S2
			f=e+2

		E1=P1.curve()
		E2=R2.curve()

		# Number of dimension 2 steps before dimension 4 gluing Am^2-->B 
		m=0
		a2_div=a2
		while a2_div%2==0:
			m+=1
			a2_div=a2_div//2

		P1_doubles,Q1_doubles,R2_doubles,S2_doubles=[P1],[Q1],[R2],[S2]
		for i in range(e):
			P1_doubles.append(2*P1_doubles[-1])
			Q1_doubles.append(2*Q1_doubles[-1])
			R2_doubles.append(2*R2_doubles[-1])
			S2_doubles.append(2*S2_doubles[-1])

		points_m=[P1_doubles[-m-2],Q1_doubles[-m-2],R2_doubles[-m-2],S2_doubles[-m-2]]
		points_4=[P1_doubles[-1],Q1_doubles[-1],R2_doubles[-1],S2_doubles[-1]]
		self.gluing_isogeny_chain=KaniGluingIsogenyChainDim4(points_m,points_4,a1,a2,q,m)

		B_Kpp=form_kernel(a1,a2,P1_doubles[0],Q1_doubles[0],R2_doubles[0],S2_doubles[0],P1_doubles[-1],Q1_doubles[-1])

		IsogenyChainDim4.__init__(self, B_Kpp, self.gluing_isogeny_chain, e, m, splitting=True, strategy=strategy)

		# Splitting
		mu=self.gluing_isogeny_chain.M_gluing_dim2[0,1]
		M_split=splitting_base_change_matrix_dim4(a1,a2,q,m,self.gluing_isogeny_chain.M_product_dim2,self.gluing_isogeny_chain.M_gluing_dim4[range(8),range(4)],mu)

		self.N_split=base_change_theta_dim4(M_split.inverse(),self.gluing_isogeny_chain.e4)#.inverse()

		Theta1=self.gluing_isogeny_chain.Theta1
		Theta2=self.gluing_isogeny_chain.Theta2
		self.codomain_product=ProductThetaStructureDim1To4(Theta1,Theta1,Theta2,Theta2)

	def evaluate(self,P):
		if self.swap:
			Q=TuplePoint(P[1],P[0],P[2],-P[3])
		else:
			Q=P
		FP=self.evaluate_isogeny(Q)
		FP=self.codomain_product.base_change_coords(self.N_split,FP)
		FP=self.codomain_product.to_tuple_point(FP)
		if self.swap:
			FP=TuplePoint(FP[0],-FP[1],FP[3],FP[2])
		return FP

	def __call__(self,P):
		return self.evaluate(P)


class KaniEndoHalf:
	r"""Class representing an isogeny chain E1^2*E2^2 --> E1^2*E2^2 derived from Kani's lemma 
	when given torsion point information sufficient to represent haf the chain and degrees.

	Since not enough torsion is given to compute the chain at once, the computation is divided
	into two as specified in https://eprint.iacr.org/2023/436, Section 4.3. Namely, we compute
	two isogeny chains F1: E1^2*E2^2-->C and \tilde{F2}: E1^2*E2^2-->C such that F=F2\circ F1.

	INPUT:
	- P1, Q1: Basis of E1[2**f].
	- R2, S2: image of (P1,Q1) via a q-isogeny \sigma: E1-->E2.
	- q, a1, a2, e: integers such that a1**2+a2**2+q=2**e.
	- f: integer such that f>=\ceil(e/2)+2 specifying the available 2**f-torsion in E1 and E2.

	OUTPUT: A chain of 2-isogenies in dimension 4 representing the endomorphism:
	F:=[[\alpha_1,Diag(\hat{\sigma},\hat{\sigma)],
	    [-Diag(\sigma,\sigma),\tilde{\alpha_2}]]\in\End(E_1^2\times E_2^2)
	with:
	\alpha_i:=[[a1, a2],
			   [-a2, a1]]\in\End(E_i^2).
	"""
	def __init__(self,P1,Q1,R2,S2,q,a1,a2,e,f,strategy1=None,strategy2=None):
		e1=ceil(e/2)
		e2=e-e1

		if q+a1**2+a2**2!=2**e:
			raise ValueError("Wrong parameters: q+a1^2+a2^2!=2^e")

		# Make sure that a1 and a2 are in the right order
		self.swap=False
		if a1%2==0:
			a1,a2=a2,a1
			self.swap=True

		# Number of dimension 2 steps before dimension 4 gluing Am^2-->B 
		m=0
		a2_div=a2
		while a2_div%2==0:
			m+=1
			a2_div=a2_div//2

		E1=P1.curve()
		E2=R2.curve()

		Fp2=E1.base_field()

		# Matrices of the symplectic basis B1 and B2 of E1^2*E2^2[2**f] adapted to ker(F1)
		# and ker(\tilde{F2}) respectively in the symplectic basis:
		# B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R2,0),(0,,0,0,R2)],
		# [(Q1,0,0,0),(0,Q1,0,0),(0,0,[1/q]*S2,0),(0,0,0,[1/q]*S2)]]
		M1, M2 = starting_two_symplectic_matrices(a1, a2, q, f)

		# Matrices of the symplectic basis induced by B1 and B2 on Am^2[4] and Apm^2[4]
		# (the m-th intermediary codomains of F1 and \tilde{F2} respectively) in the
		# symplectic basis of Am^2[4] and Apm^2[4] associated to the codomain theta-structures
		# after m isogeny computations in dimension 2.
		M_gluing_1=gluing_base_change_matrix_dim2_dim4_F1(a1,a2,q,m,M1)
		M_gluing_2=gluing_base_change_matrix_dim2_dim4_F2(a1,a2,q,m,M2)

		# Point doublings
		P1_doubles,Q1_doubles,R2_doubles,S2_doubles=[P1],[Q1],[R2],[S2]
		for i in range(f-2):
			P1_doubles.append(2*P1_doubles[-1])
			Q1_doubles.append(2*Q1_doubles[-1])
			R2_doubles.append(2*R2_doubles[-1])
			S2_doubles.append(2*S2_doubles[-1])


		# Base change matrix of the symplectic basis of E1*E2[4] induced by:
		# B1:=[(P1,0),(0,R2),(Q1,0),(0,[1/q]*S2)]
		# in the basis of E1*E2[4] induced by:
		# B0:=[(T1,0),(0,U1),(T2,0),(0,U2)]
		lamb=inverse_mod(q,4)

		_,_,T1,T2,MT=make_canonical(P1_doubles[-1],Q1_doubles[-1],4,preserve_pairing=True)
		_,_,U1,U2,MU=make_canonical(R2_doubles[-1],lamb*S2_doubles[-1],4,preserve_pairing=True)

		
		Z4=Integers(4)
		M0=matrix(Z4,[[MT[0,0],0,MT[1,0],0],
			[0,MU[0,0],0,MU[1,0]],
			[MT[0,1],0,MT[1,1],0],
			[0,MU[0,1],0,MU[1,1]]])

		# Theta structures in dimension 1 and 2
		Theta1=ThetaStructureDim1(E1,T1,T2)
		Theta2=ThetaStructureDim1(E2,U1,U2)

		Theta12=ProductThetaStructureDim2(Theta1,Theta2)

		e4=Fp2(weil_pairing_pari(T1,T2,4))

		# Gluing isogeny chains
		points_m=[P1_doubles[-m-2],Q1_doubles[-m-2],R2_doubles[-m-2],S2_doubles[-m-2]]# Points of order 2**(m+3)

		self.gluing_isogeny_chain1=KaniGluingIsogenyChainDim4Half(points_m, a1, a2, q, m, Theta12, M0, M1, M_gluing_1, e4, False)
		self.gluing_isogeny_chain2=KaniGluingIsogenyChainDim4Half(points_m, a1, a2, q, m, Theta12, M0, M2, M_gluing_2, e4, True)

		# Full isogeny chains
		lamb1=inverse_mod(q,2**(e1+2))
		lamb2=inverse_mod(q,2**(e2+2))
		B_Kpp1=kernel_basis(M1,e1,P1_doubles[f-e1-2],Q1_doubles[f-e1-2],R2_doubles[f-e1-2],lamb1*S2_doubles[f-e1-2])
		B_Kpp2=kernel_basis(M2,e2,P1_doubles[f-e2-2],Q1_doubles[f-e2-2],R2_doubles[f-e2-2],lamb2*S2_doubles[f-e2-2])

		self.F1=IsogenyChainDim4(B_Kpp1, self.gluing_isogeny_chain1, e1, m, splitting=False, strategy=strategy1)
		self.F2_dual=IsogenyChainDim4(B_Kpp2, self.gluing_isogeny_chain2, e2, m, splitting=False, strategy=strategy2)

		self.F2=self.F2_dual.dual()

	def evaluate(self,P):
		if self.swap:
			Q=TuplePoint(P[1],P[0],P[2],-P[3])
		else:
			Q=P
		FP=self.F1.evaluate_isogeny(Q)
		FP=self.F2.evaluate_isogeny(FP)
		if self.swap:
			FP=TuplePoint(FP[0],-FP[1],FP[3],FP[2])
		return FP

	def __call__(self,P):
		return self.evaluate(P)











		