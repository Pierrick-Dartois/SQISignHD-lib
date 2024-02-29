from sage.all import *
from utilities.discrete_log import weil_pairing_pari
from basis_change.canonical_basis_dim1 import make_canonical
from basis_change.kani_base_change import  (gluing_base_change_matrix_dim2, gluing_base_change_matrix_dim2_dim4,
	gluing_base_change_matrix_dim2_F1, gluing_base_change_matrix_dim2_F2, kernel_basis)
from basis_change.base_change_dim2 import base_change_theta_dim2, montgomery_to_theta_matrix_dim2
from basis_change.base_change_dim4 import base_change_theta_dim4
from theta_structures.Theta_dim1 import ThetaStructureDim1, ThetaPointDim1
from theta_structures.Theta_dim2 import ThetaStructureDim2, ProductThetaStructureDim2, ThetaPointDim2
from theta_structures.Tuple_point import TuplePoint
from theta_structures.Theta_dim4 import ProductThetaStructureDim2To4, ThetaPointDim4
from isogenies_dim2.isogeny_chain_dim2 import IsogenyChainDim2
from isogenies.gluing_isogeny_dim4 import GluingIsogenyDim4

class KaniGluingIsogenyChainDim4:
	def __init__(self,points_m,points_4,a1,a2,q,m,strategy_dim2=None):
		r"""

		INPUT: 
		- points_m: list of 4 points P1_m, Q1_m, R2_m, S2_m of order 2**(m+3)
		such that (P1_m,Q1_m) generates E1[2**(m+3)] and (R2_m,S2_m) is
		its image by sigma: E1 --> E2.
		- points_4: list of 4 points P1_4, Q1_4, R2_4, S2_4 of order 4 obtained by
		multiplying P1_m, Q1_m, R2_m, S2_m by 2**(m+1).
		- a1, a2, q: integers such that a1**2+a2**2+q=2**e.
		- m: 2-adic valuation of a2.

		OUTPUT: Composition of the m+1 first isogenies in the isogeny chained
		E1^2*E1^2 --> E1^2*E2^2 parametrized by a1, a2 and sigma via Kani's lemma.
		"""

		P1_m, Q1_m, R2_m, S2_m = points_m
		P1_4, Q1_4, R2_4, S2_4 = points_4

		E1=P1_m.curve()
		E2=R2_m.curve()

		Fp2=E1.base_field()

		lamb=inverse_mod(q,4)

		_,_,T1,T2,MT=make_canonical(P1_4,Q1_4,4,preserve_pairing=True)
		_,_,U1,U2,MU=make_canonical(R2_4,lamb*S2_4,4,preserve_pairing=True)

		Z4=Integers(4)
		M0=matrix(Z4,[[MT[0,0],0,MT[1,0],0],
			[0,MU[0,0],0,MU[1,0]],
			[MT[0,1],0,MT[1,1],0],
			[0,MU[0,1],0,MU[1,1]]])

		self.M_product_dim2=M0

		# Theta structures in dimension 1 and 2
		Theta1=ThetaStructureDim1(E1,T1,T2)
		Theta2=ThetaStructureDim1(E2,U1,U2)

		Theta12=ProductThetaStructureDim2(Theta1,Theta2)

		self.Theta1=Theta1
		self.Theta2=Theta2
		self.Theta12=Theta12

		# Gluing base change in dimension 2
		M1=gluing_base_change_matrix_dim2(a1,a2,q)
		M10=M0*M1

		self.M_gluing_dim2=M1

		e4=Fp2(weil_pairing_pari(T1,T2,4))

		self.e4=e4

		N_dim2=base_change_theta_dim2(M10,e4)
		#N_dim2=montgomery_to_theta_matrix_dim2(Theta12.zero().coords(),N1)

		# Gluing base change in dimension 4
		mua2=-M1[3,1]
		M2=gluing_base_change_matrix_dim2_dim4(a1,a2,m,mua2)

		self.M_gluing_dim4=M2

		self.N_dim4=base_change_theta_dim4(M2,e4)

		# Kernel of the 2**m-isogeny chain in dimension 2
		a1_red=a1%(2**(m+2))
		a2_red=a2%(2**(m+2))
		B_K_dim2=[TuplePoint(2*a1_red*P1_m-2*a2_red*Q1_m,2*R2_m),TuplePoint(2*a1_red*Q1_m+2*a2_red*P1_m,2*S2_m)]

		# Computation of the 2**m-isogeny chain in dimension 2
		self._isogenies_dim2=IsogenyChainDim2(B_K_dim2,Theta12,N_dim2,m,strategy_dim2)

		# Kernel of the (m+1)-th isogeny in dimension 4 f_{m+1}: A_m^2 --> B (gluing isogeny)
		a1_red=a1%(2**(m+3))
		a2_red=a2%(2**(m+3))
		
		a1P1_m=(a1_red)*P1_m
		a2P1_m=(a2_red)*P1_m
		a1Q1_m=(a1_red)*Q1_m
		a2Q1_m=(a2_red)*Q1_m

		OE2=E2(0)

		B_K_dim4=[TuplePoint(a1P1_m,a2P1_m,R2_m,OE2),TuplePoint(a1Q1_m,a2Q1_m,S2_m,OE2),
				TuplePoint(-a2P1_m,a1P1_m,OE2,R2_m),TuplePoint(-a2Q1_m,a1Q1_m,OE2,S2_m)]
		L_K_dim4=B_K_dim4+[B_K_dim4[0]+B_K_dim4[1]]

		L_K_dim4=[[self._isogenies_dim2(TuplePoint(L_K_dim4[k][0],L_K_dim4[k][2])),self._isogenies_dim2(TuplePoint(L_K_dim4[k][1],L_K_dim4[k][3]))] for k in range(5)]

		# Product Theta structure on A_m^2
		self.domain_product=ProductThetaStructureDim2To4(self._isogenies_dim2._codomain,self._isogenies_dim2._codomain)
		
		# Theta structure on A_m^2 after base change
		self.domain_base_change=self.domain_product.base_change_struct(self.N_dim4)
		
		# Converting the kernel to the Theta structure domain_base_change
		L_K_dim4=[self.domain_product.product_theta_point(L_K_dim4[k][0],L_K_dim4[k][1]) for k in range(5)]
		L_K_dim4=[self.domain_base_change.base_change_coords(self.N_dim4,L_K_dim4[k]) for k in range(5)]

		# Computing the gluing isogeny in dimension 4
		self._gluing_isogeny_dim4=GluingIsogenyDim4(self.domain_base_change,L_K_dim4,[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,1,0,0)])

		# Translates for the evaluation of the gluing isogeny in dimension 4
		self.L_trans=[2*B_K_dim4[k] for k in range(2)]
		self.L_trans_ind=[1,2] # Corresponds to multi indices (1,0,0,0) and (0,1,0,0)

		self._codomain=self._gluing_isogeny_dim4._codomain

	def evaluate(self,P):
		if not isinstance(P, TuplePoint):
			raise TypeError("KaniGluingIsogenyChainDim4 isogeny expects as input a TuplePoint on the domain product E1^2 x E2^2")

		# Translating P
		L_P_trans=[P+T for T in self.L_trans]
		
		# dim4 --> dim2 x dim2
		eval_P=[TuplePoint(P[0],P[2]),TuplePoint(P[1],P[3])]
		eval_L_P_trans=[[TuplePoint(Q[0],Q[2]),TuplePoint(Q[1],Q[3])] for Q in L_P_trans]

		# evaluating through the dimension 2 isogenies
		eval_P=[self._isogenies_dim2(eval_P[0]),self._isogenies_dim2(eval_P[1])]
		eval_L_P_trans=[[self._isogenies_dim2(Q[0]),self._isogenies_dim2(Q[1])] for Q in eval_L_P_trans]

		# Product Theta structure and base change
		eval_P=self.domain_product.product_theta_point(eval_P[0],eval_P[1])
		eval_P=self.domain_base_change.base_change_coords(self.N_dim4,eval_P)

		eval_L_P_trans=[self.domain_product.product_theta_point(Q[0],Q[1]) for Q in eval_L_P_trans]
		eval_L_P_trans=[self.domain_base_change.base_change_coords(self.N_dim4,Q) for Q in eval_L_P_trans]

		return self._gluing_isogeny_dim4.special_image(eval_P,eval_L_P_trans,self.L_trans_ind)

	def __call__(self,P):
		return self.evaluate(P)

class KaniGluingIsogenyChainDim4Half:
	def __init__(self, points_m, a1, a2, q, m, Theta12, M_product_dim2, M_start_dim4, M_gluing_dim4, e4, dual=False,strategy_dim2=None):#points_m,points_4,a1,a2,q,m,precomputed_data=None,dual=False,strategy_dim2=None):
		r"""

		INPUT: 
		- points_m: list of 4 points P1_m, Q1_m, R2_m, S2_m of order 2**(m+3)
		such that (P1_m,Q1_m) generates E1[2**(m+3)] and (R2_m,S2_m) is
		its image by sigma: E1 --> E2.
		- points_4: list of 4 points P1_4, Q1_4, R2_4, S2_4 of order 4 obtained by
		multiplying P1_m, Q1_m, R2_m, S2_m by 2**(m+1).
		- a1, a2, q: integers such that a1**2+a2**2+q=2**e.
		- m: 2-adic valuation of a2.

		OUTPUT: Composition of the m+1 first isogenies in the isogeny chained
		E1^2*E1^2 --> E1^2*E2^2 parametrized by a1, a2 and sigma via Kani's lemma.
		"""

		P1_m, Q1_m, R2_m, S2_m = points_m

		E1=P1_m.curve()
		E2=R2_m.curve()

		Fp2=E1.base_field()

		self.M_product_dim2 = M_product_dim2

		self.Theta12=Theta12

		self.e4=e4

		# Gluing base change in dimension 2
		if not dual:
			M1=gluing_base_change_matrix_dim2_F1(a1,a2,q)
		else:
			M1=gluing_base_change_matrix_dim2_F2(a1,a2,q)
		
		M10=M_product_dim2*M1

		self.M_gluing_dim2=M1

		self.e4=e4

		N_dim2=base_change_theta_dim2(M10,e4)
		#N_dim2=montgomery_to_theta_matrix_dim2(Theta12.zero().coords(),N1)

		# Gluing base change in dimension 4

		self.M_gluing_dim4 = M_gluing_dim4

		self.N_dim4 = base_change_theta_dim4(M_gluing_dim4, e4)

		# Kernel of the 2**m-isogeny chain in dimension 2
		a1_red=a1%(2**(m+2))
		a2_red=a2%(2**(m+2))
		if not dual:
			B_K_dim2=[TuplePoint(2*a1_red*P1_m-2*a2_red*Q1_m,2*R2_m),TuplePoint(2*a1_red*Q1_m+2*a2_red*P1_m,2*S2_m)]
		else:
			B_K_dim2=[TuplePoint(2*a1_red*P1_m+2*a2_red*Q1_m,-2*R2_m),TuplePoint(2*a1_red*Q1_m-2*a2_red*P1_m,-2*S2_m)]

		# Computation of the 2**m-isogeny chain in dimension 2
		self._isogenies_dim2=IsogenyChainDim2(B_K_dim2,Theta12,N_dim2,m,strategy_dim2)

		# Kernel of the (m+1)-th isogeny in dimension 4 f_{m+1}: A_m^2 --> B (gluing isogeny)
		lamb=inverse_mod(q,2**(m+3))
		B_K_dim4=kernel_basis(M_start_dim4,m+1,P1_m,Q1_m,R2_m,lamb*S2_m)
		L_K_dim4=B_K_dim4+[B_K_dim4[0]+B_K_dim4[1]]

		L_K_dim4=[[self._isogenies_dim2(TuplePoint(L_K_dim4[k][0],L_K_dim4[k][2])),self._isogenies_dim2(TuplePoint(L_K_dim4[k][1],L_K_dim4[k][3]))] for k in range(5)]

		# Product Theta structure on A_m^2
		self.domain_product=ProductThetaStructureDim2To4(self._isogenies_dim2._codomain,self._isogenies_dim2._codomain)
		
		# Theta structure on A_m^2 after base change
		self.domain_base_change=self.domain_product.base_change_struct(self.N_dim4)
		
		# Converting the kernel to the Theta structure domain_base_change
		L_K_dim4=[self.domain_product.product_theta_point(L_K_dim4[k][0],L_K_dim4[k][1]) for k in range(5)]
		L_K_dim4=[self.domain_base_change.base_change_coords(self.N_dim4,L_K_dim4[k]) for k in range(5)]

		# Computing the gluing isogeny in dimension 4
		self._gluing_isogeny_dim4=GluingIsogenyDim4(self.domain_base_change,L_K_dim4,[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,1,0,0)])

		# Translates for the evaluation of the gluing isogeny in dimension 4
		self.L_trans=[2*B_K_dim4[k] for k in range(2)]
		self.L_trans_ind=[1,2] # Corresponds to multi indices (1,0,0,0) and (0,1,0,0)

		self._codomain=self._gluing_isogeny_dim4._codomain

	def evaluate(self,P):
		if not isinstance(P, TuplePoint):
			raise TypeError("KaniGluingIsogenyChainDim4 isogeny expects as input a TuplePoint on the domain product E1^2 x E2^2")

		# Translating P
		L_P_trans=[P+T for T in self.L_trans]
		
		# dim4 --> dim2 x dim2
		eval_P=[TuplePoint(P[0],P[2]),TuplePoint(P[1],P[3])]
		eval_L_P_trans=[[TuplePoint(Q[0],Q[2]),TuplePoint(Q[1],Q[3])] for Q in L_P_trans]

		# evaluating through the dimension 2 isogenies
		eval_P=[self._isogenies_dim2(eval_P[0]),self._isogenies_dim2(eval_P[1])]
		eval_L_P_trans=[[self._isogenies_dim2(Q[0]),self._isogenies_dim2(Q[1])] for Q in eval_L_P_trans]

		# Product Theta structure and base change
		eval_P=self.domain_product.product_theta_point(eval_P[0],eval_P[1])
		eval_P=self.domain_base_change.base_change_coords(self.N_dim4,eval_P)

		eval_L_P_trans=[self.domain_product.product_theta_point(Q[0],Q[1]) for Q in eval_L_P_trans]
		eval_L_P_trans=[self.domain_base_change.base_change_coords(self.N_dim4,Q) for Q in eval_L_P_trans]

		return self._gluing_isogeny_dim4.special_image(eval_P,eval_L_P_trans,self.L_trans_ind)

	def __call__(self,P):
		return self.evaluate(P)

	def dual(self):
		domain = self._codomain.hadamard()
		codomain_base_change = self.domain_base_change
		codomain_product = self.domain_product
		N_dim4 = self.N_dim4.inverse()
		isogenies_dim2 = self._isogenies_dim2.dual()
		splitting_isogeny_dim4 = self._gluing_isogeny_dim4.dual()

		return KaniSplittingIsogenyChainDim4(domain, codomain_base_change, codomain_product, N_dim4, isogenies_dim2, splitting_isogeny_dim4)

class KaniSplittingIsogenyChainDim4:
	def __init__(self, domain, codomain_base_change, codomain_product, N_dim4, isogenies_dim2, splitting_isogeny_dim4):
		self._domain = domain
		self.codomain_base_change = codomain_base_change
		self.codomain_product = codomain_product
		self.N_dim4 = N_dim4
		self._isogenies_dim2 = isogenies_dim2
		self._splitting_isogeny_dim4 = splitting_isogeny_dim4

	def evaluate(self,P):
		if not isinstance(P, ThetaPointDim4):
			raise TypeError("KaniSplittingIsogenyChainDim4 isogeny expects as input a ThetaPointDim4")

		Q = self._splitting_isogeny_dim4(P)
		Q = self.codomain_product.base_change_coords(self.N_dim4, Q)
		Q1, Q2 = self.codomain_product.to_theta_points(Q)
		Q1, Q2 = self._isogenies_dim2._domain(Q1.hadamard()), self._isogenies_dim2._domain(Q2.hadamard()) 

		Q1 = self._isogenies_dim2(Q1)
		Q2 = self._isogenies_dim2(Q2)

		return TuplePoint(Q1[0],Q2[0],Q1[1],Q2[1])

	def __call__(self,P):
		return self.evaluate(P)