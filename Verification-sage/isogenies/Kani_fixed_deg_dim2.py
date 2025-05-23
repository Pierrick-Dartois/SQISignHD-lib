from sage.all import *

from utilities.norm_equation import solve_norm_equation
from utilities.supersingular import torsion_basis_2f_frobenius
from utilities.discrete_log import weil_pairing_pari
from basis_change.canonical_basis_dim1 import make_canonical
from basis_change.kani_base_change import fixed_deg_gluing_matrix, fixed_deg_splitting_matrix
from basis_change.base_change_dim4 import base_change_theta_dim4
from theta_structures.Tuple_point import TuplePoint
from theta_structures.Theta_dim1 import ThetaStructureDim1, ThetaPointDim1
from theta_structures.Theta_dim4 import ProductThetaStructureDim1To4, ThetaPointDim4
from theta_structures.theta_helpers_dim4 import product_to_theta_points_dim4_dim2
from isogenies.gluing_isogeny_dim4 import GluingIsogenyDim4
from isogenies.isogeny_chain_dim4 import IsogenyChainDim4
from isogenies.Kani_gluing_isogeny_chain_dim4 import KaniFixedDegDim2Gluing
from time import time

from basis_change.base_change_dim4 import is_symplectic_matrix_dim4

class GluingFixedDegDim2(GluingIsogenyDim4):
	def __init__(self,K_8,Theta_product,Theta_gluing,N_gluing):
		self.Theta_product = Theta_product
		self.Theta_gluing = Theta_gluing
		self.N_gluing = N_gluing

		# Theta conversion and kernel translates

		K_8_and_trans = K_8 + [K_8[0]+K_8[1], K_8[0]+K_8[2], K_8[0]+K_8[3], K_8[1]+K_8[2], K_8[1]+K_8[3], K_8[2]+K_8[3]]

		Th_K_8_and_trans = [self.Theta_product(T) for T in K_8_and_trans]
		Th_K_8_and_trans = [self.Theta_gluing.base_change_coords(self.N_gluing,T) for T in Th_K_8_and_trans]

		# Evaluation translates
		self.L_trans = [2*T for T in K_8_and_trans]+[2*(K_8_and_trans[4]+K_8_and_trans[2]),
		2*(K_8_and_trans[4]+K_8_and_trans[3]),2*(K_8_and_trans[5]+K_8_and_trans[3]),2*(K_8_and_trans[7]+K_8_and_trans[3]),
		2*(K_8_and_trans[4]+K_8_and_trans[9])]
		self.L_trans_ind = [1,2,4,8,3,5,9,6,10,12,7,11,13,14,15] # Contains the weight of (1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,1,0,0),
		#(1,0,1,0),(1,0,0,1),(0,1,1,0),(0,1,0,1),(0,0,1,1),(1,1,1,0),(1,1,0,1),(1,0,1,1),(0,1,1,1) and (1,1,1,1) respectively

		# Gluing isogeny construction
		L_multind=[(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(1,1,0,0),(1,0,1,0),(1,0,0,1),(0,1,1,0),(0,1,0,1),(0,0,1,1)]
		GluingIsogenyDim4.__init__(self,self.Theta_gluing,Th_K_8_and_trans,L_multind)

	def evaluate(self,P):
		if not isinstance(P, TuplePoint):
			raise TypeError("GluingFixedDegDim2 isogeny expects as input a TuplePoint on the domain product E^4")

		# Translating P
		L_P_trans = [P+T for T in self.L_trans]

		# Theta conversion
		Th_P = self.Theta_product(P)
		Th_L_P_trans = [self.Theta_product(T) for T in L_P_trans]

		Th_P=self.Theta_gluing.base_change_coords(self.N_gluing,Th_P)
		Th_L_P_trans = [self.Theta_gluing.base_change_coords(self.N_gluing,T) for T in Th_L_P_trans]

		return self.special_image(Th_P,Th_L_P_trans,self.L_trans_ind)

	def __call__(self,P):
		return self.evaluate(P)

class KaniFixedDegDim2(IsogenyChainDim4):
	r"""
	Computes a u-isogeny in dimension 2 
	Phi_u: E^2 ---> Au,
	where u is odd and close to sqrt(p).

	INPUT:
	- u: odd integer close to sqrt(p).
	- E: supersingular elliptic curve defined over F_{p^2}.
	- basis: list or tuple of points P, Q forming a
	basis of E[2^{e-1}] such that pi(P)=P and pi(Q)=-Q,
	where p=c*2**e-1 and pi is the p-th power Frobenius of E
	(optional, default = None).
	- strategy: optimal strategy to compute the isogeny chain
	(optional, default = None).

	OUTPUT: An efficient representation of Phi_u.
	"""

	def __init__(self,u,E,basis=None,strategy=None):
		
		Fp2=E.base_field()
		p=Fp2.characteristic()
		self.E=E

		print("Representing u in the quaternions.\n")
		t1=time()
		margin=100
		res=solve_norm_equation(u,p,margin)
		while res is None:
			print("Representing u failed: increasing margin.\n")
			margin*=2
			res=solve_norm_equation(u,p,margin)

		a,b,c,d,f=res
		t2=time()
		print("Representing u done in {} s.\n".format(t2-t1))


		if basis is None:
			P, Q = torsion_basis_2f_frobenius(E,f+2)
		else:
			e = valuation(p+1,2)
			P = basis[0]*2**(e-f-3)
			Q = basis[1]*2**(e-f-3)

		## Gluing isogeny

		# Number of steps before 4-dimensional gluing
		m = min(valuation(a+b,2),valuation(a-b,2))

		# Entry points
		P_mp3 = 2**(f-m-1)*P
		Q_mp3 = 2**(f-m-1)*Q

		# Instanciating gluing
		self.gluing_isogeny_chain=KaniFixedDegDim2Gluing(P_mp3,Q_mp3,a,b,c,d,u,f,m)

		## Isogeny chain

		# Kernel
		B_Kpp = [TuplePoint(u*P,E(0),(a+b)*P,(c+d)*P),
		TuplePoint(E(0),u*P,(d-c)*P,(a-b)*P),
		TuplePoint((u-2**f)*Q,E(0),(a-b)*Q,(c-d)*Q),
		TuplePoint(E(0),(u-2**f)*Q,(-c-d)*Q,(a+b)*Q)]

		# Instantiate isogeny chain
		IsogenyChainDim4.__init__(self, B_Kpp, self.gluing_isogeny_chain, f, m, splitting=True, strategy=strategy)

		## Splitting
		M_splitting = fixed_deg_splitting_matrix(u)
		self.N_splitting = base_change_theta_dim4(M_splitting,e4)
		self.codomain_product = self._isogenies[-1]._codomain.base_change_struct(self.N_splitting)

		# Compute theta structure of Au
		Au_null_point,_ = product_to_theta_points_dim4_dim2(self.codomain_product.null_point().coords())
		self.Au=ThetaStructureDim2(Au_null_point)

	def evaluate(self,P):
		FP0 = self.evaluate_isogeny(TuplePoint(P[0],self.E(0),self.E(0),P[1]))
		FP0 = self.codomain_product.base_change_coords(self.N_splitting,FP0)

		# Separate coordinates from the product.
		PhiuP, _ = product_to_theta_points_dim4_dim2(FP0.coords())
		PhiuP = self.Au(PhiuP)

		return PhiuP

	def __call__(self,P):
		return self.evaluate(P)

