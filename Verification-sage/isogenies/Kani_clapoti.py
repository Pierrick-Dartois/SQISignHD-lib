from sage.all import *
import itertools

from basis_change.kani_base_change import clapoti_cob_splitting_matrix
from basis_change.canonical_basis_dim1 import make_canonical
from basis_change.base_change_dim2 import base_change_theta_dim2, montgomery_to_theta_matrix_dim2
from basis_change.base_change_dim4 import base_change_theta_dim4
from theta_structures.Theta_dim1 import ThetaStructureDim1, ThetaPointDim1
from theta_structures.Theta_dim2 import ThetaStructureDim2, ProductThetaStructureDim2, ThetaPointDim2
from theta_structures.Tuple_point import TuplePoint
from theta_structures.Theta_dim4 import ProductThetaStructureDim1To4, ProductThetaStructureDim2To4
from theta_structures.montgomery_theta import null_point_to_montgomery_coeff, theta_point_to_montgomery_point
from theta_structures.theta_helpers_dim4 import product_to_theta_points_dim4
from utilities.discrete_log import weil_pairing_pari
from utilities.supersingular import torsion_basis_to_Fp_rational_point
from isogenies_dim2.isogeny_dim2 import ThetaIsogenyDim2
from isogenies_dim2.gluing_isogeny_dim2 import GluingThetaIsogenyDim2
from isogenies.gluing_isogeny_dim4 import GluingIsogenyDim4
from isogenies.isogeny_dim4 import IsogenyDim4
from isogenies.Kani_gluing_isogeny_chain_dim4 import KaniClapotiGluing
from isogenies.isogeny_chain_dim4 import IsogenyChainDim4

class KaniClapotiIsog(IsogenyChainDim4):
	r"""Class representing the 4-dimensional isogeny obtained via Kani's lemma F: Eu^2*Ev^2 --> Ea^2*A 
	where Ea=[\mf{a}]*E is the result of the ideal class group action by \mf{a} when given relevant 
	constants and torsion point information.

	INPUT: 
	- Pu, Qu = phi_u(P, Q)\in Eu;
	- Pv, Qv = phi_v*phi_{ck}*\hat{\phi}_{bk}(P, Q)\in Ev;
	- gu, xu, yu, gv, xv, yv, Nbk, Nck, e: positive integers;
	where:
		* gu(xu^2+yu^2)Nbk+gv(xv^2+yv^2)Nck=2^e;
		* gcd(u*Nbk,v*Nck)=1 with u:=gu(xu^2+yu^2) and v:=gv(xv^2+yv^2);
		* xu and xv are odd and yu and yv are even;
		* \mf{b}=\mf{be}*\mf{bk} is a product of ideals of norms Nbe and Nbk respectively,
		where Nbe is a product of small Elkies primes;
		* \mf{c}=\mf{ce}*\mf{ck} is a product of ideals of norms Nce and Nck respectively,
		where Nbe is a product of small Elkies primes;
		* phi_{bk}: E --> E1 and phi_{ck}: E --> E2 are induced by the action of 
		ideals \mf{bk} and \mf{ck} respectively;
		* <P,Q>=E_1[2^{e+2}];
		* phi_u: E1 --> Eu and phi_v: E2 --> Ev are gu and gv-isogenies respectively.

	OUTPUT: F: Eu^2*Ev^2 --> Ea^2*A is the isogeny:
	
	F := [[Phi_{bp}*\tilde{Phi}_u, Phi_{cp}*\tilde{Phi}_v],
		 [-Psi, \tilde{Phi}]]
	
	obtained from the Kani isogeny diamond:

	A --------------------Phi------------------> Ev^2
	^                                             ^
	|                                             |
	|                                           Phi_v
	|                                             |
	Psi                                          E2^2
	|                                             ^
	|                                             |
	|                                     \tilde{Phi}_{ck}
	|                                             |
	Eu^2 --\tilde{Phi}_{u}--> E1^2 --Phi_{bk}--> Ea^2 

	where Phi_{bk}:=Diag(phi_{bk},phi_{bk}), Phi_{ck}:=Diag(phi_{ck},phi_{ck}),

	Phi_u := [[xu, -yu],
			[yu, xu]] * Diag(phi_u,phi_u)
	
	Phi_v := [[xv, -yv],
			[yv, xv]] * Diag(phi_v,phi_v)
	"""

	def __init__(self,points,integers,strategy=None):
		gu,xu,yu,gv,xv,yv,Nbk,Nck,e = integers
		Pu,Qu,Pv,Qv = points
		if gu*(xu**2+yu**2)*Nbk+gv*(xv**2+yv**2)*Nck!=2**e:
			raise ValueError("Wrong parameters: gu(xu^2+yu^2)Nbk + gv(xv^2+yv^2)Nck != 2^e")
		if gcd(ZZ(gu*(xu**2+yu**2)*Nbk),ZZ(gv*(xv**2+yv**2)*Nck))!=1:
			raise ValueError("Non coprime parameters: gcd(gu(xu^2+yu^2)Nbk, gv(xv^2+yv^2)Nck) != 1")
		if xu%2==0:
			xu,yu=yu,xu
		if xv%2==0:
			xv,yv=yv,xv

		self.Eu = Pu.curve()
		self.Ev = Pv.curve()
		Fp2 = parent(Pu[0])
		Fp = Fp2.base_ring()

		# Number of dimension 2 steps before dimension 4 gluing Am^2-->B 
		m=valuation(xv*yu-xu*yv,2)
		integers=[gu,xu,yu,gv,xv,yv,Nbk,Nck,e,m]

		points_mp3=[(2**(e-m-1))*P for P in points]
		points_mp2=[2*P for P in points_mp3]
		points_4=[(2**m)*P for P in points_mp2]

		self.Ru_Fp = torsion_basis_to_Fp_rational_point(self.Eu,points_4[0],points_4[1],4)

		self.gluing_isogeny_chain = KaniClapotiGluing(points_mp3,points_mp2,points_4,integers, coerce=Fp)

		xuNbk = xu*Nbk
		yuNbk = yu*Nbk
		two_ep2 = 2**(e+2)
		inv_Nbk = inverse_mod(Nbk,two_ep2)
		u = gu*(xu**2+yu**2)
		inv_u = inverse_mod(u,4)
		lambxu = ((1-2**e*inv_u*inv_Nbk)*xu)%two_ep2
		lambyu = ((1-2**e*inv_u*inv_Nbk)*yu)%two_ep2
		xv_Nbk = (xv*inv_Nbk)%two_ep2
		yv_Nbk = (yv*inv_Nbk)%two_ep2


		B_Kpp = [TuplePoint(xuNbk*Pu,yuNbk*Pu,xv*Pv,yv*Pv),
		TuplePoint(-yuNbk*Pu,xuNbk*Pu,-yv*Pv,xv*Pv),
		TuplePoint(lambxu*Qu,lambyu*Qu,xv_Nbk*Qv,yv_Nbk*Qv),
		TuplePoint(-lambyu*Qu,lambxu*Qu,-yv_Nbk*Qv,xv_Nbk*Qv)]

		IsogenyChainDim4.__init__(self, B_Kpp, self.gluing_isogeny_chain, e, m, splitting=True, strategy=strategy)

		# Splitting
		M_split = clapoti_cob_splitting_matrix(integers)

		self.N_split = base_change_theta_dim4(M_split,self.gluing_isogeny_chain.e4)

		self.codomain_product = self._isogenies[-1]._codomain.base_change_struct(self.N_split)

		# Extracting the group action image Ea=[\mathfrak{a}]*E from the codomain Ea^2*E'^2
		self.theta_null_Ea, self.theta_null_Ep, self.Ea, self.Ep = self.extract_montgomery_curve()



	def extract_montgomery_curve(self):

		# Computing the theta null point of Ea
		null_point=self.codomain_product.zero()
		Fp2=parent(null_point[0])
		Fp = Fp2.base_ring()
		for i3, i4 in itertools.product([0,1],repeat=2):
			if null_point[4*i3+8*i4]!=0:
				i30=i3
				i40=i4
				theta_Ea_0=Fp(null_point[4*i3+8*i4])
				theta_Ea_1=Fp(null_point[1+4*i3+8*i4])
				break
		for i1, i2 in itertools.product([0,1],repeat=2):
			if null_point[i1+2*i2]!=0:
				i10=i1
				i20=i2
				theta_Ep_0=Fp(null_point[i1+2*i2])
				theta_Ep_1=Fp(null_point[i1+2*i2+4])
				break

		# Sanity check: is the codomain of F a product of the form Ea^2*E'^2 ?
		theta_Ea=[Fp(theta_Ea_0),Fp(theta_Ea_1)]
		theta_Ep=[Fp(theta_Ep_0),Fp(theta_Ep_1)]

		theta_Ea2Ep2=[0 for i in range(16)]
		for i1,i2,i3,i4 in itertools.product([0,1],repeat=4):
			theta_Ea2Ep2[i1+2*i2+4*i3+8*i4]=theta_Ea[i1]*theta_Ea[i2]*theta_Ep[i3]*theta_Ep[i4]
		theta_Ea2Ep2=self.codomain_product(theta_Ea2Ep2)

		assert theta_Ea2Ep2.is_zero()

		A_Ep = null_point_to_montgomery_coeff(theta_Ep_0,theta_Ep_1)
		Ep = EllipticCurve([0,A_Ep,0,1,0])

		## ## Recovering Ea over Fp and not Fp2
		## self.find_Fp_rational_theta_struct_Ea()

		## theta_Ea = self.iso_Ea(theta_Ea)
		## A_Ea = null_point_to_montgomery_coeff(theta_Ea[0],theta_Ea[1])

		## # Sanity check : the curve Ea should be defined over Fp
		## # assert A_Ea[1] == 0

		## # Twisting Ea if necessary: if A_Ea+2 is not a square in Fp, then we take the twist (A_Ea --> -A_Ea)
		## p=self.Eu.base_field().characteristic()
		## self.twist = False
		## if (A_Ea+2)**((p-1)//2)==-1:
		## 	A_Ea = -A_Ea
		## 	self.twist = True

		## Ea = EllipticCurve([0,A_Ea,0,1,0])

		A = null_point_to_montgomery_coeff(theta_Ea_0, theta_Ea_1)
		Ab = null_point_to_montgomery_coeff(theta_Ea_0+theta_Ea_1, theta_Ea_0-theta_Ea_1)
		Acan = min([A, -A, Ab, -Ab])
		Acan = A
		if (Acan == A or Acan == -A):
			# 'Id' corresponds to the point on the twist
			self.iso_type = 'Id'
		else:
			# 'Hadamard' corresponds to the point on the curve
			self.iso_type = 'Hadamard'
		if ((self.iso_type == 'Hadamard' and not (Acan+2).is_square()) or (self.iso_type == 'Id' and (Acan+2).is_square())):
			Acan=-Acan
		if (Acan == A or Acan == Ab):
			self.twist = False
		else:
			self.twist = True
		Ea = EllipticCurve([0,Acan,0,1,0])

		# Find the dual null point
		return theta_Ea, theta_Ep, Ea, Ep

	def eval_rational_point_4_torsion(self):
		T = TuplePoint(self.Ru_Fp,self.Eu(0),self.Ev(0),self.Ev(0))

		FPu_4 = self.evaluate_isogeny(T)
		FPu_4=self.codomain_product.base_change_coords(self.N_split,FPu_4)
		FPu_4=product_to_theta_points_dim4(FPu_4)

		return FPu_4[0]

	def find_Fp_rational_theta_struct_Ea(self):
		Pa = self.eval_rational_point_4_torsion()

		HPa = (Pa[0]+Pa[1],Pa[0]-Pa[1])
		i = self.Eu.base_field().gen()
		self.i = i
		iHPa = (Pa[0]+i*Pa[1],Pa[0]-i*Pa[1])

		if Pa[0]==0 or Pa[1]==0:
			self.iso_type='Id'
		elif HPa[0]==0 or HPa[1]==0:
			self.iso_type='Hadamard'
		elif iHPa[0]==0 or iHPa[1]==0:
			self.iso_type='iHadamard'
		else:
			raise ValueError("A rational theta point should be mapped to (0:1) or (1:0) after change of theta coordinates on Ea.")

	def iso_Ea(self,P):
		# Change of theta coordinates to obtain Fp-rational theta coordinates on Ea

		if self.iso_type == 'Id':
			return P
		elif self.iso_type == 'Hadamard':
			return (P[0]+P[1],P[0]-P[1])
		else:
			return (P[0]+self.i*P[1],P[0]-self.i*P[1])


	def evaluate(self,P):
		FP=self.evaluate_isogeny(P)
		FP=self.codomain_product.base_change_coords(self.N_split,FP)
		
		FP=product_to_theta_points_dim4(FP)
		FP=TuplePoint([theta_point_to_montgomery_point(self.Ea,self.theta_null_Ea,self.iso_Ea(FP[0]),self.twist),
			theta_point_to_montgomery_point(self.Ea,self.theta_null_Ea,self.iso_Ea(FP[1]),self.twist),
			theta_point_to_montgomery_point(self.Ep,self.theta_null_Ep,FP[2]),
			theta_point_to_montgomery_point(self.Ep,self.theta_null_Ep,FP[3])])

		return FP

	def __call__(self,P):
		return self.evaluate(P)











