from sage.all import *
from sage.structure.element import get_coercion_model
import itertools
import warnings

from theta_structures.theta_helpers_dim4 import (hadamard, 
	batch_inversion, 
	product_theta_point_dim4, 
	product_to_theta_points_dim4,
	product_theta_point_dim2_dim4,
	product_to_theta_points_dim4_dim2, 
	act_point,
	squared)
from theta_structures.Theta_dim1 import ThetaStructureDim1, ThetaPointDim1
from theta_structures.Tuple_point import TuplePoint
from basis_change.base_change_dim4 import apply_base_change_theta_dim4, random_symplectic_matrix, base_change_theta_dim4

cm = get_coercion_model()


class ThetaStructureDim4:
	def __init__(self,null_point,null_point_dual=None,inv_null_point_dual_sq=None):
		r"""
		INPUT:
		- null_point: theta-constants.
		- inv_null_point_dual_sq: inverse of the squares of dual theta-constants, if provided 
		(meant to prevent duplicate computation, since this data is already computed when the 
		codomain of an isogeny is computed).
		"""
		if not len(null_point) == 16:
			raise ValueError("Entry null_point should have 16 coordinates.")

		self._base_ring = cm.common_parent(*(c.parent() for c in null_point))
		self._point = ThetaPointDim4
		self._null_point = self._point(self, null_point)
		self._null_point_dual=null_point_dual
		self._inv_null_point=None
		self._inv_null_point_dual_sq=inv_null_point_dual_sq

	def null_point(self):
		"""
		"""
		return self._null_point

	def null_point_dual(self):
		if self._null_point_dual==None:
			self._null_point_dual=hadamard(self._null_point.coords())
		return self._null_point_dual

	def base_ring(self):
		"""
		"""
		return self._base_ring

	def zero(self):
		"""
		"""
		return self.null_point()

	def zero_dual(self):
		return self.null_point_dual()

	def __repr__(self):
		return f"Theta structure over {self.base_ring()} with null point: {self.null_point()}"

	def __call__(self,coords):
		return self._point(self,coords)

	def act_null(self,I,J):
		r"""
		Point of 2-torsion.

		INPUT:
		- I, J: two 4-tuples of indices in {0,1}.

		OUTPUT: the action of (I,\chi_J) on the theta null point given by:
		(I,\chi_J)*P=(\chi_J(I+K)^{-1}P[I+K])_K
		"""
		return self.null_point().act_point(I,J)

	def is_K2(self,B):
		r"""
		Given a symplectic decomposition A[2]=K_1\oplus K_2 canonically 
		induced by the theta-null point, determines if B is the canonical
		basis of K_2 given by act_nul(0,\delta_i)_{1\leq i\leq 4}.

		INPUT:
		- B: Basis of 4 points of 2-torsion.

		OUTPUT: Boolean True if and only if B is the canonical basis of K_2.
		"""
		I0=(0,0,0,0)
		if B[0]!=self.act_null(I0,(1,0,0,0)):
			return False
		if B[1]!=self.act_null(I0,(0,1,0,0)):
			return False
		if B[2]!=self.act_null(I0,(0,0,1,0)):
			return False
		if B[3]!=self.act_null(I0,(0,0,0,1)):
			return False
		return True

	def base_change_struct(self,N):
		null_coords=self.null_point().coords()
		new_null_coords=apply_base_change_theta_dim4(N,null_coords)
		return ThetaStructureDim4(new_null_coords)

	def base_change_coords(self,N,P):
		coords=P.coords()
		new_coords=apply_base_change_theta_dim4(N,coords)
		return self.__call__(new_coords)

	#@cached_method
	def _arithmetic_precomputation(self):
		r"""
		Initializes the precomputation containing the inverse of the theta-constants in standard
		and dual (Hadamard transformed) theta-coordinates. Assumes no theta-constant is zero.
		"""
		O=self.null_point()
		if all([O[k]!=0 for k in range(16)]):
			self._inv_null_point=batch_inversion(O.coords())
			if self._inv_null_point_dual_sq==None:
				U_chi_0_sq=hadamard(squared(O.coords()))
				if all([U_chi_0_sq[k]!=0 for k in range(16)]):	
					self._inv_null_point_dual_sq=batch_inversion(U_chi_0_sq)
					self._arith_base_change=False
				else:
					self._arith_base_change=True
					self._arithmetic_base_change()
					print("Zero dual theta constants.\nRandom symplectic base change is being used for duplication.\nDoublings are more costly than expected.")
		else:
			self._arith_base_change=True
			self._arithmetic_base_change()
			print("Zero theta constants.\nRandom symplectic base change is being used for duplication.\nDoublings are more costly than expected.")

	def _arithmetic_base_change(self,max_iter=50):
		i=self._base_ring.gen()
		count=0
		O=self.null_point()
		while count<max_iter:
			count+=1
			M=random_symplectic_matrix(4)
			N=base_change_theta_dim4(M,i)

			NO=apply_base_change_theta_dim4(N,O)
			NU_chi_0_sq=hadamard(squared(NO))

			if all([NU_chi_0_sq[k]!=0 for k in range(16)]) and all([NO[k]!=0 for k in range(16)]):
				self._arith_base_change_matrix=N
				self._arith_base_change_matrix_inv=N.inverse()
				self._inv_null_point=batch_inversion(NO)
				self._inv_null_point_dual_sq=batch_inversion(NU_chi_0_sq)
				break

	def has_suitable_doubling(self):
		O=self.null_point()
		UO=hadamard(O.coords())
		if all([O[k]!=0 for k in range(16)]) and all([UO[k]!=0 for k in range(16)]):
			return True
		else:
			return False

	def hadamard(self):
		return ThetaStructureDim4(self.null_point_dual())

	def hadamard_change_coords(self,P):
		new_coords=hadamard(P)
		return self.__call__(new_coords)


class ProductThetaStructureDim1To4(ThetaStructureDim4):
	def __init__(self,*args):
		r"""Defines the product theta structure at level 2 of 4 elliptic curves.

		Input: Either
		- 4 theta structures of dimension 1: T0, T1, T2, T3;
		- 4 elliptic curves: E0, E1, E2, E3.
		- 4 elliptic curves E0, E1, E2, E3 and their respective canonical 4-torsion basis B0, B1, B2, B3.
		"""
		if len(args)==4:
			theta_structures=list(args)
			for k in range(4):
				if not isinstance(theta_structures[k],ThetaStructureDim1):
					try:
						theta_structures[k]=ThetaStructureDim1(theta_structures[k])
					except:
						pass
		elif len(args)==8:
			theta_structures=[ThetaStructureDim1(args[k],args[4+k][0],args[4+k][1]) for k in range(4)]
		else:
			raise ValueError("4 or 8 arguments expected but {} were given.\nYou should enter a list of 4 elliptic curves or ThetaStructureDim1\nor a list of 4 elliptic curves with a 4-torsion basis for each of them.".format(len(args)))

		self._theta_structures=theta_structures

		null_point=product_theta_point_dim4(theta_structures[0].zero().coords(),theta_structures[1].zero().coords(),
			theta_structures[2].zero().coords(),theta_structures[3].zero().coords())

		ThetaStructureDim4.__init__(self,null_point)

	def product_theta_point(self,theta_points):
		return self._point(self,product_theta_point_dim4(theta_points[0].coords(),theta_points[1].coords(),
			theta_points[2].coords(),theta_points[3].coords()))

	def __call__(self,point):
		if isinstance(point,TuplePoint):
			theta_points=[]
			theta_structures=self._theta_structures
			for i in range(4):
				theta_points.append(theta_structures[i](point[i]))
			return self.product_theta_point(theta_points)
		else:
			return self._point(self,point)

	def to_theta_points(self,P):
		theta_coords=product_to_theta_points_dim4(P)
		theta_points=[self._theta_structures[i](theta_coords[i]) for i in range(4)]
		return theta_points

	def to_tuple_point(self,P):
		theta_points=self.to_theta_points(P)
		montgomery_points=[self._theta_structures[i].to_montgomery_point(theta_points[i]) for i in range(4)]
		return TuplePoint(montgomery_points)

class ProductThetaStructureDim2To4(ThetaStructureDim4):
	def __init__(self,theta1,theta2):
		self._theta_structures=(theta1,theta2)

		null_point=product_theta_point_dim2_dim4(theta1.zero().coords(),theta2.zero().coords())

		ThetaStructureDim4.__init__(self,null_point)

	def product_theta_point(self,P1,P2):
		return self._point(self,product_theta_point_dim2_dim4(P1.coords(),P2.coords()))

	def __call__(self,point):
		return self._point(self,point)

	def to_theta_points(self,P):
		theta_coords=product_to_theta_points_dim4_dim2(P)
		theta_points=[self._theta_structures[i](theta_coords[i]) for i in range(2)]
		return theta_points

class ThetaPointDim4:
    def __init__(self, parent, coords):
        """
        """
        if not isinstance(parent, ThetaStructureDim4):
            raise ValueError("Entry parent should be a ThetaStructureDim4 object.")

        self._parent = parent
        self._coords = tuple(coords)

    def parent(self):
        """
        """
        return self._parent
    
    def theta(self):
        """
        """
        return self.parent()
    
    def coords(self):
        """
        """
        return self._coords
    
    def is_zero(self):
        """
        """
        return self == self.parent().zero()

    def __eq__(self, other):
    	P=self.coords()
    	Q=other.coords()
    	
    	k0=0
    	while k0<15 and P[k0]==0:
    		k0+=1

    	for l in range(16):
    		if P[l]*Q[k0]!=Q[l]*P[k0]:
    			return False
    	return True

    def __repr__(self):
    	return f"Theta point with coordinates: {self.coords()}"

    def __getitem__(self,i):
    	return self._coords[i]

    def scale(self,lamb):
    	if lamb==0:
    		raise ValueError("Entry lamb should be non-zero.")

    	P=self.coords()
    	return self._parent([lamb*x for x in P])

    def act_point(self,I,J):
    	r"""
		Translation by a point of 2-torsion.

		Input:
		- I, J: two 4-tuples of indices in {0,1}.

		Output: the action of (I,\chi_J) on P given by:
		(I,\chi_J)*P=(\chi_J(I+K)^{-1}P[I+K])_K
		"""
    	return self._parent(act_point(self._coords,I,J))

    def double(self):
    	## This formula is projective.
    	## Works only when theta constants are non-zero.
    	P=self.coords()
    	if self.parent()._inv_null_point==None or self.parent()._inv_null_point_dual_sq==None:
    		self.parent()._arithmetic_precomputation()
    	inv_O,inv_U_chi_0_sq=self.parent()._inv_null_point,self.parent()._inv_null_point_dual_sq

    	if self.parent()._arith_base_change:
    		P=apply_base_change_theta_dim4(self.parent()._arith_base_change_matrix,P)

    	U_chi_P=squared(hadamard(squared(P)))
    	for chi in range(16):
    		U_chi_P[chi]*=inv_U_chi_0_sq[chi]

    	theta_2P = list(hadamard(U_chi_P))
    	for i in range(16):
    		theta_2P[i] *= inv_O[i]

    	if self.parent()._arith_base_change:
    		theta_2P=apply_base_change_theta_dim4(self.parent()._arith_base_change_matrix_inv,theta_2P)

    	return self._parent(theta_2P)

    def double_iter(self,n):
    	## Computes 2**n*self
    	Q=self
    	for i in range(n):
    		Q=Q.double()
    	return Q