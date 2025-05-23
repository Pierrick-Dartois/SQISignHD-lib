from sage.all import *
from utilities.supersingular import compute_linearly_independent_point
from theta_structures.montgomery_theta import (montgomery_point_to_theta_point, theta_point_to_montgomery_point, torsion_to_theta_null_point)


class ThetaStructureDim1:
	def __init__(self,E,P=None,Q=None):
		self.E=E

		a_inv=E.a_invariants()
		
		A =a_inv[1]
		if a_inv != (0,A,0,1,0):
			raise ValueError("The elliptic curve E is not in the Montgomery model.")

		if Q==None:
			y2=A-2
			y=y2.sqrt()
			Q=E([-1,y,1])
			P=compute_linearly_independent_point(E,Q,4)
		else:
			if Q[0]!=-1:
				raise ValueError("You should enter a canonical 4-torsion basis. Q[0] should be -1.")

		self.P=P
		self.Q=Q

		self._base_ring=E.base_ring()

		self._point=ThetaPointDim1
		self._null_point=self._point(self,torsion_to_theta_null_point(P))

	def null_point(self):
		"""
		"""
		return self._null_point

	def base_ring(self):
		"""
		"""
		return self._base_ring

	def zero(self):
		"""
		"""
		return self.null_point()

	def elliptic_curve(self):
		return self.E

	def torsion_basis(self):
		return (self.P,self.Q)

	def __call__(self,coords):
		r"""
		Input: either a tuple or list of 2 coordinates or an elliptic curve point.

		Output: the corresponding theta point for the self theta structure.
		"""
		if isinstance(coords,tuple):
			return self._point(self,coords)
		elif isinstance(coords,list):
			return self._point(self,coords)
		else:
			return self._point(self,montgomery_point_to_theta_point(self.null_point().coords(),coords))

	def __repr__(self):
		return f"Theta structure on {self.elliptic_curve()} with null point: {self.null_point()} induced by the 4-torsion basis {self.torsion_basis()}"

	def to_montgomery_point(self,P):
		return theta_point_to_montgomery_point(self.E,self.null_point().coords(),P.coords())


class ThetaPointDim1:
    def __init__(self, parent, coords):
        """
        """
        if not isinstance(parent, ThetaStructureDim1):
            raise ValueError("Entry parent should be a ThetaStructureDim1 object.")

        self._parent = parent
        self._coords = tuple(coords)


    def coords(self):
    	return self._coords

    def __repr__(self):
    	return f"Theta point with coordinates: {self.coords()}"





