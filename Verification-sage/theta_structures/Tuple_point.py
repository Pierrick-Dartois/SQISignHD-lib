from sage.all import *
from utilities.discrete_log import weil_pairing_pari

class TuplePoint:
	def __init__(self,*args):
		if len(args)==1:
			self._points=list(args[0])
		else:
			self._points=list(args)

	def points(self):
		return self._points

	def parent_curves(self):
		return [x.curve() for x in self._points]

	def parent_curve(self,i):
		return self._points[i].curve()

	def n_points(self):
		return len(self._points)

	def is_zero(self):
		return all([self._points[i]==0 for i in range(self.n_points())])

	def __repr__(self):
		return str(self._points)

	def __getitem__(self,i):
		return self._points[i]

	def __setitem__(self,i,P):
		self._points[i]=P

	def __eq__(self,other):
		n_self=self.n_points()
		n_other=self.n_points()
		return n_self==n_other and all([self._points[i]==other._points[i] for i in range(n_self)])

	def __add__(self,other):
		n_self=self.n_points()
		n_other=self.n_points()
		
		if n_self!=n_other:
			raise ValueError("Cannot add TuplePoint of distinct lengths {} and {}.".format(n_self,n_other))

		points=[]
		for i in range(n_self):
			points.append(self._points[i]+other._points[i])
		return self.__class__(points)

	def __sub__(self,other):
		n_self=self.n_points()
		n_other=self.n_points()
		
		if n_self!=n_other:
			raise ValueError("Cannot substract TuplePoint of distinct lengths {} and {}.".format(n_self,n_other))

		points=[]
		for i in range(n_self):
			points.append(self._points[i]-other._points[i])
		return self.__class__(points)

	def __neg__(self):
		n_self=self.n_points()
		points=[]
		for i in range(n_self):
			points.append(-self._points[i])
		return self.__class__(points)

	def __mul__(self,m):
		n_self=self.n_points()
		points=[]
		for i in range(n_self):
			points.append(m*self._points[i])
		return self.__class__(points)

	def __rmul__(self,m):
		return self*m

	def double_iter(self,n):
		result=self
		for i in range(n):
			result=2*result
		return result

	def weil_pairing(self,other,n):
		n_self=self.n_points()
		n_other=self.n_points()
		
		if n_self!=n_other:
			raise ValueError("Cannot compute the Weil pairing of TuplePoint of distinct lengths {} and {}.".format(n_self,n_other))

		zeta=1
		for i in range(n_self):
			zeta*=weil_pairing_pari(self._points[i],other._points[i],n)

		return zeta







		