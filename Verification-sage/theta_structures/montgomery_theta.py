from sage.all import *

def torsion_to_theta_null_point(P):
	r=P[0]
	s=P[2]
	return (r+s,r-s)

def montgomery_point_to_theta_point(O,P):
	if P[0]==P[2]==0:
		return O
	else:
		a,b=O
		return (a*(P[0]-P[2]),b*(P[0]+P[2]))

def theta_point_to_montgomery_point(E,O,P,twist=False):
	a,b=O

	x=a*P[1]+b*P[0]
	z=a*P[1]-b*P[0]

	if twist:
		x=-x

	if z==0:
		return E(0)
	else:
		x=x/z

		a_inv=E.a_invariants()
		
		A =a_inv[1]
		if a_inv != (0,A,0,1,0):
			raise ValueError("The elliptic curve E is not in the Montgomery model.")

		y2=x*(x**2+A*x+1)
		if not is_square(y2):
			raise ValueError("The Montgomery point is not on the base field.")
		else:
			y=y2.sqrt()
			return E([x,y,1])

def lift_kummer_montgomery_point(E,x,z=1):
	if z==0:
		return E(0)
	elif z!=1:
		x=x/z

	a_inv=E.a_invariants()
		
	A =a_inv[1]
	if a_inv != (0,A,0,1,0):
		raise ValueError("The elliptic curve E is not in the Montgomery model.")

	y2=x*(x**2+A*x+1)
	if not is_square(y2):
		raise ValueError("The Montgomery point is not on the base field.")
	else:
		y=y2.sqrt()
		return E([x,y,1])

def null_point_to_montgomery_coeff(a,b):
	return -2*(a**4+b**4)/(a**4-b**4)






