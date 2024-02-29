from sage.all import *

def batch_inversion(L):
	r"""Does n inversions in 3(n-1)M+1I.

	Input:
	- L: list of elements to invert.

	Output:
	- [1/x for x in L]
	"""
	# Given L=[a0,...,an]
	# Computes multiples=[a0, a0.a1, ..., a0...an]
	multiples=[L[0]]
	for ai in L[1:]:
		multiples.append(multiples[-1]*ai)

	# Computes inverses=[1/(a0...an),...,1/a0]
	inverses=[1/multiples[-1]]
	for i in range(1,len(L)):
		inverses.append(inverses[-1]*L[-i])

	# Finally computes [1/a0,...,1/an]
	result=[inverses[-1]]
	for i in range(2,len(L)+1):
		result.append(inverses[-i]*multiples[i-2])
	return result

def product_theta_point(theta1,theta2):
    a,b=theta1
    c,d=theta2
    return [a*c,b*c,a*d,b*d]