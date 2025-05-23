from sage.all import *

def solve_norm_equation(u,p,margin=100):
	r""" 
	INPUT:
	- u: odd integer close to sqrt(p).
	- p: prime numner.
	- magin: integer to fix the bounds (optional, default = 100).

	OUTPUT:
	- a, b, c, d, f such that a**2+c**2+p*(b**2+d**2)=u*(2**f-u) with a, b and d odd.
	"""
	bound=ceil(margin*log(p)/log(2)*p)
	f=0
	two_f=1
	N=u*(two_f-u)
	while N<bound:
		f+=1
		two_f*=2
		N=u*(two_f-u)

	Q0=BinaryQF([1,0,1])
	bound_bd=floor(sqrt(margin*log(p)/log(2))/2)
	for b in range(bound_bd):
		if b%2==1:
			for d in range(bound_bd):
				if d%2==1:
					M=N-p*(b**2+d**2)
					if M%4==1 and is_pseudoprime(M):
						a,c=Q0.solve_integer(M,algorithm='cornacchia')
						if a%2==0:
							a,c=c,a
						return a,b,c,d,f




	

