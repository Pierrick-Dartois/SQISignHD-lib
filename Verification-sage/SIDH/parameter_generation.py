from sage.all import *

d_params={'p434':{'e2':216,'e3':137},
'p503':{'e2':250,'e3':159},
'p610':{'e2':305,'e3':192},
'p751':{'e2':372,'e3':239}}

L_primes=prime_range(1000)

def lazy_factor(N):
	r"""
	We return the factorization of N if it is easy 
	to factor and can be decomposed as a sum of two squares 
	and False otherwise.
	"""

	M=N
	L_fact=[]
	for p in L_primes:
		e_p=0
		while M%p==0:
			M=M//p
			e_p+=1
		if e_p>0:
			if p==2 or p%4==1 or (p%4==3 and e_p%2==0):
				L_fact.append((p,e_p))
			else:
				return False
	if is_prime(M) and M%4==1:
		L_fact.append((M,1))
		return L_fact
	else:
		return False

def lazy_cornacchia(N):
	r"""
	Returns integers a1, a2 such that N=a1**2+a2**2 if they exist and are easy to find.
	Returns False otherwise.
	"""
	Q0=BinaryQF([1,0,1])

	L_fact=lazy_factor(N)
	if L_fact:
		L_a1a2=[]
		sq_factor=1
		for fact in L_fact:
			p,e_p=fact
			if e_p%2==0:
				sq_factor *= p**(e_p//2)
			else:
				a1,a2=Q0.solve_integer(p, algorithm='cornacchia')
				L_a1a2.append((a1,a2,e_p))

		gauss_integer=1
		for elt in L_a1a2:
			a1,a2,e_p=elt
			gauss_integer*=(a1+a2*I)**e_p
		return sq_factor*gauss_integer[0],sq_factor*gauss_integer[1]
	else:
		return False

def find_embedding_params(e2,e3):
	r"""
	Given exponents e2 and e3 of the accessible 2**e2 and 3**e3-torsion
	respectively, this function finds integers e, a1, a2 such that 
	2**e=3**(e3+eps)+a1**2+a2**2 with eps=1-e3%2.
	Adding eps increases greatly the probability to find a1 and a2 as 
	it quarantees that 2**e-3**(e3+eps)=1 mod 4.
	"""

	e=e2
	if e3%2==1:
		q=3**e3
	else:
		q=3**(e3+1)
	a=2**e-q
	a1a2=lazy_cornacchia(a)
	while not a1a2:
		e+=1
		a=2**e-q
		a1a2=lazy_cornacchia(a)
	a1,a2=a1a2

	assert q+a1**2+a2**2==2**e
	return e,a1,a2

def read_params_SIKE(file):
	d_params_ret={}
	with open(file,'r',encoding='utf-8') as f:
		char=f.read()
		L=char.split('}')
		for i in range(len(L)):
			if i==0:
				ID=L[0][2:6]
				Li=L[0][10:].split(',')
			elif len(L[i])>0:
				ID=L[i][3:7]
				Li=L[i][11:].split(',')
			else:
				continue
			d_params_ret[ID]={}
			for x in Li:
				Lx=x.split(':')
				d_params_ret[ID][Lx[0].split(' ')[-1][1:-1]]=int(Lx[1])
		return d_params_ret


if __name__=="__main__":
	d_params_save={}
	for param in d_params:
		e,a1,a2=find_embedding_params(d_params[param]['e2'],d_params[param]['e3'])
		d_params_save[param]={'e2':d_params[param]['e2'],'e3':d_params[param]['e3'],'e':e,'a1':a1,'a2':a2}
	with open('SIDH/parameters_SIKE_NIST.txt','w',encoding='utf-8') as f:
		f.write(str(d_params_save))
	d_params_ret=read_params_SIKE('SIDH/parameters_SIKE_NIST.txt')

