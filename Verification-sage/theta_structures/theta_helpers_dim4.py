from sage.all import *
import itertools

## Index management
@cached_function
def index_to_multindex(k):
	r"""
	Input: 
	- k: integer between 0 and 15.

	Output: binary decomposition of k.
	"""
	L_ind=[]
	l=k
	for i in range(4):
		L_ind.append(l%2)
		l=l//2
	return tuple(L_ind)

@cached_function
def multindex_to_index(*args):
	r"""
	Input: 4 elements i0,i1,i2,i3 in {0,1}.

	Output: k=i0+2*i1+4*i2+8*i3.
	"""
	if len(args)==4:
		i0,i1,i2,i3=args
	else:
		i0,i1,i2,i3=args[0]
	return i0+2*i1+4*i2+8*i3

@cached_function
def scal_prod(i,j):
	r"""
	Input: Two integers i and j in {0,...,15}.

	Output: Scalar product of the bits of i and j mod 2.
	"""
	return (int(i)&int(j)).bit_count()%2

def act_point(P,I,J):
	r"""
	Input:
	- P: a point with 16 coordinates.
	- I, J: two 4-tuples of indices in {0,1}.

	Output: the action of (I,\chi_J) on P given by:
	(I,\chi_J)*P=(\chi_J(I+K)^{-1}P[I+K])_K
	"""
	Q=[]
	i=multindex_to_index(I)
	j=multindex_to_index(J)
	for k in range(16):
		ipk=i^k
		Q.append((-1)**scal_prod(ipk,j)*P[ipk])
	return Q

## Product of theta points
def product_theta_point_dim4(P0,P1,P2,P3):
	# Computes the product theta coordinates of a product of 4 elliptic curves.
	P=[0 for k in range(16)]
	for i0, i1, i2, i3 in itertools.product([0,1],repeat=4):
		P[multindex_to_index(i0,i1,i2,i3)]=P0[i0]*P1[i1]*P2[i2]*P3[i3]
	return P

def product_theta_point_dim2_dim4(P0,P1):
	# Computes the product theta coordinates of a product of 2 abelian surfaces.
	P=[0 for k in range(16)]
	for i0, i1, i2, i3 in itertools.product([0,1],repeat=4):
		P[multindex_to_index(i0,i1,i2,i3)]=P0[i0+2*i1]*P1[i2+2*i3]
	return P

## 4-dimensional product Theta point to 1-dimensional Theta points
def product_to_theta_points_dim4(P):
	Fp2=P[0].parent()
	d_Pi={0:None,1:None,2:None,3:None}
	d_index_ratios={0:None,1:None,2:None,3:None}# Index of numertors/denominators to compute the theta points.
	for k in range(4):
		is_zero=True# theta_1(Pk)=0
		for I in itertools.product([0,1],repeat=3):
			J=list(I)
			J.insert(k,1)
			j=multindex_to_index(*J)
			if P[j]!=0:
				is_zero=False
				d_index_ratios[k]=[j^(2**k),j]
				break
		if is_zero:
			d_Pi[k]=(Fp2(1),Fp2(0))
	L_num=[]
	L_denom=[]
	d_index_num_denom={}
	for k in range(4):
		if d_Pi[k]==None:# Point has non-zero coordinate theta_1
			d_index_num_denom[k]=len(L_num)
			L_num.append(P[d_index_ratios[k][0]])
			L_denom.append(P[d_index_ratios[k][1]])
	L_denom=batch_inversion(L_denom)
	for k in range(4):
		if d_Pi[k]==None:
			d_Pi[k]=(L_num[d_index_num_denom[k]]*L_denom[d_index_num_denom[k]],Fp2(1))
	return d_Pi[0],d_Pi[1],d_Pi[2],d_Pi[3]

## 4-dimensional product Theta point to 2-dimensional Theta points
def product_to_theta_points_dim4_dim2(P):
	Fp2=P[0].parent()
	k0=0# Index of the denominator k0=multindex_to_index(I0,J0) and
	# we divide by theta_{I0,J0}=theta_{I0}*theta_{J0}!=0
	while k0<=15 and P[k0]==0:
		k0+=1
	i0, j0 = k0%4, k0//4
	inv_theta_k0=1/P[k0]
	theta_P1=[]
	theta_P2=[]
	for i in range(4):
		if i==i0:
			theta_P1.append(1)
		else:
			theta_P1.append(P[i+4*j0]*inv_theta_k0)
	for j in range(4):
		if j==j0:
			theta_P2.append(1)
		else:
			theta_P2.append(P[i0+4*j]*inv_theta_k0)
	return theta_P1, theta_P2





## Usual theta transformations
@cached_function
def inv_16(F):
	return 1/F(16)

def hadamard2(x,y):
    return (x+y, x-y)

def hadamard4(x,y,z,t):
    x,y=hadamard2(x,y)
    z,t=hadamard2(z,t)
    return (x+z, y+t, x-z, y-t)

def hadamard8(a,b,c,d,e,f,g,h):
    a,b,c,d=hadamard4(a,b,c,d)
    e,f,g,h=hadamard4(e,f,g,h)
    return (a+e, b+f, c+g, d+h, a-e, b-f, c-g, d-h)

def hadamard16(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p):
    a,b,c,d,e,f,g,h=hadamard8(a,b,c,d,e,f,g,h)
    i,j,k,l,m,n,o,p=hadamard8(i,j,k,l,m,n,o,p)
    return (a+i, b+j, c+k, d+l, e+m, f+n, g+o, h+p, a-i, b-j, c-k, d-l, e-m, f-n, g-o, h-p)

def hadamard_inline(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p):
    a,b=a+b,a-b
    c,d=c+d,c-d
    e,f=e+f,e-f
    g,h=g+h,g-h
    i,j=i+j,i-j
    k,l=k+l,k-l
    m,n=m+n,m-n
    o,p=o+p,o-p
    a,b,c,d=a+c,b+d,a-c,b-d
    e,f,g,h=e+g,f+h,e-g,f-h
    i,j,k,l=i+k,j+l,i-k,j-l
    m,n,o,p=m+o,n+p,m-o,n-p
    a,b,c,d,e,f,g,h=a+e, b+f, c+g, d+h, a-e, b-f, c-g, d-h
    i,j,k,l,m,n,o,p=i+m,j+n,k+o,l+p,i-m,j-n,k-o,l-p
    return (a+i, b+j, c+k, d+l, e+m, f+n, g+o, h+p, a-i, b-j, c-k, d-l, e-m, f-n, g-o, h-p)

def hadamard(P):
    return hadamard16(*P)
    #return hadamard_inline(*P)

def hadamard_inverse(P):
	H_inv_P=[]
	C=inv_16(P[0].parent())
	for j in range(16):
		HP.append(0)
		for k in range(16):
			if scal_prod(k,j)==0:
				H_inv_P[j]+=P[k]
			else:
				H_inv_P[j]-=P[k]
		H_inv_P[j]=H_inv_P[j]*C
	return H_inv_P

def squared(P):
	return [x**2 for x in P]

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


## Functions to handle zero theta coordinates
def find_zeros(P):
	L_ind_zeros=[]
	for i in range(16):
		if P[i]==0:
			L_ind_zeros.append(i)
	return L_ind_zeros

def find_translates(L_ind_zeros):
	L_ind_non_zero=[]
	L_ind_origin=L_ind_zeros.copy()
	
	for i in range(16):
		if i not in L_ind_zeros:
			L_ind_non_zero.append(i)

	L_trans=[]
	while L_ind_origin!=[]:
		n_target_max=0
		ind_trans_max=0
		for k in range(16):
			trans=[x^k for x in L_ind_origin]
			n_target=0
			for y in trans:
				if y in L_ind_non_zero:
					n_target+=1
			if n_target>n_target_max:
				n_target_max=n_target
				ind_trans_max=k
		L_trans.append(ind_trans_max)
		L_ind_remove=[]
		for x in L_ind_origin:
			if x^ind_trans_max in L_ind_non_zero:
				L_ind_remove.append(x)
		for x in L_ind_remove:
			L_ind_origin.remove(x)
	return L_trans




