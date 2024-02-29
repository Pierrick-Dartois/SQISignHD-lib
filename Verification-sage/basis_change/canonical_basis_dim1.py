from sage.all import *
from utilities.discrete_log import weil_pairing_pari, discrete_log_pari

def last_four_torsion(E):
	a_inv=E.a_invariants()
	A =a_inv[1]
	if a_inv != (0,A,0,1,0):
		raise ValueError("The elliptic curve E is not in the Montgomery model.")
	y2=A-2
	y=y2.sqrt()
	return E([-1,y,1])


def make_canonical(P,Q,A,preserve_pairing=False):
	r"""
	Input:
	- P,Q: a basis of E[A].
	- A: an integer divisible by 4.
	- preserve_pairing: boolean indicating if we want to preserve pairing at level 4.
	
	Output:
	- P1,Q1: basis of E[A].
	- U1,U2: basis of E[4] induced by (P1,Q1) ((A//4)*P1=U1, (A//4)*Q1=U2) such that U2[0]=-1 
	and e_4(U1,U2)=i if not preserve_pairing and e_4(U1,U2)=e_4((A//4)*P,(A//4)*Q) if preserve_pairing.
	- M: base change matrix (in row convention) from (P1,Q1) to (P,Q).
	
	We say that (U1,U2) is canonical and that (P1,Q1) induces or lies above a canonical basis. 
	"""
	E=P.curve()
	Fp2=E.base_ring()
	i=Fp2.gen()

	assert i**2==-1

	T2=last_four_torsion(E)
	V1=(A//4)*P
	V2=(A//4)*Q
	U1=V1
	U2=V2
	
	a1=discrete_log_pari(weil_pairing_pari(U1,T2,4),i,4)
	b1=discrete_log_pari(weil_pairing_pari(U2,T2,4),i,4)

	if a1%2!=0:
		c1=inverse_mod(a1,4)
		d1=c1*b1
		P1=P
		Q1=Q-d1*P
		U1,U2=U1,U2-d1*U1
		M=matrix(ZZ,[[1,0],[d1,1]])
	else:
		c1=inverse_mod(b1,4)
		d1=c1*a1
		P1=Q
		Q1=P-d1*Q
		U1,U2=U2,U1-d1*U2
		M=matrix(ZZ,[[d1,1],[1,0]])

	if preserve_pairing:
		e4=weil_pairing_pari(V1,V2,4)
	else:
		e4=i
	
	if weil_pairing_pari(U1,U2,4)!=e4:
		U2=-U2
		Q1=-Q1
		M[0,1]=-M[0,1]
		M[1,1]=-M[1,1]
	
	assert (A//4)*P1==U1
	assert (A//4)*Q1==U2
	assert weil_pairing_pari(U1,U2,4)==e4
	assert M[0,0]*P1+M[0,1]*Q1==P
	assert M[1,0]*P1+M[1,1]*Q1==Q
	
	return P1,Q1,U1,U2,M
