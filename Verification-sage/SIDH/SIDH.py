from sage.all import *
from SIDH.parameter_generation import read_params_SIKE
from utilities.supersingular import torsion_basis
from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only
from basis_change.canonical_basis_dim1 import make_canonical
from Tests import random_walk

def SIDH_key_exchange(params):
	e2=params['e2']
	e3=params['e3']
	e=params['e']
	a1=params['a1']
	a2=params['a2']

	p=2**e2*3**e3-1

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	i=Fp2.gen()
	
	E0=EllipticCurve(Fp2,[1,0])

	# Starting curve
	E1=random_walk(E0,ZZ(2**(e2-1)*3**(e3-1)))

	NA=ZZ(2**e2)
	NB=ZZ(3**e3)

	# Basis of Alice and Bob
	PA,QA=torsion_basis(E1,NA)
	PA,QA,_,_,_=make_canonical(PA,QA,NA)# QB is above (0,0) which should not be in the kernel
	PB,QB=torsion_basis(E1,NB)

	## Alice isogeny
	# Alice's secret key
	sa=randint(0,NA-1)

	phiA, EA=isogeny_from_scalar_x_only(E1, NA, sa, basis=(PA,QA))
	phiA_PB,phiA_QB=evaluate_isogeny_x_only(phiA, PB, QB, NB, NA)

	# Alice's public key
	pubA=(EA,phiA_PB,phiA_QB)

	## Bob isogeny
	# Bob's secret key
	sb=randint(0,NB-1)

	phiB, EB=isogeny_from_scalar_x_only(E1, NB, sb, basis=(PB,QB))
	phiB_PA,phiB_QA=evaluate_isogeny_x_only(phiB, PA, QA, NA, NB)

	# Bob's public key
	pubB=(EB,phiB_PA,phiB_QA)

	## Action of Alice on EB
	psiA,EBA=isogeny_from_scalar_x_only(EB, NA, sa, basis=(phiB_PA,phiB_QA))

	## Action of Bob on EA
	psiB,EAB=isogeny_from_scalar_x_only(EA, NB, sb, basis=(phiA_PB,phiA_QB))

	## Shared secret key
	assert EAB.j_invariant()==EBA.j_invariant()

	return E1,pubA,pubB,EAB

if __name__=="__main__":
	d_params=read_params_SIKE('SIDH/parameters_SIKE_NIST.txt')
	E1,pubA,pubB,EAB=SIDH_key_exchange(d_params['p434'])



