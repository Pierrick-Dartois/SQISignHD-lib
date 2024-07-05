from sage.all import *
from SIDH.parameter_generation import read_params_SIKE
from utilities.supersingular import torsion_basis, weil_pairing_pari
from utilities.discrete_log import ell_discrete_log_pari #discrete_log_pari
from utilities.order import has_order_D
from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only, random_isogeny_x_only
from basis_change.canonical_basis_dim1 import make_canonical
from Tests import random_walk
from utilities.strategy import precompute_strategy_with_first_eval
from isogenies.Kani_endomorphism import KaniEndoHalf
from theta_structures.Tuple_point import TuplePoint
from time import time

# CLI imports
import argparse
from pathlib import Path

def SIDH_key_exchange(params):
	t1=time()
	print("=================")
	print("# SIDH protocol #")
	print("=================\n")

	print("Public parameter generation:")
	print("- Starting curve E1.")
	print("- <PA,QA>=E1[2**e2].")
	print("- <PB,QB>=E1[3**e3].")

	e2=params['e2']
	e3=params['e3']

	p=2**e2*3**e3-1

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	
	E0=EllipticCurve(Fp2,[1,0])

	# Starting curve
	E1=random_walk(E0,ZZ(2**(e2-1)*3**(e3-1)))

	NA=ZZ(2**e2)
	NB=ZZ(3**e3)

	# Basis of Alice and Bob
	PA,QA=torsion_basis(E1,NA)
	PA,QA,_,_,_=make_canonical(PA,QA,NA)# QB is above (0,0) which should not be in the kernel
	PB,QB=torsion_basis(E1,NB)

	pub_params=(E1,PA,QA,PB,QB)
	t2=time()
	print("Parameter generation time: {} s\n".format(t2-t1))

	## Alice isogeny
	print("Alice: 2**e2-isogeny phiA: E1-->EA")
	# Alice's secret key
	sa=randint(0,NA-1)
	# ker(phiA)=PA+sa*QA
	phiA, EA=isogeny_from_scalar_x_only(E1, NA, sa, basis=(PA,QA))
	phiA_PB,phiA_QB=evaluate_isogeny_x_only(phiA, PB, QB, NB, NA)

	# Alice's public key
	pubA=(EA,phiA_PB,phiA_QB)

	## Bob isogeny
	print("Bob: 3**e3-isogeny phiB: E1-->EB")
	# Bob's secret key
	sb=randint(0,NB-1)
	# ker(phiB)=PB+sb*QB
	phiB, EB=isogeny_from_scalar_x_only(E1, NB, sb, basis=(PB,QB))
	phiB_PA,phiB_QA=evaluate_isogeny_x_only(phiB, PA, QA, NA, NB)

	# Bob's public key
	pubB=(EB,phiB_PA,phiB_QA)

	## Action of Alice on EB
	print("Alice: 2**e2-isogeny psiA: EB-->EBA")
	psiA,EBA=isogeny_from_scalar_x_only(EB, NA, sa, basis=(phiB_PA,phiB_QA))

	## Action of Bob on EA
	print("Bob: 3**e3-isogeny psiB: EA-->EAB")
	psiB,EAB=isogeny_from_scalar_x_only(EA, NB, sb, basis=(phiA_PB,phiA_QB))

	## Shared secret key
	print("Shared secret key: j(EAB)==j(EBA)?")
	print(EAB.j_invariant()==EBA.j_invariant())

	t3=time()
	print("Protocol time: {} s".format(t3-t2))

	return pub_params,pubA,pubB,EAB

def SIDH_key_recovery_attack(params,pub_params,pubA,pubB,EAB=None):
	t1=time()
	print("============================")
	print("# SIDH key recovery attack #")
	print("============================\n")
	
	E1,PA,QA,PB,QB=pub_params
	EA,phiA_PB,phiA_QB=pubA
	EB,phiB_PA,phiB_QA=pubB

	e2=params['e2']
	e3=params['e3']
	e=params['e']
	a1=params['a1']
	a2=params['a2']

	# Degree of Bob's isogeny phiB to recover (up to a tweak) and torsion image
	if e3%2==1:
		q=3**e3
		phipB_PA=phiB_PA
		phipB_QA=phiB_QA
		EBp=EB
	else:
		# Tweaking phiB
		# We add a 3-isogeny to phiB in accordance with the tweak in parameter generation
		q=3**(e3+1)
		psi,EBp=random_isogeny_x_only(EB,3)
		phipB_PA,phipB_QA=evaluate_isogeny_x_only(psi, phiB_PA, phiB_QA, ZZ(2**e2), 3)

	assert q+a1**2+a2**2==2**e

	# m=max(v_2(a1),v_2(a2))
	if a1%2==0:
		ai_div=a1
	else:
		ai_div=a2

	ai_div=ai_div//2
	m=1
	while ai_div%2==0:
		ai_div=ai_div//2
		m+=1
	
	f1=ceil(e/2)
	f2=e-f1

	# Strategy precomputation
	strategy1=precompute_strategy_with_first_eval(f1,m,M=1,S=0.8,I=100)
	if f2==f1:
		strategy2=strategy1
	else:
		strategy2=precompute_strategy_with_first_eval(f2,m,M=1,S=0.8,I=100)

	t2=time()
	print("Precomputations: {} s".format(t2-t1))

	# Dimension 4 embedding of phipB=(psi*)phiB
	F=KaniEndoHalf(PA,QA,phipB_PA,phipB_QA,q,a1,a2,e,e2,strategy1,strategy2)
	t3=time()
	print("Dimension 4 embedding: {} s".format(t3-t2))


	# Evaluation of phipB on E1[3**e3] to recover ker(phiB)
	T=TuplePoint(PB,E1(0),EBp(0),EBp(0))
	FT=F(T)
	phipB_PB=FT[2]

	U=TuplePoint(QB,E1(0),EBp(0),EBp(0))
	FU=F(U)
	phipB_QB=FU[2]

	t4=time()
	print("Evaluation of Bob's isogeny phiB: {} s".format(t4-t3))

	# DL computation to find sb
	# ker(phiB)=PB+sb*QB so sb=-DL(phiB(PB),phiB(QB))
	NB=ZZ(3**e3)
	backtrack=False
	if (NB//3)*phipB_QB!=EBp(0):
		# When ker(psi)\cap<phiB(QB)>={0}
		sb=ell_discrete_log_pari(EBp,phipB_PB,phipB_QB,NB)
	else:
		# When ker(psi)\cap<phiB(QB)>!={0} (psi backtracks)
		backtrack=True
		sb=ell_discrete_log_pari(EBp,phipB_PB,phipB_QB,ZZ(NB//3))
	t5=time()
	print("Discrete log to find ker(phiB): {} s".format(t5-t4))
	

	# Find the correct sb
	if backtrack:
		# sb determined modulo NB//3=3**(e3-1) and up to a sign
		sb_candidates=[sb,-sb,sb+NB//3,-sb+NB//3,sb+2*NB//3,-sb+2*NB//3]
	else:
		# sb determined modulo NB=3**e3 up to a sign
		sb_candidates=[sb,-sb]

	found_sb=False
	for s in sb_candidates:
		phiB, EC=isogeny_from_scalar_x_only(E1, NB, s, basis=(PB,QB))
		if EC.j_invariant()==EB.j_invariant():
			found_sb=True
			sb=s
			break
	t6=time()
	print("Time to recompute phiB: {} s".format(t6-t5))
	print("Curve EB has been recovered?")
	print(found_sb)
	
	# Find the shared secret by acting on Alice's public key.
	psiB,EABp=isogeny_from_scalar_x_only(EA, NB, sb, basis=(phiA_PB,phiA_QB))

	t7=time()
	print("Computation of shared secret: {} s".format(t7-t6))
	if EAB!=None:
		print("Shared secret has been recovered?")
		print(EABp.j_invariant()==EAB.j_invariant())

	print("Total time: {} s".format(t7-t1))

	return EABp

if __name__=="__main__":
	d_params=read_params_SIKE('SIDH/parameters_SIKE_NIST.txt')

	parser = argparse.ArgumentParser()

	# To display stored parameters
	parser.add_argument("-d","--display",action="store_true")
	# To run the protocol or the protcol and the attack
	parser.add_argument("--Protocol",action="store_true")
	parser.add_argument("--Protocol_and_Attack",action="store_true")
	# Parameter
	parser.add_argument("-p")

	args = parser.parse_args()

	if args.display:
		print("============================")
		print("# SIKE available parameters #")
		print("============================\n")
		for x in d_params:
			print("Prime SIKE {}=2**{}*3**{}-1".format(x,d_params[x]['e2'],d_params[x]['e3']))
	elif args.Protocol:
		x=args.p
		print("Attack on SIDH with {}=2**{}*3**{}-1\n".format(x,d_params[x]['e2'],d_params[x]['e3']))

		pub_params,pubA,pubB,EAB=SIDH_key_exchange(d_params[x])
	elif args.Protocol_and_Attack:
		x=args.p
		print("Attack on SIDH with {}=2**{}*3**{}-1\n".format(x,d_params[x]['e2'],d_params[x]['e3']))

		pub_params,pubA,pubB,EAB=SIDH_key_exchange(d_params[x])
		print("\n")
		EABp=SIDH_key_recovery_attack(d_params[x],pub_params,pubA,pubB,EAB)



