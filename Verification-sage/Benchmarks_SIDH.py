from sage.all import *
from SIDH.parameter_generation import read_params_SIKE
from utilities.supersingular import torsion_basis, weil_pairing_pari
from utilities.discrete_log import ell_discrete_log_pari
from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only, random_isogeny_x_only
from basis_change.canonical_basis_dim1 import make_canonical
from Tests import random_walk
from utilities.strategy import precompute_strategy_with_first_eval
from isogenies.Kani_endomorphism import KaniEndoHalf
from theta_structures.Tuple_point import TuplePoint
from time import time

def SIDH_key_exchange_bench(pub_params,NA,NB,fast=True):
	E1,PA,QA,PB,QB=pub_params

	## Alice isogeny
	# Alice's secret key
	sa=randint(0,NA-1)
	# ker(phiA)=PA+sa*QA
	phiA, EA=isogeny_from_scalar_x_only(E1, NA, sa, basis=(PA,QA))
	phiA_PB,phiA_QB=evaluate_isogeny_x_only(phiA, PB, QB, NB, NA)

	# Alice's public key
	pubA=(EA,phiA_PB,phiA_QB)

	## Bob isogeny
	# Bob's secret key
	sb=randint(0,NB-1)
	# ker(phiB)=PB+sb*QB
	phiB, EB=isogeny_from_scalar_x_only(E1, NB, sb, basis=(PB,QB))
	phiB_PA,phiB_QA=evaluate_isogeny_x_only(phiB, PA, QA, NA, NB)

	# Bob's public key
	pubB=(EB,phiB_PA,phiB_QA)

	## Action of Alice on EB
	psiA,EBA=isogeny_from_scalar_x_only(EB, NA, sa, basis=(phiB_PA,phiB_QA))

	## Action of Bob on EA
	if fast:
		# Skipped with fast option
		EAB=EBA
	else:
		psiB,EAB=isogeny_from_scalar_x_only(EA, NB, sb, basis=(phiA_PB,phiA_QB))

		assert EAB.j_invariant()==EBA.j_invariant()
	
	return pubA,pubB,EAB

def SIDH_key_recovery_attack_bench(params,pub_params,pubA,pubB,EAB=None):
	E1,PA,QA,PB,QB=pub_params
	EA,phiA_PB,phiA_QB=pubA
	EB,phiB_PA,phiB_QA=pubB

	e2=params['e2']
	e3=params['e3']
	e=params['e']
	a1=params['a1']
	a2=params['a2']
	m=params['m']
	strategy1=params['strategy1']
	strategy2=params['strategy2']

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

	# Dimension 4 embedding of phipB=(psi*)phiB
	F=KaniEndoHalf(PA,QA,phipB_PA,phipB_QA,q,a1,a2,e,e2,strategy1,strategy2)

	# Evaluation of phipB on E1[3**e3] to recover ker(phiB)
	T=TuplePoint(PB,E1(0),EBp(0),EBp(0))
	FT=F(T)
	phipB_PB=FT[2]

	U=TuplePoint(QB,E1(0),EBp(0),EBp(0))
	FU=F(U)
	phipB_QB=FU[2]

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
	
	# Find the shared secret by acting on Alice's public key.
	psiB,EABp=isogeny_from_scalar_x_only(EA, NB, sb, basis=(phiA_PB,phiA_QB))

	if EAB!=None:
		assert EABp.j_invariant()==EAB.j_invariant()

	return EABp


def benchmark_attacks(params,N_iter):
	## Public parameters for all attack instances
	e2=params['e2']
	e3=params['e3']
	e=params['e']
	a1=params['a1']
	a2=params['a2']

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

	## Additional params for the attack
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

	params['m']=m
	
	f1=ceil(e/2)
	f2=e-f1

	# Strategy precomputation
	strategy1=precompute_strategy_with_first_eval(f1,m,M=1,S=0.8,I=100)
	if f2==f1:
		strategy2=strategy1
	else:
		strategy2=precompute_strategy_with_first_eval(f2,m,M=1,S=0.8,I=100)

	params['strategy1']=strategy1
	params['strategy2']=strategy2

	L_time_protocol=[]
	L_time_attack=[]

	for i in range(N_iter):
		print(i)
		t0=time()
		pubA,pubB,EAB=SIDH_key_exchange_bench(pub_params,NA,NB,fast=True)
		t1=time()
		EABp=SIDH_key_recovery_attack_bench(params,pub_params,pubA,pubB,EAB=None)
		t2=time()
		assert EABp.j_invariant()==EAB.j_invariant()
		L_time_protocol.append(t1-t0)
		L_time_attack.append(t2-t1)

	return L_time_protocol, L_time_attack


if __name__=="__main__":
	d_params=read_params_SIKE('SIDH/parameters_SIKE_NIST.txt')

	d_results={}
	N_iter=100
	for x in d_params:
		print("Benchmarking attack against SIDH {}.".format(x))
		t0=time()
		L_time_protocol, L_time_attack=benchmark_attacks(d_params[x],N_iter)
		d_results[(x,'protocol')]=L_time_protocol
		d_results[(x,'attack')]=L_time_attack
		t1=time()
		print(t1-t0)

	with open("Benchmarking_results_SIDH.csv",'w',encoding='utf-8') as f:
		for i in range(2):
			line=""
			for r in d_results:
				line+=str(r[i])+","
			f.write(line+"\n")
		for i in range(N_iter):
			line=""
			for r in d_results:
				line+=str(d_results[r][i])+","
			f.write(line+"\n")






