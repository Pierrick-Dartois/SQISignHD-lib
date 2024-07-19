from sage.all import *
#from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
from time import time

from parameters.parameter_generation import read_params
from utilities.supersingular import random_point, compute_point_order_D, torsion_basis
from isogenies.Kani_endomorphism import KaniEndo, KaniEndoHalf
from theta_structures.Tuple_point import TuplePoint
from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only
from utilities.strategy import precompute_strategy_with_first_eval, precompute_strategy_with_first_eval_and_splitting
from basis_change.canonical_basis_dim1 import make_canonical
from Tests import random_walk

from pathlib import Path

def benchmark_kani_endo(N_iter,l_B,e_A,e_B,a1,a2,f,f_A,f_B,p,m=None):

	q=l_B**e_B

	# Recovering m from a2
	if a1%2==0:
		ai_div=a1
	else:
		ai_div=a2

	ai_div=ai_div//2
	m_comp=1
	while ai_div%2==0:
		ai_div=ai_div//2
		m_comp+=1

	if m==None:
		m=m_comp

	## Setting up the base field, starting elliptic curve, torsion basis

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	i=Fp2.gen()
	
	E0=EllipticCurve(Fp2,[1,0])

	# Generating a random elliptic curve in the supersingular isogeny graph.
	N=ZZ(2**(f_A-1)*l_B**(f_B-1))
	E1=random_walk(E0,N)

	Basis_sigma=torsion_basis(E1,q)

	## Setting up other parameters
	strategy=precompute_strategy_with_first_eval_and_splitting(e_A,m,M=1,S=0.8,I=100)

	L_time_compute=[]
	L_time_eval=[]

	for k in range(N_iter):
		scalar=randint(0,q-1)

		# q-isogeny of kernel <Basis_sigma[0]+scalar*Basis_sigma[1]>
		sigma, E2 = isogeny_from_scalar_x_only(E1, ZZ(q), scalar, basis=Basis_sigma)

		P1,Q1=torsion_basis(E1,2**f_A)
		R2,S2=evaluate_isogeny_x_only(sigma, P1, Q1, ZZ(2**f_A), ZZ(q))

		t0=time()

		F=KaniEndo(P1,Q1,R2,S2,q,a1,a2,e_A,f_A,strategy)

		t1=time()

		L_time_compute.append(t1-t0)

		T=TuplePoint(P1,E1(0),E2(0),E2(0))

		t2=time()

		FT=F(T)

		t3=time()

		L_time_eval.append(t3-t2)

	return L_time_compute, L_time_eval


def benchmark_kani_endo_half(N_iter,l_B,e_A,e_B,a1,a2,f,f_A,f_B,p,m=None):

	q=l_B**e_B

	# Recovering m from a2
	if a1%2==0:
		ai_div=a1
	else:
		ai_div=a2

	ai_div=ai_div//2
	m_comp=1
	while ai_div%2==0:
		ai_div=ai_div//2
		m_comp+=1

	if m==None:
		m=m_comp

	## Setting up the base field, starting elliptic curve, torsion basis

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	i=Fp2.gen()
	
	E0=EllipticCurve(Fp2,[1,0])

	# Generating a random elliptic curve in the supersingular isogeny graph.
	N=ZZ(2**(f_A-1)*l_B**(f_B-1))
	E1=random_walk(E0,N)

	Basis_sigma=torsion_basis(E1,q)

	## Setting up other parameters

	f=ceil(e_A/2)+2

	e1=ceil(e_A/2)
	e2=e_A-e1

	strategy1=precompute_strategy_with_first_eval(e1,m,M=1,S=0.8,I=100)
	if e2==e1:
		strategy2=strategy1
	else:
		strategy2=precompute_strategy_with_first_eval(e2,m,M=1,S=0.8,I=100)


	L_time_compute=[]
	L_time_eval=[]

	for k in range(N_iter):
		scalar=randint(0,q-1)

		sigma, E2 = isogeny_from_scalar_x_only(E1, ZZ(q), scalar, basis=Basis_sigma)

		P1,Q1=torsion_basis(E1,2**f)
		R2,S2=evaluate_isogeny_x_only(sigma, P1, Q1, ZZ(2**f), ZZ(q))

		t0=time()

		F=KaniEndoHalf(P1,Q1,R2,S2,q,a1,a2,e_A,f,strategy1,strategy2)

		t1=time()

		L_time_compute.append(t1-t0)

		T=TuplePoint(P1,E1(0),E2(0),E2(0))

		t2=time()

		FT=F(T)

		t3=time()

		L_time_eval.append(t3-t2)

	return L_time_compute, L_time_eval

def benchmark_dim_one(N_iter,l_B,e_A,e_B,a1,a2,f,f_A,f_B,p,m=None):
	## Setting up the base field, starting elliptic curve, torsion basis

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	i=Fp2.gen()
	
	E0=EllipticCurve(Fp2,[1,0])

	# Generating a random elliptic curve in the supersingular isogeny graph.
	N=ZZ(2**(f_A-1)*l_B**(f_B-1))
	E1=random_walk(E0,N)

	D=ZZ(2**e_A)
	P,Q=torsion_basis(E1,D)
	P,Q,_,_,_=make_canonical(P,Q,D)# Q is above (0,0) which should not be in the kernel

	L_time_compute=[]
	L_time_eval=[]
	for k in range(N_iter):
		scalar=randint(0,2**e_A-1)

		t0=time()
		# Isogeny computation
		phi, E2 = isogeny_from_scalar_x_only(E1, D, scalar, basis=(P,Q))
		t1=time()
		# Isogeny evaluation
		L1=phi.domain()
		xP=L1(P[0])
		xphiP=phi(xP)
		phiP=xphiP.curve_point()

		t2=time()

		L_time_compute.append(t1-t0)
		L_time_eval.append(t2-t1)

	return L_time_compute, L_time_eval


if __name__=="__main__":
	## Retrieving parameter files
	target_dir=Path("parameters")

	# Dictionnary of list of params indexed by l_B (3 or 7)
	d_L_params={}

	for file in target_dir.iterdir():
		# Files containing parameters are .txt files of the form "parameters_{l_B}.txt"
		if str(file)[-3::]=="txt":
			try:
				L=str(file)[:-4].split('_')
				l_B=int(L[-1])
				d_L_params[l_B]=read_params(file)
			except:
				raise ValueError(".txt parameter file in the wrong format.")

	# Main benchmarking
	d_results_compute={}
	N_iter=100
	d_L_index={3:[1,2,3,4,9,10,17,18],7:[0,1,2,3,4,5,7,8,12,13]}
	for l_B in d_L_params:
		print("====================================")
		print("Benchmarking parameters with l_B={}.".format(l_B))
		print("====================================\n")
		for i in d_L_index[l_B]:
			print("Index i={}.".format(i))
			t0=time()
			L_time_compute, L_time_eval=benchmark_kani_endo(N_iter,l_B,*d_L_params[l_B][i])
			L_time_compute_half, L_time_eval_half=benchmark_kani_endo_half(N_iter,l_B,*d_L_params[l_B][i])
			L_time_compute_one, L_time_eval_one=benchmark_dim_one(N_iter,l_B,*d_L_params[l_B][i])
			
			d_results_compute[(l_B,i,"full","compute")]=L_time_compute
			d_results_compute[(l_B,i,"full","eval")]=L_time_eval
			d_results_compute[(l_B,i,"half","compute")]=L_time_compute_half
			d_results_compute[(l_B,i,"half","eval")]=L_time_eval_half
			d_results_compute[(l_B,i,"one","compute")]=L_time_compute_one
			d_results_compute[(l_B,i,"one","eval")]=L_time_eval_one
			t1=time()
			print("Time {} s.".format(t1-t0))

	with open("Benchmarking_results.csv",'w',encoding='utf-8') as f:
		for i in range(4):
			line=""
			for r in d_results_compute:
				line+=str(r[i])+","
			f.write(line+"\n")
		for i in range(N_iter):
			line=""
			for r in d_results_compute:
				line+=str(d_results_compute[r][i])+","
			f.write(line+"\n")

