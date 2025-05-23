from sage.all import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
from time import time

from parameters.parameter_generation import read_params, find_param, find_param_gen, save_params
from utilities.supersingular import random_point, compute_point_order_D, torsion_basis
from isogenies.Kani_endomorphism import KaniEndo, KaniEndoHalf
from theta_structures.Tuple_point import TuplePoint
from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only
from basis_change.canonical_basis_dim1 import make_canonical
from utilities.strategy import precompute_strategy_with_first_eval#, precompute_strategy_with_first_eval_and_splitting

# CLI imports
import argparse
from pathlib import Path


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


def lookup_and_display_params(index,l_B):
	if l_B not in d_L_params:
		raise ValueError("No parameter with prime l_B={l_B} has been generated.")

	if l_B==3:
		e_A,e_B,a1,a2,f,f_A,f_B,p=d_L_params[l_B][index]
		m=1
	else:
		e_A,e_B,a1,a2,f,f_A,f_B,p,m=d_L_params[l_B][index]

	print(" - Index in the list of parameters = {}".format(index))
	print(" - Prime characteristic p = {} * 2**{} * {}**{} - 1".format(f,f_A,l_B,f_B))
	print(" - Degree of the embedded isogeny sigma q = {}**{}".format(l_B,e_B))
	print(" - a1 = {}".format(a1))
	print(" - a2 = {}".format(a2))
	print(" - m = max(v_2(a1),v_2(a2)) = {}".format(m))
	print(" - Length of the 4-dimensional 2-isogeny = {}".format(e_A))

	return e_A,e_B,a1,a2,f,f_A,f_B,p,m

def display_all_params(l_B,index):
	if l_B==None:
		for l in d_L_params:
			print("===========================================")
			print("All parameters with second prime l_B={}.".format(l))
			print("===========================================\n")
	
			for index in range(len(d_L_params[l])):
				lookup_and_display_params(index,l)
				print("\n")
	elif int(l_B) not in d_L_params:
		raise ValueError("No parameter with prime l_B={} has been generated.".format(l_B))
	elif index==None:
		l=int(l_B)
		print("===========================================")
		print("All parameters with second prime l_B={}.".format(l))
		print("===========================================\n")
	
		for index in range(len(d_L_params[l])):
			lookup_and_display_params(index,l)
			print("\n")
	else:
		l=int(l_B)
		i=int(index)
		lookup_and_display_params(i,l)

## Add params
def add_params(l_B,e_A):
	if l_B%2==0:
		raise ValueError("l_B should be odd.")
	elif l_B==3:
		L_params=find_param(e_A)
	else:
		L_params=find_param_gen(2,l_B,e_A)

	filename=str(target_dir)+"/parameters_"+str(l_B)+".txt"
	if filename in [str(file) for file in target_dir.iterdir()]:
		L_params_prev=read_params(filename)
	else:
		L_params_prev=[]

	print("===========================================")
	print("New parameters with second prime l_B={}.".format(l_B))
	print("===========================================\n")

	for i in range(len(L_params)):
		if l_B==3:
			e_A,e_B,a1,a2,f,f_A,f_B,p=L_params[i]
			print(" - Index in the list of parameters = {}".format(i+len(L_params_prev)))
			print(" - Prime characteristic p = {} * 2**{} * {}**{} - 1".format(f,f_A,l_B,f_B))
			print(" - Degree of the embedded isogeny sigma q = {}**{}".format(l_B,e_B))
			print(" - a1 = {}".format(a1))
			print(" - a2 = {}".format(a2))
			print(" - Length of the 4-dimensional 2-isogeny = {}".format(e_A))
			print("\n")
		else:
			e_A,e_B,a1,a2,f,f_A,f_B,p,m=L_params[i]
			print(" - Index in the list of parameters = {}".format(i+len(L_params_prev)))
			print(" - Prime characteristic p = {} * 2**{} * {}**{} - 1".format(f,f_A,l_B,f_B))
			print(" - Degree of the embedded isogeny sigma q = {}**{}".format(l_B,e_B))
			print(" - a1 = {}".format(a1))
			print(" - a2 = {}".format(a2))
			print(" - m = max(v_2(a1),v_2(a2)) = {}".format(m))
			print(" - Length of the 4-dimensional 2-isogeny = {}".format(e_A))
			print("\n")

	save_params(L_params_prev+L_params,filename)

## Testing 4-dimensional isogenies
def random_walk(E0,N):
	P0,Q0=torsion_basis(E0,N)
	P0,Q0,_,_,_=make_canonical(P0,Q0,N)# Q0 is above (0,0) which should not be in the kernel

	lamb=randint(0,N-1)

	_, E1 = isogeny_from_scalar_x_only(E0, N, lamb, basis=(P0,Q0))
	return E1

def test_kani_endomorphism(l_B,e_A,e_B,a1,a2,f,f_A,f_B,p,m=None,primality_check=True):
	r""" Computes dimension 4 2-isogeny chains derived from Kani's lemma when the full torsion
	is available.

	INPUT:
	- index: specifies the index of the set of parameters in the list extracted from the file "parameters/parameters.txt"
	if l_B=3 or "parameters/parameters_7.txt" if l_B=7.
	- l_B: prime specifying the degree of the embedded isogeny sigma (deg(sigma)=l_B**e_B), l_B=3 or 7 (7 by default).

	OUTPUT: an object of the class KaniEndo representing a dimension 4 2-isogeny chain derived from Kani's lemma.
	"""

	t0=time()

	q=l_B**e_B

	# Recovering m from a1 or a2
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

	if q+a1**2+a2**2!=2**e_A:
		raise ValueError("Wrong choice of parameters: q+a1**2+a2**2!=2**e_A.")
	elif p!=f*2**f_A*l_B**f_B-1:
		raise ValueError("Wrong choice of parameters: p!=f*2**f_A*l_B**f_B-1.")
	elif primality_check and not is_prime(p):
		raise ValueError("Wrong choice of parameters: p is not prime.")
	elif f_A<e_A+2:
		raise ValueError("Wrong choice of parameters: f_A<e_A+2.")
	elif f_B<e_B:
		raise ValueError("Wrong choice of parameters: f_B<e_B.")
	elif m!=m_comp:
		raise ValueError("Wrong choice of parameters: m!=max(v_2(a1),v_2(a2)).")

	print("Testing KaniEndo with parameters:")
	print(" - Prime characteristic p = {} * 2**{} * {}**{} - 1".format(f,f_A,l_B,f_B))
	print(" - Degree of the embedded isogeny sigma q = {}**{}".format(l_B,e_B))
	print(" - a1 = {}".format(a1))
	print(" - a2 = {}".format(a2))
	print(" - m = max(v_2(a1),v_2(a2)) = {}".format(m))
	print(" - Length of the dimension 4 2-isogeny = {}".format(e_A))

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	i=Fp2.gen()
	
	E0=EllipticCurve(Fp2,[1,0])

	t1=time()
	print("Setup: {} s".format(t1-t0))

	# Generating a random elliptic curve in the supersingular isogeny graph.
	N=ZZ(2**(f_A-1)*l_B**(f_B-1))
	E1=random_walk(E0,N)

	t2=time()
	print("Random walk: {} s".format(t2-t1))

	Basis_sigma=torsion_basis(E1,q)
	scalar=randint(0,q-1)

	# q-isogeny of kernel <Basis_sigma[0]+scalar*Basis_sigma[1]>
	sigma, E2 = isogeny_from_scalar_x_only(E1, ZZ(q), scalar, basis=Basis_sigma)

	t3=time()
	print("Generation of sigma: {} s".format(t3-t2))

	P1,Q1=torsion_basis(E1,2**f_A)
	R2,S2=evaluate_isogeny_x_only(sigma, P1, Q1, ZZ(2**f_A), ZZ(q))

	t4=time()
	print("Generation and evaluation of the torsion basis: {} s".format(t4-t3))

	strategy=precompute_strategy_with_first_eval(e_A,m,M=1,S=0.8,I=100)

	t5=time()
	print("Strategy computation: {} s".format(t5-t4))

	F=KaniEndo(P1,Q1,R2,S2,q,a1,a2,e_A,f_A,strategy)

	t6=time()
	print("Dimension 4 endomorphism: {} s".format(t6-t5))

	T=TuplePoint(P1,E1(0),E2(0),E2(0))

	try:
		FT=F(T)
	except:
		return F

	t7=time()

	print("Is evaluation correct?\n{}".format((FT[0][0]==(a1*P1)[0])&(FT[1][0]==(-a2*P1)[0])&(FT[2][0]==(-R2)[0])&(FT[3]==E2(0))))

	print("Time evaluation: {} s\n".format(t7-t6))

	return F
 

def test_kani_endomorphism_half(l_B,e_A,e_B,a1,a2,f,f_A,f_B,p,m=None,primality_check=True):
	r""" Computes dimension 4 2-isogeny chains derived from Kani's lemma when only half of the full torsion
	is available. 

	Since not enough torsion is given to compute the chain F at once, the computation is divided
	into two as specified in https://eprint.iacr.org/2023/436, Section 4.3. Namely, we compute
	two isogeny chain F1: E1^2*E2^2-->C and \tilde{F2}: E1^2*E2^2-->C such that F=F2 \circ F1.

	NB: In practice, the full torsion is available with the tested sets of parameters but only half is used in this
	function.

	INPUT:
	- index: specifies the index of the set of parameters in the list extracted from the file "parameters/parameters.txt"
	if l_B=3 or "parameters/parameters_7.txt" if l_B=7.
	- l_B: prime specifying the degree of the embedded isogeny $\sigma$ ($\deg(\sigma)=l_B^{e_B}$), l_B=3 or 7 (7 by default).

	OUTPUT: an object of the class KaniEndoHalf representing a dimension 4 2-isogeny chain derived from Kani's lemma.
	"""
	t0=time()

	q=l_B**e_B

	# Recovering m from a1 or a2
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
	
	if q+a1**2+a2**2!=2**e_A:
		raise ValueError("Wrong choice of parameters: q+a1**2+a2**2!=2**e_A.")
	elif p!=f*2**f_A*l_B**f_B-1:
		raise ValueError("Wrong choice of parameters: p!=f*2**f_A*l_B**f_B-1.")
	elif primality_check and not is_prime(p):
		raise ValueError("Wrong choice of parameters: p is not prime.")
	elif f_A<e_A+2:
		raise ValueError("Wrong choice of parameters: f_A<e_A+2.")
	elif f_B<e_B:
		raise ValueError("Wrong choice of parameters: f_B<e_B.")
	elif m!=m_comp:
		raise ValueError("Wrong choice of parameters: m!=max(v_2(a1),v_2(a2)).")

	print("Testing KaniEndoHalf with parameters:")
	print(" - Prime characteristic p = {} * 2**{} * {}**{} - 1".format(f,f_A,l_B,f_B))
	print(" - Degree of the embedded isogeny sigma q = {}**{}".format(l_B,e_B))
	print(" - a1 = {}".format(a1))
	print(" - a2 = {}".format(a2))
	print(" - m = max(v_2(a1),v_2(a2)) = {}".format(m))
	print(" - Length of the dimension 4 2-isogeny = {}".format(e_A))
	print(" - Used available torsion = 2**{}".format(ceil(e_A/2)+2))

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)

	i=Fp2.gen()

	E0=EllipticCurve(Fp2,[1,0])

	t1=time()
	print("Setup: {} s".format(t1-t0))

	N=ZZ(2**(f_A-1)*l_B**(f_B-1))
	E1=random_walk(E0,N)

	t2=time()
	print("Random walk: {} s".format(t2-t1))

	Basis_sigma=torsion_basis(E1,q)
	scalar=randint(0,q-1)

	# q-isogeny of kernel <Basis_sigma[0]+scalar*Basis_sigma[1]>
	sigma, E2 = isogeny_from_scalar_x_only(E1, ZZ(q), scalar, basis=Basis_sigma)

	t3=time()
	print("Generation of sigma: {} s".format(t3-t2))

	f=ceil(e_A/2)+2

	P1,Q1=torsion_basis(E1,2**f)
	R2,S2=evaluate_isogeny_x_only(sigma, P1, Q1, ZZ(2**f), ZZ(q))

	t4=time()
	print("Generation and evaluation of the torsion basis: {} s".format(t4-t3))

	e1=ceil(e_A/2)
	e2=e_A-e1

	strategy1=precompute_strategy_with_first_eval(e1,m,M=1,S=0.8,I=100)
	if e2==e1:
		strategy2=strategy1
	else:
		strategy2=precompute_strategy_with_first_eval(e2,m,M=1,S=0.8,I=100)

	t5=time()
	print("Computation of strategies: {} s".format(t5-t4))

	F=KaniEndoHalf(P1,Q1,R2,S2,q,a1,a2,e_A,f,strategy1,strategy2)

	t6=time()
	print("Dimension 4 endomorphism: {} s".format(t6-t5))

	T=TuplePoint(P1,E1(0),E2(0),E2(0))

	FT=F(T)

	t7=time()

	print("Is evaluation correct?\n{}".format((FT[0][0]==(a1*P1)[0])&(FT[1][0]==(-a2*P1)[0])&(FT[2][0]==(-R2)[0])&(FT[3]==E2(0))))
	print("Time evaluation: {} s\n".format(t7-t6))

	return F


## CLI (command line interface)
if __name__=="__main__":
	parser = argparse.ArgumentParser()

	# To display stored parameters
	parser.add_argument("-d","--display",action="store_true")
	# To add new parameters
	parser.add_argument("--add_params",action="store_true")
	# To run algorithms
	parser.add_argument("--no_primality_check",action="store_true")
	parser.add_argument("--KaniEndo",action="store_true")
	parser.add_argument("--KaniEndoHalf",action="store_true")
	# To set parameters 
	parser.add_argument("-l_B")
	parser.add_argument("-e_A")
	parser.add_argument("-e_B")
	parser.add_argument("-a1")
	parser.add_argument("-a2")
	parser.add_argument("-f")
	parser.add_argument("-f_A")
	parser.add_argument("-f_B")
	parser.add_argument("-p")
	parser.add_argument("-m")
	parser.add_argument("-i","--index")
	

	args = parser.parse_args()
	
	if args.display:
		display_all_params(args.l_B,args.index)
	if args.add_params:
		add_params(int(args.l_B),int(args.e_A))
	else:
		if args.KaniEndo:
			if args.l_B==None:
				for l_B in d_L_params:
					print("===========================================================")
					print("Testing Kani endomorphism computation (class KaniEndo) when\nthe embedded isogeny has degree deg(sigma) = {}**e_B.".format(l_B))
					print("===========================================================\n")

					n=len(d_L_params[l_B])
					if l_B==3:
						i_min=1
					else:
						i_min=0
					for i in range(i_min,n):
						test_kani_endomorphism(l_B,*d_L_params[l_B][i],primality_check=False)
			elif args.e_A!=None and args.e_B!=None and args.a1!=None and args.a2!=None and args.f!=None and args.f_A!=None and args.f_B!=None and args.p!=None:
				l_B,e_A,e_B,a1,a2,f,f_A,f_B,p=int(args.l_B),int(args.e_A),int(args.e_B),int(args.a1),int(args.a2),int(args.f),int(args.f_A),int(args.f_B),int(args.p)
				if args.m==None:
					m=None
				else:
					m=int(args.m)
				test_kani_endomorphism(l_B,e_A,e_B,a1,a2,f,f_A,f_B,p,m,not args.no_primality_check)
			elif args.index!=None:
				l_B=int(args.l_B)
				index=int(args.index)
				test_kani_endomorphism(l_B,*d_L_params[l_B][index],primality_check=False)
			else:
				l_B=int(args.l_B)
				print("===========================================================")
				print("Testing Kani endomorphism computation (class KaniEndo) when\nthe embedded isogeny has degree deg(sigma) = {}**e_B.".format(l_B))
				print("===========================================================\n")

				n=len(d_L_params[l_B])
				if l_B==3:
					i_min=1
				else:
					i_min=0
				for i in range(i_min,n):
					test_kani_endomorphism(l_B,*d_L_params[l_B][i],primality_check=False)
		if args.KaniEndoHalf:
			if args.l_B==None:
				for l_B in d_L_params:
					print("===============================================================================")
					print("Testing Kani endomorphism with half the torsion computation (class KaniEndoHalf)\nwhen the embedded isogeny has degree deg(sigma) = {}**e_B.".format(l_B))
					print("===============================================================================\n")

					n=len(d_L_params[l_B])
					if l_B==3:
						i_min=1
					else:
						i_min=0
					for i in range(i_min,n):
						test_kani_endomorphism_half(l_B,*d_L_params[l_B][i],primality_check=False)
			elif args.e_A!=None and args.e_B!=None and args.a1!=None and args.a2!=None and args.f!=None and args.f_A!=None and args.f_B!=None and args.p!=None:
				l_B,e_A,e_B,a1,a2,f,f_A,f_B,p=int(args.l_B),int(args.e_A),int(args.e_B),int(args.a1),int(args.a2),int(args.f),int(args.f_A),int(args.f_B),int(args.p)
				if args.m==None:
					m=None
				else:
					m=int(args.m)
				test_kani_endomorphism_half(l_B,e_A,e_B,a1,a2,f,f_A,f_B,p,m,not args.no_primality_check)
			elif args.index!=None:
				l_B=int(args.l_B)
				index=int(args.index)
				test_kani_endomorphism_half(l_B,*d_L_params[l_B][index],primality_check=False)
			else:
				l_B=int(args.l_B)
				print("===============================================================================")
				print("Testing Kani endomorphism with half the torsion computation (class KaniEndoHalf)\nwhen the embedded isogeny has degree deg(sigma) = {}**e_B.".format(l_B))
				print("===============================================================================\n")

				n=len(d_L_params[l_B])
				if l_B==3:
					i_min=1
				else:
					i_min=0
				for i in range(i_min,n):
					test_kani_endomorphism_half(l_B,*d_L_params[l_B][i],primality_check=False)
