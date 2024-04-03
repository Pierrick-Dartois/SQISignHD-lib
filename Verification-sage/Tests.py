from sage.all import *
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite
from time import time

from parameters.parameter_generation import read_params
from utilities.supersingular import random_point, compute_point_order_D, torsion_basis
from isogenies.Kani_endomorphism import KaniEndo, KaniEndoHalf
from theta_structures.Tuple_point import TuplePoint
from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only
from basis_change.canonical_basis_dim1 import make_canonical
from utilities.strategy import precompute_strategy_with_first_eval, precompute_strategy_with_first_eval_and_splitting


filename_3="parameters/parameters.txt"
L_params_3=read_params(filename_3)

filename_7="parameters/parameters_7.txt"
L_params_7=read_params(filename_7)



def random_walk(E0,N):
	P0,Q0=torsion_basis(E0,N)
	P0,Q0,_,_,_=make_canonical(P0,Q0,N)# Q0 is above (0,0) which should not be in the kernel

	lamb=randint(1,N-1)

	_, E1 = isogeny_from_scalar_x_only(E0, N, lamb, basis=(P0,Q0))
	return E1

def test_kani_endomorphism(index,l_B=7):
	r""" Computes dimension 4 2-isogeny chains derived from Kani's lemma when the full torsion
	is available.

	INPUT:
	- index: specifies the index of the set of parameters in the list extracted from the file "parameters/parameters.txt"
	if l_B=3 or "parameters/parameters_7.txt" if l_B=7.
	- l_B: prime specifying the degree of the embedded isogeny sigma (deg(sigma)=l_B**e_B), l_B=3 or 7 (7 by default).

	OUTPUT: an object of the class KaniEndo representing a dimension 4 2-isogeny chain derived from Kani's lemma.
	"""

	t0=time()

	if l_B==7:
		e_A,e_B,a1,a2,f,f_A,f_B,p,m=L_params_7[index]
	elif l_B==3:
		e_A,e_B,a1,a2,f,f_A,f_B,p=L_params_3[index]
		m=1
	else:
		raise ValueError("Last parameter l_B should be 3 or 7.")

	print("Testing KaniEndo with parameters:")
	print(" - Prime characteristic p = {} * 2**{} * {}**{} - 1".format(f,f_A,l_B,f_B))
	print(" - Degree of the embedded isogeny sigma q = {}**{}".format(l_B,e_B))
	print(" - a1 = {}".format(a1))
	print(" - a2 = {}".format(a2))
	print(" - Length of the dimension 4 2-isogeny = {}".format(e_A))

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	i=Fp2.gen()

	q=l_B**e_B

	E0=EllipticCurve(Fp2,[1,0])

	t1=time()
	print("Parameter import: {} s".format(t1-t0))

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

	strategy=precompute_strategy_with_first_eval_and_splitting(e_A,m,M=1,S=0.8,I=100)

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

	print("Time evaluation: {} s".format(t7-t6))

	return F

def test_kani_endomorphism_half(index,l_B=7):
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

	if l_B==7:
		e_A,e_B,a1,a2,f,f_A,f_B,p,m=L_params_7[index]
	elif l_B==3:
		e_A,e_B,a1,a2,f,f_A,f_B,p=L_params_3[index]
		m=1
	else:
		raise ValueError("Last parameter l_B should be 3 or 7.")

	print("Testing KaniEndoHalf with parameters:")
	print(" - Prime characteristic p = {} * 2**{} * {}**{} - 1".format(f,f_A,l_B,f_B))
	print(" - Degree of the embedded isogeny sigma q = {}**{}".format(l_B,e_B))
	print(" - a1 = {}".format(a1))
	print(" - a2 = {}".format(a2))
	print(" - Length of the dimension 4 2-isogeny = {}".format(e_A))
	print(" - Used available torsion = 2**{}".format(ceil(e_A/2)+2))

	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)

	i=Fp2.gen()

	q=l_B**e_B

	E0=EllipticCurve(Fp2,[1,0])

	t1=time()
	print("Parameter import: {} s".format(t1-t0))

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
	print("Time evaluation: {} s".format(t7-t6))

	return F

test_endomorphism_3=False
if test_endomorphism_3:
	print("===========================================================")
	print("Testing Kani endomorphism computation (class KaniEndo) when\nthe embedded isogeny has degree deg(sigma) = 3**{*}.")
	print("===========================================================\n")
	for i in range(1,len(L_params_3)):
		F=test_kani_endomorphism(i,3)
		print("\n")

test_endomorphism_7=False
if test_endomorphism_7:
	print("===========================================================")
	print("Testing Kani endomorphism computation (class KaniEndo) when\nthe embedded isogeny has degree deg(sigma) = 7**{*}.")
	print("===========================================================\n")
	for i in range(len(L_params_7)):
		F=test_kani_endomorphism(i)
		print("\n")

test_half_endomorphism_3=False
if test_half_endomorphism_3:
	print("===============================================================================")
	print("Testing Kani endomorphism with half the torsion computation (class KaniEndoHalf)\nwhen the embedded isogeny has degree deg(sigma) = 3**{*}.")
	print("===============================================================================\n")
	for i in range(1,len(L_params_3)):
		F=test_kani_endomorphism_half(i,3)
		print("\n")

test_half_endomorphism_7=False
if test_half_endomorphism_7:
	print("===============================================================================")
	print("Testing Kani endomorphism with half the torsion computation (class KaniEndoHalf)\nwhen the embedded isogeny has degree deg(sigma) = 7**{*}.")
	print("===============================================================================\n")
	for i in range(len(L_params_7)):
		F=test_kani_endomorphism_half(i)
		print("\n")









