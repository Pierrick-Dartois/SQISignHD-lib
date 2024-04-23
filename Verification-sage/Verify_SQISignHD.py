from sage.all import *
from time import time

from montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only
from utilities.discrete_log import weil_pairing_pari, discrete_log_pari
from utilities.order import has_order_D
from utilities.supersingular import compute_point_order_D
from theta_structures.Tuple_point import TuplePoint
from isogenies.Kani_endomorphism import KaniEndoHalf
from utilities.strategy import precompute_strategy_with_first_eval

# CLI imports
import argparse
from pathlib import Path

# Public parameters

f=126# Power of 2
fp=78# Power of 3
c=13
p=c*2**f*3**fp-1
e=142
Q0=BinaryQF([1,0,1])

Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
i=Fp2.gen()

def to_Fp2(char):
	S = "".join(char.split()) # remove white spaces
	L=S.split('+')
	return int(L[0],16)+i*int(L[1][2:],16)

def parse_sign_data(row):
	# To parse a row in the sample signature file "SQISignHD_data/SQISignHD_executions.txt"
	S = "".join(row.split()) # remove white spaces
	L=S.split(',')
	d_sign_str={}
	for s in L:
		p = s.split('=')
		# print(p)
		d_sign_str[p[0]]=p[1]

	d_sign={}
	d_sign['q']=int(d_sign_str['q'])
	d_sign['e']=int(d_sign_str['e'])
	d_sign['A_EA']=to_Fp2(d_sign_str['A_EA'])
	d_sign['xP3A']=to_Fp2(d_sign_str['xP3A'])
	d_sign['xQ3A']=to_Fp2(d_sign_str['xQ3A'])
	d_sign['xP3AmQ3A']=to_Fp2(d_sign_str['xP3AmQ3A'])
	d_sign['ker_phi_vect']=[int(d_sign_str['ker_phi_vect[0]']),int(d_sign_str['ker_phi_vect[1]'])]
	d_sign['xPA']=to_Fp2(d_sign_str['xPA'])
	d_sign['xQA']=to_Fp2(d_sign_str['xQA'])
	d_sign['xPAmQA']=to_Fp2(d_sign_str['xPAmQA'])
	d_sign['A_E1']=to_Fp2(d_sign_str['A_E1'])
	d_sign['xP1']=to_Fp2(d_sign_str['xP1'])
	d_sign['xQ1']=to_Fp2(d_sign_str['xQ1'])
	d_sign['xP1mQ1']=to_Fp2(d_sign_str['xP1mQ1'])
	d_sign['j_EA']=to_Fp2(d_sign_str['j_EA'])
	d_sign['j_E1']=to_Fp2(d_sign_str['j_E1'])
	d_sign['M_sigma']=matrix(ZZ,[[int(d_sign_str['M_sigma[00]']),int(d_sign_str['M_sigma[01]'])],[int(d_sign_str['M_sigma[10]']),int(d_sign_str['M_sigma[11]'])]])
	return d_sign

def parse_one_signature(file_open):
	# To parse a new generated signature
	with open(file_open,'r',encoding='utf-8') as f:
		L=f.readlines()
		row_1=L[26]
		row_2=L[29]
		L_row_1=row_1.split(",")
		row=L_row_1[0]+","+L_row_1[1]+","+row_2
	return parse_sign_data(row)

def read_exec(file_open):
	L_ret=[]
	with open(file_open,'r',encoding='utf-8') as f:
		L=f.readlines()
		for x in L:
			L_ret.append(parse_sign_data(x))
	return L_ret



def j_inv_montgomery(A):
	return 256*(A**2-3)**3/(A**2-4)

def is_on_curve(A,xP):
	return is_square(xP*(xP**2+A*xP+1))

def x_to_points(E,xP,xQ,xPmQ):
	a_inv=E.a_invariants()
	
	A = a_inv[1]
	if a_inv != (0,A,0,1,0):
		raise ValueError("The elliptic curve E is not in the Montgomery model.")

	yP2=xP*(xP**2+A*xP+1)
	if not is_square(yP2):
		raise ValueError("The Montgomery point is not on the base field.")
	else:
		yP=yP2.sqrt()
		P=E([xP,yP])

	yQ2=xQ*(xQ**2+A*xQ+1)
	if not is_square(yQ2):
		raise ValueError("The Montgomery point is not on the base field.")
	else:
		yQ=yQ2.sqrt()
		Q=E([xQ,yQ])

	PmQ=P-Q
	if PmQ[0]==xPmQ:
		return (P,Q)
	else:
		Q=-Q
		PmQ=P-Q
		if PmQ[0]==xPmQ:
			return (P,Q)
		else:
			raise ValueError("xPmQ inconsistent with xP and xQ.")

# Deprecated
def has_order(l,e,P):
	E=P.curve()
	O=E(0)
	Q=l**(e-1)*P
	return (Q!=O)&(l*Q==O)

# Deprecated
def has_order_fact(L_fact,P):
	for x in L_fact:
		Q=P
		for y in L_fact:
			if y!=x:
				Q=(y[0]**y[1])*Q
		if not has_order(y[0],y[1],Q):
			return False
	return True

def is_sum_two_squares(n):
	L_fact=list(factor(n))
	for x in L_fact:
		if x[0]%4!=1:
			return False
	return True


def verify(d_sign):
	t0=time()
	q=d_sign['q']
	e = d_sign['e']
	assert is_prime(2**e-q)&((2**e-q)%8==5)
	#print("Is q 2^e-good?\n{}".format(is_prime(2**e-q)&((2**e-q)%8==5)))# Needs congruent to 5 mod 8 to facilitate verification

	j_EA=j_inv_montgomery(d_sign['A_EA'])
	assert j_EA==d_sign['j_EA']
	#print("Is the public key j-invariant correct?\n{}".format(j_EA==d_sign['j_EA']))

	j_E1=j_inv_montgomery(d_sign['A_E1'])
	assert j_E1==d_sign['j_E1']
	#print("Is the commitment j-invariant correct?\n{}".format(j_E1==d_sign['j_E1']))

	EA=EllipticCurve(Fp2,[0,d_sign['A_EA'],0,1,0])

	E1=EllipticCurve(Fp2,[0,d_sign['A_E1'],0,1,0])

	P3A,Q3A=x_to_points(EA,d_sign['xP3A'],d_sign['xQ3A'],d_sign['xP3AmQ3A'])
	PA,QA=x_to_points(EA,d_sign['xPA'],d_sign['xQA'],d_sign['xPAmQA'])
	P1,Q1=x_to_points(E1,d_sign['xP1'],d_sign['xQ1'],d_sign['xP1mQ1'])

	two_power=ZZ(2**f)
	three_power=ZZ(3**fp)

	for point in [P3A,Q3A]:
		assert has_order_D(point,three_power)
		#print("Is point {} of order {}?\n{}".format(point,three_power,has_order_D(point,three_power)))
	for point in [PA,QA,P1,Q1]:
		assert has_order_D(point,two_power)
		#print("Is point {} of order {}?\n{}".format(point,two_power,has_order_D(point,two_power)))

	t1=time()

	print("Signature data parsing time {} s".format(t1-t0))


	phi,E2=isogeny_from_scalar_x_only(EA,three_power, d_sign['ker_phi_vect'][1], basis=(P3A,Q3A))

	P2,Q2=evaluate_isogeny_x_only(phi, PA, QA, two_power, three_power)
	P20=P2

	R1=d_sign['M_sigma'][0,0]*P1+d_sign['M_sigma'][1,0]*Q1
	S1=d_sign['M_sigma'][0,1]*P1+d_sign['M_sigma'][1,1]*Q1
	R10=R1

	#print("Does the basis (R1,S1) have correct Weil pairing?\n{}".format(weil_pairing_pari(P2,Q2,two_power)**q==weil_pairing_pari(R1,S1,two_power)))

	assert weil_pairing_pari(P2,Q2,two_power)**q==weil_pairing_pari(R1,S1,two_power)

	a1,a2=Q0.solve_integer(2**e-q)

	# Minimal torsion requirement
	f1=ceil(e/2)+2

	# Correcting the torsion
	N=2**(f-f1)
	P2=N*P2
	Q2=N*Q2
	R1=N*R1
	S1=N*S1
	t2=time()
	print("Challenge computation time {} s".format(t2-t1))

	# Strategies
	m=0
	if a2%2==0:
		ai_div=a2
	else:
		ai_div=a1
	while ai_div%2==0:
		m+=1
		ai_div=ai_div//2

	e1=ceil(e/2)
	e2=e-e1

	strategy1=precompute_strategy_with_first_eval(e1,m,M=1,S=0.8,I=100)
	if e2==e1:
		strategy2=strategy1
	else:
		strategy2=precompute_strategy_with_first_eval(e2,m,M=1,S=0.8,I=100)

	t3=time()
	print("Time to compute the strategies {} s".format(t3-t2))

	F=KaniEndoHalf(P2,Q2,R1,S1,q,a1,a2,e,f1,strategy1,strategy2)
	t4=time()
	print("Endomorphism F=F2*F1 computation time {} s".format(t4-t3))

	C1=F.F1._isogenies[-1]._codomain
	C2=F.F2_dual._isogenies[-1]._codomain
	HC2=C2.hadamard()
	print("Do F1 and F2_dual have the same codomain?\n{}".format(C1.zero()==HC2.zero()))
	t5=time()
	print("Codomain matching verification time {} s".format(t5-t4))

	#Point to evaluate
	Q=compute_point_order_D(E2,two_power*three_power)
	t6=time()
	print("Time to find a point of order 2^f*3^fp {} s".format(t6-t5))

	T=TuplePoint(Q,E2(0),E1(0),E1(0))

	FT=F(T)
	t7=time()
	print("Time point evaluation {} s".format(t7-t6))

	print("Is the point evaluation correct ?\n{}".format(((FT[0]==a1*Q)|(FT[0]==-a1*Q))&((FT[1]==-a2*Q)|(FT[1]==a2*Q))&(FT[3]==E1(0))))
	t8=time()
	print("Total verification time {} s".format(t8-t1))

	return P3A,Q3A,PA,QA,P1,Q1,P2,Q2,R1,S1,q,Q,FT

## CLI (command line interface)
if __name__=="__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-vs","--verify_samples",action="store_true")
	parser.add_argument("-i","--index")
	parser.add_argument("-vo","--verify_one_signature")

	args = parser.parse_args()

	if args.verify_one_signature!=None:
		file_open = Path(args.verify_one_signature)
		data=parse_one_signature(file_open)
		print("===========================================================")
		print("Verifying one SQISignHD signature with parameters:")
		print(" - Prime characteristic p = {} * 2**{} * 3**{} - 1".format(c,f,fp))
		print(" - Length of the 4-dimensional 2-isogeny chain = {}".format(e))
		print(" - Used available torsion = 2**{}".format(ceil(e/2)+2))
		print("===========================================================\n")
		verify(data)
		print("\n")
	elif args.verify_samples:
		L_exec=read_exec("SQISignHD_data/SQISignHD_executions.txt")
		if args.index==None:
			print("===============================================================")
			print("Testing {} instances of SQISignHD verification with parameters:".format(len(L_exec)))
			print(" - Prime characteristic p = {} * 2**{} * 3**{} - 1".format(c,f,fp))
			print(" - Length of the 4-dimensional 2-isogeny chain = {}".format(e))
			print(" - Used available torsion = 2**{}".format(ceil(e/2)+2))
			print("===============================================================\n")
			k=0
			for data in L_exec:
				print("Test {}".format(k))
				verify(data)
				print("\n")
				k+=1
		else:
			i=int(args.index)
			print("============================================================")
			print("Testing SQISignHD verification of sample {} with parameters:".format(i))
			print(" - Prime characteristic p = {} * 2**{} * 3**{} - 1".format(c,f,fp))
			print(" - Length of the 4-dimensional 2-isogeny chain = {}".format(e))
			print(" - Used available torsion = 2**{}".format(ceil(e/2)+2))
			print("============================================================\n")
			verify(L_exec[i])

