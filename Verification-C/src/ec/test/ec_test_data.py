from sage.all import GF, EllipticCurve, randint
from ec_test_data_functions import random_point, int_to_digit_t

import argparse
from pathlib import Path

def generate_arithmetic_data(f,l_A,f_A,l_B,f_B,n_words):
	p=f*l_A**f_A*l_B**f_B-1
	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	
	# Elliptic curve y**2=x**3+x
	# A=0 and C=1 with the convention y**2=x**3+A*x**2+C*x
	E0=EllipticCurve(Fp2,[1,0])

	P = random_point(E0)
	Q = random_point(E0)
	twoP = 2*P
	PpQ = P+Q
	PmQ = P-Q

	m = randint(0,p-1)
	mP = m*P
	PpmQ=P+m*Q

	xP = P[0]
	xQ = Q[0]
	xtwoP = twoP[0]
	xPpQ = PpQ[0]
	xPmQ = PmQ[0]
	xmP = mP[0]
	xPpmQ = PpmQ[0]

	# List of outputs in hex format (each item is a list of hex of length 64 representing output.re or output.im, 
	# except for m which is not in Fp2 but in ZZ)
	# A, C
	L_data=[int_to_digit_t(0,n_words),int_to_digit_t(0,n_words),int_to_digit_t(1,n_words),int_to_digit_t(0,n_words),int_to_digit_t(m,n_words)]

	for x in [xP,xQ,xtwoP,xPpQ,xPmQ,xmP,xPpmQ]:
		L_data.append(int_to_digit_t(int(x[0]),n_words))
		L_data.append(int_to_digit_t(int(x[1]),n_words))

	return L_data

def save_data(L,filename):
	with open(filename,'w',encoding='utf-8') as f:
		for x in L:
			f.write(str(x)+"\n")


# CLI used in the Makefile
if __name__=="__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-arith",action="store_true")
	parser.add_argument("path")
	parser.add_argument("-f")
	parser.add_argument("-l_A")
	parser.add_argument("-f_A")
	parser.add_argument("-l_B")
	parser.add_argument("-f_B")
	parser.add_argument("-n_words")
	
	args = parser.parse_args()
	
	f=int(args.f)
	f_A=int(args.f_A)
	l_B=int(args.l_B)
	f_B=int(args.f_B)
	n_words=int(args.n_words)

	if args.arith:
		if args.l_A==None:
			L_arith_data=generate_arithmetic_data(f,2,f_A,l_B,f_B,n_words)
			save_data(L_arith_data,args.path)
		else:
			l_A=int(args.l_A)
			L_arith_data=generate_arithmetic_data(f,l_A,f_A,l_B,f_B,n_words)
			save_data(L_arith_data,args.path)





