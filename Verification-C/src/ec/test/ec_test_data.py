from sage.all import GF, EllipticCurve, randint
from ec_test_data_functions import random_point, int_to_digit_t

import argparse
from pathlib import Path

def generate_arithmetic_params(f,l_A,f_A,l_B,f_B,n_words,filename):
	p=f*l_A**f_A*l_B**f_B-1
	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)
	
	# Elliptic curve y**2=x**3+x
	# A=0 and C=1 with the convention y**2=x**3+A*x**2+C*x
	E0=EllipticCurve(Fp2,[1,0])

	P = random_point(E0)
	Q = random_point(E0)
	k = randint(0,p-1)

	d_points={'P':P, 'Q':Q, 'twoP':2*P, 'PpQ':P+Q, 'PmQ':P-Q,'kP':k*P, 'PpkQ': P+k*Q}

	with open(filename,'w',encoding='utf-8') as f:
		L_file=filename[:-2].split("/")
		f.write("#ifndef "+L_file[-1].upper()+"_H\n")
		f.write("#define "+L_file[-1].upper()+"_H\n")
		f.write("\n")
		f.write("#include \"../include/ec.h\"\n")
		f.write("\n")

		L_k=int_to_digit_t(k,n_words)
		str_k='{'
		for i in range(n_words):
			str_k+=L_k[i]
			if i<n_words-1:
				str_k+=','
			else:
				str_k+='}'
		f.write("digit_t k[NWORDS_FIELD]="+str_k+";\n")

		f.write("ec_point_t AC = {{{0x0"+",0x0"*(n_words-1)+"},{0x0"+",0x0"*(n_words-1)+"}}, {{0x1"+",0x0"*(n_words-1)+"},{0x0"+",0x0"*(n_words-1)+"}}};\n")
		f.write("\n")

		for point in d_points:
			L_point_re=int_to_digit_t(int(d_points[point][0][0]),n_words)
			L_point_im=int_to_digit_t(int(d_points[point][0][1]),n_words)
			line="ec_point_t "+point+" = {{{"
			for i in range(n_words-1):
				line+=L_point_re[i]+","
			line+=L_point_re[-1]+"},{"
			for i in range(n_words-1):
				line+=L_point_im[i]+","
			line+=L_point_im[-1]+"}}, "
			line+="{{0x1"+",0x0"*(n_words-1)+"},{0x0"+",0x0"*(n_words-1)+"}}};\n"
			f.write(line)
			f.write("\n")

		for point in d_points:
			x_point_re=int_to_digit_t(int(d_points[point][0][0]),n_words)
			x_point_im=int_to_digit_t(int(d_points[point][0][1]),n_words)
			y_point_re=int_to_digit_t(int(d_points[point][1][0]),n_words)
			y_point_im=int_to_digit_t(int(d_points[point][1][1]),n_words)
			
			line="jac_point_t jac_"+point+" = {{{"
			
			for i in range(n_words-1):
				line+=x_point_re[i]+","
			line+=x_point_re[-1]+"},{"
			for i in range(n_words-1):
				line+=x_point_im[i]+","
			line+=x_point_im[-1]+"}}, {{"

			for i in range(n_words-1):
				line+=y_point_re[i]+","
			line+=y_point_re[-1]+"},{"
			for i in range(n_words-1):
				line+=y_point_im[i]+","
			line+=y_point_im[-1]+"}}, "
			
			line+="{{0x1"+",0x0"*(n_words-1)+"},{0x0"+",0x0"*(n_words-1)+"}}};\n"
			f.write(line)
			f.write("\n")
		f.write("#endif")

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
			L_arith_data=generate_arithmetic_params(f,2,f_A,l_B,f_B,n_words,args.path)
		else:
			l_A=int(args.l_A)
			L_arith_data=generate_arithmetic_params(f,l_A,f_A,l_B,f_B,n_words,args.path)





