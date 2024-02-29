from sage.all import *

def find_prime(e_A,e_B):
	f_A=e_A+2
	f_B=e_B
	f=1
	while True:
		p=f*2**f_A*3**f_B-1
		if is_prime(p):
			while f%2==0:
				f=f//2
				f_A+=1
			while f%3==0:
				f=f//3
				f_B+=1
			return f,f_A,f_B,p
		else:
			f+=1

def find_param(e_A):
	A=2**e_A
	e_B_min=floor(e_A*log(2.)/(2*log(3.)))
	e_B_max=floor(e_A*log(2.)/log(3.))
	Q0=BinaryQF([1,0,1])
	L_params=[]
	for e_B in range(e_B_min,e_B_max+1):
		B=3**e_B
		C=A-B
		try:
			x,y=Q0.solve_integer(C)
			f,f_A,f_B,p=find_prime(e_A,e_B)
			L_params.append([e_A,e_B,x,y,f,f_A,f_B,p])
		except:
			continue
	return L_params

def find_prime_gen(l_A,l_B,e_A,e_B):
	f_A,f_B=e_A,e_B
	f=4
	while True:
		p=f*l_A**f_A*l_B**f_B-1
		if is_prime(p):
			while f%l_A==0:
				f=f//l_A
				f_A+=1
			while f%l_B==0:
				f=f//l_B
				f_B+=1
			return f,f_A,f_B,p
		else:
			f=f+4

def find_param_gen(l_A,l_B,e_A):
	A=l_A**e_A
	e_B_min=floor(e_A*log(l_A)/(2*log(l_B)))
	e_B_max=floor(e_A*log(l_A)/log(l_B))
	Q0=BinaryQF([1,0,1])
	L_params=[]
	for e_B in range(e_B_min,e_B_max+1):
		B=l_B**e_B
		C=A-B
		try:
			x,y=Q0.solve_integer(C)
			if x%2==0:
				x,y=y,x
			m=0
			z=y
			while z%2==0:
				z=z//2
				m+=1
			f,f_A,f_B,p=find_prime_gen(l_A,l_B,e_A,e_B)
			L_params.append([e_A,e_B,x,y,f,f_A,f_B,p,m])
		except:
			continue
	return L_params



def generate_params(*args):
	L_params=[]
	for e_A in args:
		L_params+=find_param(e_A)
	return L_params

def save_params(L_params,file_save):
	with open(file_save,'w',encoding='utf-8') as f:
		for x in L_params:
			f.write(str(x)+"\n")

def parse_row(row):
	L=row.split(',')
	L_ret=[int(L[0][1:])]
	n=len(L)
	for i in range(1,n-1):
		L_ret.append(int(L[i]))
	L_ret.append(int(L[-1][:-2]))
	return L_ret

def read_params(file_open):
	L_ret=[]
	with open(file_open,'r',encoding='utf-8') as f:
		L=f.readlines()
		for x in L:
			L_ret.append(parse_row(x))
	return L_ret








