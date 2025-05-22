from sage.all import *
proof.all(False) 

import linecache
#from Cryptodome.Hash import SHAKE256
from time import time

from Theta_dim4.Theta_dim4_sage.pkg.utilities.basis_from_hints import torsion_basis_2f_from_hint
from Theta_dim4.Theta_dim4_sage.pkg.utilities.discrete_log import weil_pairing_pari, discrete_log_pari
from Theta_dim4.Theta_dim4_sage.pkg.montgomery_isogenies.isogenies_x_only import isogeny_from_scalar_x_only, evaluate_isogeny_x_only
from Theta_dim4.Theta_dim4_sage.pkg.theta_structures.Tuple_point import TuplePoint
from Theta_dim4.Theta_dim4_sage.pkg.isogenies.Kani_endomorphism import KaniEndoHalf
from Theta_dim4.Theta_dim4_sage.pkg.utilities.strategy import precompute_strategy_with_first_eval

# CLI imports
import argparse
from pathlib import Path

public_params = {1:{'p':5*2**248-1,'c':5,'e':248},3:{'p':65*2**376-1,'c':65,'e':376},5:{'p':27*2**500-1,'c':27,'e':500}}


class SQIsignHD:
	def __init__(self,lvl):
		if lvl not in [1,3,5]:
			raise ValueError("Level (parameter -lvl) should be 1, 3 or 5.")

		self.p = public_params[lvl]['p']
		self.c = public_params[lvl]['c']
		self.e = public_params[lvl]['e']
		self.Fp2 = GF(self.p**2,'i',modulus=[1,0,1],proof=False)
		self.i = self.Fp2.gen()

		self.NQR_TABLE = self.read_gf_table("Data/NQR_TABLE_lvl"+str(lvl)+".txt")
		self.Z_NQR_TABLE = self.read_gf_table("Data/Z_NQR_TABLE_lvl"+str(lvl)+".txt")

		for x in self.NQR_TABLE:
			assert(x.is_square() == False)
		for x in self.Z_NQR_TABLE:
			assert(x.is_square() and (x-1).is_square() == False)

		self.lvl = lvl
		self.lamb = 96+32*lvl
		self.n_bytes = 2*self.lamb//8
		self.f = self.lamb + ceil(log(2*self.lamb)/log(2))
		self.r = ceil(self.f/2)+2
		self.Q0 = BinaryQF([1,0,1])


	def read_gf_table(self,file):
		L=[]
		for j in range(20):
			line = linecache.getline(file,j+1)
			M = line.split('+')
			x = int(M[0][:-1],16)+self.i*int(M[1][3:],16)
			L.append(x)
		return L

	def read_public_key(self,file,number):
		# Starting from 0
		n_line=number*3+1
	
		# A_pk
		line = linecache.getline(file,n_line)
		L=line.split('=')
		L=L[1][1:-1].split('+')
		A_pk=int(L[0][:-1],16)+self.i*int(L[1][3:],16)

		# hints
		line = linecache.getline(file,n_line+1)
		L=line.split('=')
		hint_pk_P = int(L[1][1:-1])

		line = linecache.getline(file,n_line+2)
		L=line.split('=')
		hint_pk_Q = int(L[1][1:-1])

		return A_pk,hint_pk_P,hint_pk_Q

	def read_signature(self,file,number):
		# number starting from 0
		n_line=number*11+1

		# A_com
		line = linecache.getline(file,n_line)
		L=line.split('=')
		L=L[1][1:-1].split('+')
		A_com=int(L[0][:-1],16)+self.i*int(L[1][3:],16)

		# a
		line = linecache.getline(file,n_line+1)
		L=line.split('=')
		a = int('0x'+L[1][1:-1],16)

		# b
		line = linecache.getline(file,n_line+2)
		L=line.split('=')
		b = int('0x'+L[1][1:-1],16)

		# c_or_d
		line = linecache.getline(file,n_line+3)
		L=line.split('=')
		c_or_d = int('0x'+L[1][1:-1],16)

		# q
		line = linecache.getline(file,n_line+4)
		L=line.split('=')
		q = int('0x'+L[1][1:-1],16)

		# h_com_P
		line = linecache.getline(file,n_line+5)
		L=line.split('=')
		hint_com_P = int(L[1][1:-1])

		# h_com_Q
		line = linecache.getline(file,n_line+6)
		L=line.split('=')
		hint_com_Q = int(L[1][1:-1])

		# h_chal_P
		line = linecache.getline(file,n_line+7)
		L=line.split('=')
		hint_chal_P = int(L[1][1:-1])

		# h_chal_Q
		line = linecache.getline(file,n_line+8)
		L=line.split('=')
		hint_chal_Q = int(L[1][1:-1])

		# vec_chal
		vec_chal = [0,0]
		line = linecache.getline(file,n_line+9)
		L=line.split('=')
		vec_chal[0] = int('0x'+L[1][1:-1],16)
		line = linecache.getline(file,n_line+10)
		L=line.split('=')
		vec_chal[1] = int('0x'+L[1][1:-1],16)


		return A_com, a, b, c_or_d, q, hint_com_P, hint_com_Q, hint_chal_P, hint_chal_Q, vec_chal

class SQIsignHD_verif:
	def __init__(self,pp,number):
		self.params = pp

		file = "Data/Public_keys_lvl"+str(pp.lvl)+".txt"
		self.A_pk, self.h_pk_P, self.h_pk_Q = pp.read_public_key(file,number)

		file = "Data/Signatures_lvl"+str(pp.lvl)+".txt"
		self.A_com, self.a, self.b, self.c_or_d, self.q, self.h_com_P, self.h_com_Q,\
		self.h_chal_P, self.h_chal_Q, self.vec_chal = pp.read_signature(file,number)


	def recover_pk_and_com(self):
		self.E_pk = EllipticCurve([0,self.A_pk,0,1,0])
		self.P_pk, self.Q_pk = torsion_basis_2f_from_hint(self.E_pk,self.h_pk_P,self.h_pk_Q,self.params.NQR_TABLE,self.params.Z_NQR_TABLE)

		self.E_com = EllipticCurve([0,self.A_com,0,1,0])
		self.P_com, self.Q_com = torsion_basis_2f_from_hint(self.E_com,self.h_com_P,self.h_com_Q,self.params.NQR_TABLE,self.params.Z_NQR_TABLE)

	def recover_chal(self):
		rescale = ZZ(2**(self.params.e-self.params.lamb))
		deg = ZZ(2**(self.params.lamb))

		B_pk_lamb = (rescale*self.P_pk, rescale*self.Q_pk)
		phi_chal, self.E_chal = isogeny_from_scalar_x_only(self.E_pk, deg, self.vec_chal, B_pk_lamb)

		self.P_chal, self.Q_chal = torsion_basis_2f_from_hint(self.E_chal,self.h_chal_P,self.h_chal_Q,self.params.NQR_TABLE,self.params.Z_NQR_TABLE)

	def image_response(self):
		rescale = ZZ(2**(self.params.e-self.params.r))
		R_com = rescale*self.P_com
		S_com = rescale*self.Q_com
		R_chal = rescale*self.P_chal
		S_chal = rescale*self.Q_chal
		order = ZZ(2**self.params.r)

		w_com = weil_pairing_pari(R_com, S_com, order)
		w_chal = weil_pairing_pari(R_chal, S_chal, order)

		k = discrete_log_pari(w_com, w_chal, order)

		if self.a%2==1:
			self.c = self.c_or_d
			self.d = (inverse_mod(ZZ(self.a),order)*(k*self.q+self.b*self.c))%order
		else:
			self.d = self.c_or_d
			self.c = (inverse_mod(ZZ(self.b),order)*(self.a*self.d-k*self.q))%order

		self.R_com = R_com
		self.S_com = S_com
		self.phi_rsp_R_com = self.a*R_chal + self.b*S_chal
		self.phi_rsp_S_com = self.c*R_chal + self.d*S_chal

	def compute_HD(self):
		N = 2**self.params.f-self.q
		assert is_prime(N)&(N%4==1)
		a1, a2 = self.params.Q0.solve_integer(N)

		# Strategies
		m=0
		if a2%2==0:
			ai_div=a2
		else:
			ai_div=a1
		while ai_div%2==0:
			m+=1
			ai_div=ai_div//2

		e1=ceil(self.params.f/2)
		e2=self.params.f-e1

		strategy1=precompute_strategy_with_first_eval(e1,m,M=1,S=0.8,I=100)
		if e2==e1:
			strategy2=strategy1
		else:
			strategy2=precompute_strategy_with_first_eval(e2,m,M=1,S=0.8,I=100)

		self.F = KaniEndoHalf(self.R_com,self.S_com,self.phi_rsp_R_com,self.phi_rsp_S_com,
			self.q,a1,a2,self.params.f,self.params.r,strategy1,strategy2)

		self.a1, self.a2 = a1, a2

	def verify_middle_codomain(self):
		C1=self.F.F1._isogenies[-1]._codomain
		C2=self.F.F2_dual._isogenies[-1]._codomain
		HC2=C2.hadamard()

		return C1.zero()==HC2.zero()

	def verify_HD_image(self):
		T=TuplePoint(self.P_com,self.E_com(0),self.E_chal(0),self.E_chal(0))
		FT=self.F(T)

		a1P_com = self.a1*self.P_com
		a2P_com = self.a2*self.P_com

		return ((FT[0]==a1P_com)|(FT[0]==-a1P_com))&((FT[1]==a2P_com)|(FT[1]==-a2P_com))&(FT[3]==self.E_chal(0))

	def verify(self,verbose=True):
		if verbose:
			print("\nStarting verification:")
			t0 = time()
		self.recover_pk_and_com()
		if verbose:
			t1 = time()
			print("Commitment and public key basis generation time {} s".format(t1-t0))
			t1 = time()
		self.recover_chal()
		if verbose:
			t2 = time()
			print("Challenge and challenge basis computation time {} s".format(t2-t1))
			t2 = time()
		self.image_response()
		if verbose:
			t3 = time()
			print("Response image recovery time {} s".format(t3-t2))
			t3 = time()
		self.compute_HD()
		if verbose:
			t4 = time()
			print("4-dimensional isogeny computation time {} s".format(t4-t3))
			t5 = time()
		matching_codomain = self.verify_middle_codomain()
		if verbose:
			t6 = time()
			print("Do F1 and F2_dual have the same codomain?\n{}".format(matching_codomain))
			print("Codomain matching verification time {} s".format(t6-t5))
			t6 = time()
		correct_image = self.verify_HD_image()
		if verbose:
			t7 = time()
			print("Is the point evaluation correct ?\n{}".format(correct_image))
			print("Time image verification {} s".format(t7-t6))
			print("\nEnd verification. Total verification time {} s\n".format(t7-t0))

		return matching_codomain&correct_image

if __name__=="__main__":
	parser = argparse.ArgumentParser()

	parser.add_argument("-lvl","--level")
	parser.add_argument("-test","--test",action="store_true")
	parser.add_argument("-i","--index")
	parser.add_argument("-n","--n_samples")
	parser.add_argument("-bench","--benchmark",action="store_true")

	args = parser.parse_args()

	lvl = int(args.level)
	pp = SQIsignHD(lvl)


	if args.test:
		print("##########################################")
		print("# Testing level {} SQIsignHD verification #".format(lvl))
		print("# p = {}*2**{} - 1                       #".format(pp.c,pp.e))
		print("##########################################\n")
		if args.index!=None:
			i = int(args.index)
			n_samples = 1
			if i>=100:
				raise ValueError("Index (parameter -i) should be <=99.")
			print("Verifying response number {}".format(i))
			verif = SQIsignHD_verif(pp,i)
			t0 = time()
			test = verif.verify(verbose=True)
			t1 = time()
			dt = t1-t0
		else:
			if args.n_samples!=None:
				n_samples = int(args.n_samples)
				if n_samples>100:
					raise ValueError("The number of samples (parameter -n) should be <=100.")
			else:
				n_samples = 100

			test = True
			dt = 0
			for i in range(n_samples):
				print("Verifying response number {}".format(i))
				verif = SQIsignHD_verif(pp,i)
				t0 = time()
				test = test&verif.verify(verbose=True)
				t1 = time()
				dt += t1-t0

	if args.benchmark:
		print("###############################################")
		print("# Benchmarking level {} SQIsignHD verification #".format(lvl))
		print("# p = {}*2**{} - 1                            #".format(pp.c,pp.e))
		print("###############################################\n")

		if args.n_samples!=None:
			n_samples = int(args.n_samples)
			if n_samples>100:
				raise ValueError("The number of samples (parameter -n) should be <=100.")
		else:
			n_samples = 100

		test = True
		dt = 0
		for i in range(n_samples):
			print("Verifying response number {}".format(i))
			verif = SQIsignHD_verif(pp,i)
			t0 = time()
			test = test&verif.verify(verbose=False)
			t1 = time()
			dt += t1-t0

	if test:
		print("\nAll tests passed.")
		print("\nAverage verification time {} s".format(dt/n_samples))
	else:
		print("\nError: some tests failed.")
		print("\nAverage verification time {} s".format(dt/n_samples))
		
