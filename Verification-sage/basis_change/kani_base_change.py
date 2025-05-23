from sage.all import *
from basis_change.canonical_basis_dim1 import make_canonical
from basis_change.base_change_dim2 import complete_symplectic_matrix_dim2, is_symplectic_matrix_dim2
from basis_change.base_change_dim4 import complete_symplectic_matrix_dim4, is_symplectic_matrix_dim4, bloc_decomposition
from theta_structures.Tuple_point import TuplePoint


def base_change_canonical_dim2(P1,P2,R1,R2,q,f):
	r"""

	Input:
	- P1, P2: basis of E1[2**f].
	- R1, R2: images of P1, P2 by \sigma: E1 --> E2.
	- q: degree of \sigma.
	- f: log_2(order of P1 and P2).

	Output:
	- P1_doubles: list of 2**i*P1 for i in {0,...,f-2}.
	- P2_doubles: list of 2**i*P2 for i in {0,...,f-2}.
	- R1_doubles: list of 2**i*R1 for i in {0,...,f-2}.
	- R2_doubles: list of 2**i*R2 for i in {0,...,f-2},
	- T1, T2: canonical basis of E1[4].
	- U1, U2: canonical basis of E2[4].
	- M0: base change matrix of the symplectic basis 2**(f-2)*B1 of E1*E2[4] given by:
	B1:=[[(P1,0),(0,R1)],[(P2,0),(0,lamb*R2)]]
	where lamb is the modular inverse of q mod 2**f, so that:
	e_{2**f}(P1,P2)=e_{2**f}(R1,lamb*R2).
	in the canonical symplectic basis:
	B0:=[[(T1,0),(0,U1)],[(T2,0),(0,U2)]].
	"""
	lamb=inverse_mod(q,4)

	P1_doubles=[P1]
	P2_doubles=[P2]
	R1_doubles=[R1]
	R2_doubles=[R2]

	for i in range(f-2):
		P1_doubles.append(2*P1_doubles[-1])
		P2_doubles.append(2*P2_doubles[-1])
		R1_doubles.append(2*R1_doubles[-1])
		R2_doubles.append(2*R2_doubles[-1])

	# Constructing canonical basis of E1[4] and E2[4].
	_,_,T1,T2,MT=make_canonical(P1_doubles[-1],P2_doubles[-1],4,preserve_pairing=True)
	_,_,U1,U2,MU=make_canonical(R1_doubles[-1],lamb*R2_doubles[-1],4,preserve_pairing=True)

	Z4=Integers(4)
	M0=matrix(Z4,[[MT[0,0],0,MT[1,0],0],
		[0,MU[0,0],0,MU[1,0]],
		[MT[0,1],0,MT[1,1],0],
		[0,MU[0,1],0,MU[1,1]]])

	return P1_doubles,P2_doubles,R1_doubles,R2_doubles,T1,T2,U1,U2,M0

def gluing_base_change_matrix_dim2(a1,a2,q):
	r"""Computes the symplectic base change matrix of a symplectic basis (*,B_K4) of E1*E2[4]
	given by the kernel of the dimension 2 gluing isogeny:
	B_K4=2**(f-2)[([a1]P1-[a2]P2,R1),([a1]P2+[a2]P1,R2)]
	in the basis $2**(f-2)*B1$ given by:
	B1:=[[(P1,0),(0,R1)],[(P2,0),(0,lamb*R2)]]
	where:
	- lamb is the inverse of q modulo 2**f.
	- (P1,P2) is the canonical basis of E1[2**f].
	- (R1,R2) is the image of (P1,P2) by sigma.

	Input:
	- a1, q: integers.

	Output:
	- M: symplectic base change matrix of (*,B_K4) in 2**(f-2)*B1.
	"""

	Z4=Integers(4)

	mu=inverse_mod(a1,4)

	A=matrix(Z4,[[0,mu],
		[0,0]])
	B=matrix(Z4,[[0,0],
		[-1,-ZZ(mu*a2)]])

	C=matrix(Z4,[[ZZ(a1),ZZ(a2)],
		[1,0]])
	D=matrix(Z4,[[-ZZ(a2),ZZ(a1)],
		[0,ZZ(q)]])

	#M=complete_symplectic_matrix_dim2(C, D, 4)
	M=block_matrix([[A,C],[B,D]])

	assert is_symplectic_matrix_dim2(M)

	return M

# ============================================== #
#     Functions for the class KaniClapotiIsog    #
# ============================================== #

def clapoti_cob_matrix_dim2(integers):
	gu,xu,yu,gv,xv,yv,Nbk,Nck,e,m = integers

	xu = ZZ(xu)
	xv = ZZ(xv)
	Nbk = ZZ(Nbk)
	Nck = ZZ(Nck)
	u = ZZ(gu*(xu**2+yu**2))
	v = ZZ(gv*(xv**2+yv**2))
	mu = inverse_mod(u,4)
	suv = xu*xv+yu*yv
	inv_Nbk = inverse_mod(Nbk,4)
	inv_gugvNcksuv = inverse_mod(gu*gv*Nck*suv,4)

	Z4=Integers(4)

	M=matrix(Z4,[[0,0,u*Nbk,0],
		[0,inv_Nbk*inv_gugvNcksuv,gu*suv,0],
		[-inv_Nbk*mu,0,0,gu*Nbk*u],
		[0,0,0,gu*gv*Nbk*Nck*suv]])

	assert is_symplectic_matrix_dim2(M)

	return M

def clapoti_cob_matrix_dim2_dim4(integers):
	gu,xu,yu,gv,xv,yv,Nbk,Nck,e,m = integers

	xu = ZZ(xu)
	yu = ZZ(yu)
	xv = ZZ(xv)
	yv = ZZ(yv)
	gu = ZZ(gu)
	gv = ZZ(gv)
	Nbk = ZZ(Nbk)
	Nck = ZZ(Nck)
	u = ZZ(gu*(xu**2+yu**2))
	v = ZZ(gv*(xv**2+yv**2))
	suv = xu*xv+yu*yv
	duv = xv*yu-xu*yv
	duv_2m = duv//2**m
	mu = inverse_mod(u,4)
	nu = inverse_mod(v,4)
	sigmauv = inverse_mod(suv,4)
	inv_guNbk = inverse_mod(gu*Nbk,4)
	lamb = nu*gu*gv*Nbk*suv
	mu1 = ZZ(mu*gu**2*gv*suv*Nbk*Nck*duv_2m)
	mu2 = ZZ(duv_2m*gu*sigmauv*(Nbk*u*yu+gv*xv*Nck*duv))
	mu3 = ZZ(duv_2m*gu*sigmauv*(Nbk*u*xu-gv*yv*Nck*duv))

	Z4=Integers(4)

	M=matrix(Z4,[[gu*xu,-gu*yu,0,0,0,0,mu2,mu3],
		[0,0,lamb*xv,-lamb*yv,mu1*yu,mu1*xu,0,0],
		[gu*yu,gu*xu,0,0,0,0,-mu3,mu2],
		[0,0,lamb*yv,lamb*xv,-mu1*xu,mu1*yu,0,0],
		[0,0,0,0,mu*xu,-mu*yu,0,0],
		[0,0,0,0,0,0,inv_guNbk*xv*sigmauv,-inv_guNbk*yv*sigmauv],
		[0,0,0,0,mu*yu,mu*xu,0,0],
		[0,0,0,0,0,0,inv_guNbk*yv*sigmauv,inv_guNbk*xv*sigmauv]])

	assert is_symplectic_matrix_dim4(M)

	return M

def clapoti_cob_splitting_matrix(integers):
	gu,xu,yu,gv,xv,yv,Nbk,Nck,e,m = integers

	v=ZZ(gv*(xv**2+yv**2))
	vNck=ZZ(v*Nck)
	inv_vNck=inverse_mod(vNck,4)

	Z4=Integers(4)

	M=matrix(Z4,[[0,0,0,0,-1,0,0,0],
		[0,0,0,0,0,-1,0,0],
		[0,0,vNck,0,0,0,0,0],
		[0,0,0,vNck,0,0,0,0],
		[1,0,-vNck,0,0,0,0,0],
		[0,1,0,-vNck,0,0,0,0],
		[0,0,0,0,1,0,inv_vNck,0],
		[0,0,0,0,0,1,0,inv_vNck]])

	assert is_symplectic_matrix_dim4(M)

	return M

# =============================================== #
#     Functions for the class KaniFixedDegDim2    #
# =============================================== #

def fixed_deg_gluing_matrix_Phi1(u,a,b,c,d):
	u,a,b,c,d = ZZ(u),ZZ(a),ZZ(b),ZZ(c),ZZ(d)

	mu = inverse_mod(u,4)
	inv_cmd = inverse_mod(c-d,4)

	Z4 = Integers(4)

	M = matrix(Z4,[[0,0,u,0],
		[0,inv_cmd,c+d,0],
		[-mu,0,0,(d**2-c**2)*mu],
		[0,0,0,c-d]])

	assert is_symplectic_matrix_dim2(M)

	return M

def fixed_deg_gluing_matrix_Phi2(u,a,b,c,d):
	u,a,b,c,d = ZZ(u),ZZ(a),ZZ(b),ZZ(c),ZZ(d)

	mu = inverse_mod(u,4)
	inv_cpd = inverse_mod(c+d,4)

	Z4 = Integers(4)

	M = matrix(Z4,[[0,0,u,0],
		[0,-inv_cpd,d-c,0],
		[-mu,0,0,(d**2-c**2)*mu],
		[0,0,0,-(c+d)]])

	assert is_symplectic_matrix_dim2(M)

	return M

def fixed_deg_gluing_matrix_dim4(u,a,b,c,d,m):
	u,a,b,c,d = ZZ(u),ZZ(a),ZZ(b),ZZ(c),ZZ(d)

	mu = inverse_mod(u,4)
	nu = ZZ((-mu**2)%4)
	amb_2m = ZZ((a-b)//2**m)
	apb_2m = ZZ((a+b)//2**m)
	u2pc2md2_2m = ZZ((u**2+c**2-d**2)//2**m)
	inv_cmd = inverse_mod(c-d,4)
	inv_cpd = inverse_mod(c+d,4)


	Z4 = Integers(4)

	M = matrix(Z4,[[1,0,0,0,0,0,-u2pc2md2_2m,-apb_2m*(c+d)],
		[0,0,(c+d)*(c-d)*nu,(a-b)*(c-d)*nu,0,amb_2m*(c-d),0,0],
		[0,1,0,0,0,0,amb_2m*(c-d),-u2pc2md2_2m],
		[0,0,-(a+b)*(c+d)*nu,(c+d)*(c-d)*nu,-apb_2m*(c+d),0,0,0],
		[0,0,0,0,1,0,0,0],
		[0,0,0,0,0,0,1,(a+b)*inv_cmd],
		[0,0,0,0,0,1,0,0],
		[0,0,0,0,0,0,(b-a)*inv_cpd,1]])

	assert is_symplectic_matrix_dim4(M)

	return M

def fixed_deg_gluing_matrix(u,a,b,c,d):
	r"""
	Deprecated.
	"""

	mu = inverse_mod(u,4)
	nu = (-mu**2)%4

	Z4 = Integers(4)

	M = matrix(Z4,[[0,0,0,0,ZZ(u),0,0,0],
		[0,0,0,0,0,ZZ(u),0,0],
		[0,0,ZZ(nu*(a+b)),ZZ(nu*(d-c)),ZZ(a+b),ZZ(d-c),0,0],
		[0,0,ZZ(nu*(c+d)),ZZ(nu*(a-b)),ZZ(c+d),ZZ(a-b),0,0],
		[ZZ(-mu),0,0,0,0,0,ZZ(u),0],
		[0,ZZ(-mu),0,0,0,0,0,ZZ(u)],
		[0,0,0,0,0,0,ZZ(a-b),ZZ(-c-d)],
		[0,0,0,0,0,0,ZZ(c-d),ZZ(a+b)]])

	assert is_symplectic_matrix_dim4(M)

	return M

def fixed_deg_splitting_matrix(u):

	mu = inverse_mod(u,4)

	Z4 = Integers(4)

	M = matrix(Z4,[[0,0,0,0,-1,0,0,0],
		[0,0,0,0,0,-1,0,0],
		[0,0,ZZ(-u),0,0,0,0,0],
		[0,0,0,ZZ(-u),0,0,0,0],
		[1,0,ZZ(-mu),0,0,0,0,0],
		[0,1,0,ZZ(-mu),0,0,0,0],
		[0,0,0,0,ZZ(mu),0,ZZ(mu),0],
		[0,0,0,0,0,ZZ(mu),0,ZZ(mu)]])

	assert is_symplectic_matrix_dim4(M)

	return M
	

# ========================================================== #
#     Functions for the class KaniEndo (one isogeny chain)   #
# ========================================================== #

def gluing_base_change_matrix_dim2_dim4(a1,a2,m,mua2):
	r"""Computes the symplectic base change matrix of a symplectic basis (*,B_K4) of Am*Am[4]
	given by the kernel of the dimension 4 gluing isogeny Am*Am-->B:
	
	B_K4=2**(e-m)[(Phi([a1]P1,sigma(P1)),Phi([a2]P1,0)),(Phi([a1]P2,sigma(P2)),Phi([a2]P2,0)),
	(Phi(-[a2]P1,0),Phi([a1]P1,sigma(P1))),(Phi(-[a2]P2,0),Phi([a2]P2,sigma(P2)))]
	
	in the basis associated to the product theta-structure of level 2 of Am*Am:
	
	B:=[(S1,0),(S2,0),(0,S1),(0,S2),(T1,0),(T2,0),(0,T1),(0,T2)]
	
	where:
	- (P1,P2) is the canonical basis of E1[2**f].
	- (R1,R2) is the image of (P1,P2) by sigma.
	- Phi is the 2**m-isogeny E1*E2-->Am (m first steps of the chain in dimension 2).
	- S1=[2**e]Phi([lamb]P2,[a]sigma(P1)+[b]sigma(P2)).
	- S2=[2**e]Phi([mu]P1,[c]sigma(P1)+[d]sigma(P2)).
	- T1=[2**(e-m)]Phi([a1]P1-[a2]P2,sigma(P1)).
	- T2=[2**(e-m)]Phi([a1]P2+[a2]P1,sigma(P2)).
	- (S1,S2,T1,T2) is induced by the image by Phi of a symplectic basis of E1*E2[2**(m+2)] lying
	above the symplectic basis of E1*E2[4] outputted by gluing_base_change_matrix_dim2.

	INPUT:
	- a1, a2: integers.
	- m: integer (number of steps in dimension 2).
	- mua2: product mu*a2.

	OUTPUT:
	- M: symplectic base change matrix of (*,B_K4) in B.
	"""
	a1a2_2m=ZZ(a1*a2//2**m)
	a22_2m=ZZ(a2**2//2**m)

	Z4=Integers(4)

	C=matrix(Z4,[[-a1a2_2m,a22_2m,a22_2m,a1a2_2m],
		[-a22_2m,-a1a2_2m,-a1a2_2m,a22_2m],
		[-a22_2m,-a1a2_2m,-a1a2_2m,a22_2m],
		[a1a2_2m,-a22_2m,-a22_2m,-a1a2_2m]])

	D=matrix(Z4,[[1,0,0,0],
		[mua2,1,0,-mua2],
		[0,0,1,0],
		[0,mua2,mua2,1]])

	M=complete_symplectic_matrix_dim4(C,D,4)

	assert is_symplectic_matrix_dim4(M)

	return M

def splitting_base_change_matrix_dim4(a1,a2,q,m,M0,A_B,mu=None):
	r"""
	Let F be the endomorphism of E1^2*E2^2 given by Kani's lemma. Write: 
			E1^2*E2^2 -- Phi x Phi --> Am^2 -- G --> E1^2*E2^2,
	where Phi: E1 x E1 --> Am is a 2**m-isogeny in dimension 2. 
	Let (U_1,...,U_4,V_1,...,V_4) be a symplectic basis of Am^2[2**(e-m+2)] 
	such that V_i=Phi x Phi(W_i), where W_1,...,W_4 have order 2**(e+2), lie over ker(F) 
	and generate an isotropic subgroup:
	W_1=([a1]P1-[2^e/a1]P1,[a2]P1,R2,0)
	W_2=([a1]Q1,[a2]Q1,S2,0)
	W_3=(-[a2]P1,[a1]P1,0,R2)
	W_4=(-[a2]Q1,[a1]Q1,0,S2),
	with (P1,Q1), a basis of E1[2**(e+2)] and (R2,S2) its image via 
	sigma: E1 --> E2. Then B:=([2^(e-m)]G(U_1),...,[2^(e-m)]G(U_4),G(V_1),...,G(V_4))
	is a symplectic basis of E1^2*E2^2[4].

	We assume that ([2^(e-m)]U_1,...,[2^(e-m)]U_4) is the symplectic complement of 
	([2^(e-m)]V_1,...,[2^(e-m)]V_4) that has been outputted by 
	gluing_base_change_matrix_dim2_dim4 for the gluing isogeny on Am^2 
	(first 2-isogeny of G). This function computes the base change matrix of B
	in the symplectic basis of E1^2*E2^2[4]:
	B0=[(T1,0,0,0),(0,T1,0,0),(0,0,T2,0),(0,0,0,T2),(U1,0,0,0),(0,U1,0,0),
	(0,0,U2,0),(0,0,0,U2)]
	associated to the product Theta structure on E1^2*E2^2.

	INPUT:
	- a1,a2,q: integers defining F (q=deg(sigma)).
	- m: 2-adic valuation of a2.
	- M0: base change matrix of the symplectic basis 2**e*B1 of E1*E2[4] 
	given by:
	B1:=[[(P1,0),(0,R2)],[(Q1,0),(0,lamb*S2)]]
	in the canonical symplectic basis:
	B0:=[[(T1,0),(0,T2)],[(U1,0),(0,U2)]],
	where lamb is the modular inverse of q mod 2**(e+2), so that:
	e_{2**(e+2)}(P1,P2)=e_{2**(e+2)}(R1,lamb*R2).
	- A_B: 4 first columns (left part) of the symplectic matrix outputted by
	gluing_base_change_matrix_dim2_dim4.
	- mu, a, b, c, d: integers defining the product Theta structure of Am^2
	given by the four torsion basis [2**(e-m)]*B1 of Am, where:
	B1=[[2**m]Phi([2**(m+1)]P2,[a]sigma(P1)+[b]sigma(P2)),
	[2**m]Phi([mu]P1,[2**(m+1)]sigma(P1)+[d]sigma(P2)),
	Phi([a1]P1-[a2]P2,sigma(P1)),
	Phi([a1]P2+[a2]P1,sigma(P2))].
	Only mu is given.

	OUTPUT: The desired base change matrix.
	"""
	Z4=Integers(4)

	a2_2m=ZZ(a2//2**m)
	a12_q_2m=ZZ((a1**2+q)//2**m)

	inv_q=inverse_mod(q,4)
	inv_a1=inverse_mod(a1,4)

	lamb=ZZ(2**(m+1))
	if mu==None:
		mu=ZZ((1-2**(m+1)*q)*inv_a1)
	a=ZZ(2**(m+1)*a2*inv_q)
	b=ZZ(-(1+2**(m+1)*a1)*inv_q)
	c=ZZ(2**(m+1))
	d=ZZ(-mu*a2*inv_q)

	# Matrix of the four torsion basis of E1^2*E2^2[4] given by
	# ([2^(e-m)]G(B1[0],0),[2^(e-m)]G(B1[1],0),[2^(e-m)]G(0,B1[0]),[2^(e-m)]G(0,B1[1]),
	# G(B1[2],0),G(B1[3],0),G(0,B1[2]),G(0,B1[3])) in the basis induced by 
	# [2**e](P1,Q1,R2,[1/q]S2)
	M1=matrix(Z4,[[a*q,mu*a1+c*q,0,mu*a2,a12_q_2m,a1*a2_2m,a1*a2_2m,a2*a2_2m],
		[0,-mu*a2,a*q,mu*a1+c*q,-a1*a2_2m,-a2*a2_2m,a12_q_2m,a1*a2_2m],
		[a1*a,a1*c-mu,-a*a2,-c*a2,0,-a2_2m,-a2_2m,0],
		[a2*a,a2*c,a*a1,c*a1-mu,a2_2m,0,0,-a2_2m],
		[lamb*a1+b*q,d*q,lamb*a2,0,-a1*a2_2m,a12_q_2m,-a2*a2_2m,a1*a2_2m],
		[-lamb*a2,0,lamb*a1+b*q,d*q,a2*a2_2m,-a1*a2_2m,-a1*a2_2m,a12_q_2m],
		[(a1*b-lamb)*q,a1*d*q,-b*a2*q,-a2*d*q,a2_2m*q,0,0,-a2_2m*q],
		[a2*b*q,a2*d*q,(b*a1-lamb)*q,a1*d*q,0,a2_2m*q,a2_2m*q,0]])
	#A,B,C,D=bloc_decomposition(M1)
	#if B.transpose()*A!=A.transpose()*B:
		#print("B^T*A!=A^T*B")
	#if C.transpose()*D!=D.transpose()*C:
		#print("C^T*D!=D^T*C")
	#if A.transpose()*D-B.transpose()*C!=identity_matrix(4):
		#print(A.transpose()*D-B.transpose()*C)
		#print("A^T*D-B^T*C!=I")
	#print(M1)
	#print(M1.det())

	# Matrix of ([2^e]G(U_1),...,[2^e]G(U_4)) in the basis induced by
	# [2**e](P1,Q1,R2,[1/q]S2)
	M_left=M1*A_B
	#print(A_B)
	#print(M_left)

	# Matrix of (G(V_1),...,G(V_4)) in the basis induced by [2**e](P1,Q1,R2,[1/q]S2)
	M_right=matrix(Z4,[[0,0,0,0],
		[a2*inv_a1,0,1,0],
		[inv_a1,0,0,0],
		[0,0,0,0],
		[0,1,0,-a2*inv_a1],
		[0,0,0,0],
		[0,0,0,0],
		[0,0,0,q*inv_a1]])

	# Matrix of the basis induced by [2**e](P1,Q1,R2,[1/q]S2) in the basis 
	# B0 (induced by T1, U1, T2, U2)
	MM0=matrix(Z4,[[M0[0,0],0,M0[0,1],0,M0[0,2],0,M0[0,3],0],
				  [0,M0[0,0],0,M0[0,1],0,M0[0,2],0,M0[0,3]],
				  [M0[1,0],0,M0[1,1],0,M0[1,2],0,M0[1,3],0],
				  [0,M0[1,0],0,M0[1,1],0,M0[1,2],0,M0[1,3]],
				  [M0[2,0],0,M0[2,1],0,M0[2,2],0,M0[2,3],0],
				  [0,M0[2,0],0,M0[2,1],0,M0[2,2],0,M0[2,3]],
				  [M0[3,0],0,M0[3,1],0,M0[3,2],0,M0[3,3],0],
				  [0,M0[3,0],0,M0[3,1],0,M0[3,2],0,M0[3,3]]])

	M=MM0*block_matrix(1,2,[M_left,M_right])

	#A,B,C,D=bloc_decomposition(M)

	#M=complete_symplectic_matrix_dim4(C,D)

	#print(M.det())
	#print(M)

	A,B,C,D=bloc_decomposition(M)
	if B.transpose()*A!=A.transpose()*B:
		print("B^T*A!=A^T*B")
	if C.transpose()*D!=D.transpose()*C:
		print("C^T*D!=D^T*C")
	if A.transpose()*D-B.transpose()*C!=identity_matrix(4):
		print("A^T*D-B^T*C!=I")
	assert is_symplectic_matrix_dim4(M)

	return M

# ============================================================================ #
#     Functions for the class KaniEndoHalf (isogeny chain decomposed in two)   #
# ============================================================================ #

def complete_kernel_matrix_F1(a1,a2,q,f):
	r"""Computes the symplectic base change matrix of a symplectic basis of the form (*,B_Kp1) 
	in the symplectic basis of E1^2*E2^2[2**f] given by:
	B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1)],
	[(P2,0,0,0),(0,P2,0,0),(0,0,lamb*R2,0),(0,0,0,lamb*R2)]]
	where:
	- B_Kp1 is a basis of an isotropic subgroup of E1^2*E2^2[2**f] lying above ker(F1).
	By convention B_Kp1=[(\tilde{\alpha}_1(P1,0),\Sigma(P1,0)),
	(\tilde{\alpha}_1(P2,0),\Sigma(P2,0)),
	(\tilde{\alpha}_1(0,P1),\Sigma(0,P1)),
	(\tilde{\alpha}_1(0,P2),\Sigma(0,P2))]
	- lamb is the inverse of q modulo 2**f.
	- (P1,P2) is the canonical basis of E1[2**f].
	- (R1,R2) is the image of (P1,P2) by sigma.

	Input:
	- a1, a2, q: Integers such that q+a1**2+a2**2=2**e.
	- f: integer determining the accessible 2-torsion in E1 (E1[2**f]).

	Output:
	- M: symplectic base change matrix of (*,B_Kp1) in B1.
	"""
	N=2**f
	ZN=Integers(N)

	C=matrix(ZN,[[a1,0,-a2,0],
		[a2,0,a1,0],
		[1,0,0,0],
		[0,0,1,0]])

	D=matrix(ZN,[[0,a1,0,-a2],
		[0,a2,0,a1],
		[0,q,0,0],
		[0,0,0,q]])

	assert C.transpose()*D==D.transpose()*C

	M=complete_symplectic_matrix_dim4(C,D,N)

	assert is_symplectic_matrix_dim4(M)

	return M

def complete_kernel_matrix_F2_dual(a1,a2,q,f):
	r"""Computes the symplectic base change matrix of a symplectic basis of the form (*,B_Kp2) 
	in the symplectic basis of E1^2*E2^2[2**f] given by:
	B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1)],
	[(P2,0,0,0),(0,P2,0,0),(0,0,lamb*R2,0),(0,0,0,lamb*R2)]]
	where:
	- B_Kp2 is a basis of an isotropic subgroup of E1^2*E2^2[2**f] lying above ker(\tilde{F2}).
	By convention B_Kp2=[(\alpha_1(P1,0),-\Sigma(P1,0)),
	(\alpha_1(P2,0),-\Sigma(P2,0)),
	(\alpha_1(0,P1),-\Sigma(0,P1)),
	(\alpha_1(0,P2),-\Sigma(0,P2))]. 
	- lamb is the inverse of q modulo 2**f.
	- (P1,P2) is the canonical basis of E1[2**f].
	- (R1,R2) is the image of (P1,P2) by sigma.

	Input:
	- a1, a2, q: Integers such that q+a1**2+a2**2=2**e.
	- f: integer determining the accessible 2-torsion in E1 (E1[2**f]).

	Output:
	- M: symplectic base change matrix of (*,B_Kp2) in B1.
	"""
	N=2**f
	ZN=Integers(N)

	C=matrix(ZN,[[a1,0,a2,0],
		[-a2,0,a1,0],
		[-1,0,0,0],
		[0,0,-1,0]])

	D=matrix(ZN,[[0,a1,0,a2],
		[0,-a2,0,a1],
		[0,-q,0,0],
		[0,0,0,-q]])



	M=complete_symplectic_matrix_dim4(C,D,N)

	assert is_symplectic_matrix_dim4(M)

	return M

def matrix_F_dual(a1,a2,q,f):
	r""" Computes the matrix of \tilde{F}(B1) in B1, where:
	B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1)],
	[(P2,0,0,0),(0,P2,0,0),(0,0,lamb*R2,0),(0,0,0,lamb*R2)]]
	as defined in complete_kernel_matrix_F2_dual.

	Input:
	- a1, a2, q: Integers such that q+a1**2+a2**2=2**e.
	- f: integer determining the accessible 2-torsion in E1 (E1[2**f]).

	Output:
	- M: symplectic base change matrix of \tilde{F}(B1) in B1.
	"""
	N=2**f
	ZN=Integers(N)

	M=matrix(ZN,[[a1,-a2,-q,0,0,0,0,0],
		[a2,a1,0,-q,0,0,0,0],
		[1,0,a1,a2,0,0,0,0],
		[0,1,-a2,a1,0,0,0,0],
		[0,0,0,0,a1,-a2,-1,0],
		[0,0,0,0,a2,a1,0,-1],
		[0,0,0,0,q,0,a1,a2],
		[0,0,0,0,0,q,-a2,a1]])

	return M

def matrix_F(a1,a2,q,f):
	r""" Computes the matrix of F(B1) in B1, where:
	B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1)],
	[(P2,0,0,0),(0,P2,0,0),(0,0,lamb*R2,0),(0,0,0,lamb*R2)]]
	as defined in complete_kernel_matrix_F1.

	Input:
	- a1, a2, q: Integers such that q+a1**2+a2**2=2**e.
	- f: integer determining the accessible 2-torsion in E1 (E1[2**f]).

	Output:
	- M: symplectic base change matrix of \tilde{F}(B1) in B1.
	"""
	N=2**f
	ZN=Integers(N)

	M=matrix(ZN,[[a1,a2,q,0,0,0,0,0],
		[-a2,a1,0,q,0,0,0,0],
		[-1,0,a1,-a2,0,0,0,0],
		[0,-1,a2,a1,0,0,0,0],
		[0,0,0,0,a1,a2,1,0],
		[0,0,0,0,-a2,a1,0,1],
		[0,0,0,0,-q,0,a1,-a2],
		[0,0,0,0,0,-q,a2,a1]])

	return M

def starting_two_symplectic_matrices(a1,a2,q,f):
	r"""
	Computes the matrices of two symplectic basis of E1^2*E2^2[2**f] given 
	by (*,B_Kp1) and (*,B_Kp2) in the basis
	B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1)],
	[(P2,0,0,0),(0,P2,0,0),(0,0,lamb*R2,0),(0,0,0,lamb*R2)]]
	as defined in complete_kernel_matrix_F1.

	Input:
	- a1, a2, q: Integers such that q+a1**2+a2**2=2**e.
	- f: integer determining the accessible 2-torsion in E1 (E1[2**f]).

	Output:
	- M1, M2: the symplectic base change matrices of (*,B_Kp1) and (*,B_Kp2) in B1.
	"""
	M1_0=complete_kernel_matrix_F1(a1,a2,q,f)
	MatF=matrix_F(a1,a2,q,f)

	# Matrix of an isotropic subgroup of E1^2*E2^2[2**f] lying above ker(\tilde{F2}).
	Block_right2=MatF*M1_0[:,[0,1,2,3]]

	N=ZZ(2**f)

	C=Block_right2[[0,1,2,3],:]
	D=Block_right2[[4,5,6,7],:]

	assert C.transpose()*D==D.transpose()*C

	# Matrix of the resulting symplectic basis (*,B_Kp2)
	M2=complete_symplectic_matrix_dim4(C,D,N)

	MatF_dual=matrix_F_dual(a1,a2,q,f)

	Block_right1=MatF_dual*M2[:,[0,1,2,3]]

	C=Block_right1[[0,1,2,3],:]
	D=Block_right1[[4,5,6,7],:]

	A=M1_0[[0,1,2,3],[0,1,2,3]]
	B=M1_0[[4,5,6,7],[0,1,2,3]]

	assert C.transpose()*D==D.transpose()*C
	assert B.transpose()*A==A.transpose()*B

	# Matrix of the resulting symplectic basis (*,B_Kp1)
	M1=block_matrix(1,2,[M1_0[:,[0,1,2,3]],-Block_right1])

	assert is_symplectic_matrix_dim4(M1)

	A,B,C,D=bloc_decomposition(M1)
	a2_div=a2
	m=0
	while a2_div%2==0:
		m+=1
		a2_div=a2_div//2
	for j in range(4):
		assert (-D[0,j]*a1-C[0,j]*a2-D[2,j])%2**m==0
		assert (C[0,j]*a1-D[0,j]*a2+C[2,j]*q)%2**m==0
		assert (-D[1,j]*a1-C[1,j]*a2-D[3,j])%2**m==0
		assert (C[1,j]*a1-D[1,j]*a2+C[3,j]*q)%2**m==0

	return M1, M2

def gluing_base_change_matrix_dim2_F1(a1,a2,q):
	r"""Computes the symplectic base change matrix of a symplectic basis (*,B_K4) of E1*E2[4]
	given by the kernel of the dimension 2 gluing isogeny:
	B_K4=2**(f-2)[([a1]P1-[a2]P2,R1),([a1]P2+[a2]P1,R2)]
	in the basis $2**(f-2)*B1$ given by:
	B1:=[[(P1,0),(0,R1)],[(P2,0),(0,[1/q]*R2)]]
	where:
	- lamb is the inverse of q modulo 2**f.
	- (P1,P2) is the canonical basis of E1[2**f].
	- (R1,R2) is the image of (P1,P2) by sigma.

	Input:
	- a1, q: integers.

	Output:
	- M: symplectic base change matrix of (*,B_K4) in 2**(f-2)*B1.
	"""
	return gluing_base_change_matrix_dim2(a1,a2,q)

def gluing_base_change_matrix_dim2_dim4_F1(a1,a2,q,m,M1):
	r"""Computes the symplectic base change matrix of the symplectic basis Bp of Am*Am[4] induced
	by the image of the symplectic basis (x_1, ..., x_4, y_1, ..., y_4) of E1^2*E2^2[2**(e1+2)]
	adapted to ker(F1)=[4]<y_1, ..., y_4> in the basis associated to the product theta-structure 
	of level 2 of Am*Am:
	
	B:=[(S1,0),(S2,0),(0,S1),(0,S2),(T1,0),(T2,0),(0,T1),(0,T2)]
	
	where:
	- (P1,Q1) is the canonical basis of E1[2**f].
	- (R2,S2) is the image of (P1,P2) by sigma.
	- Phi is the 2**m-isogeny E1*E2-->Am (m first steps of the chain in dimension 2).
	- S1=[2**e]Phi([lamb]Q1,[a]sigma(P1)+[b]sigma(Q1)).
	- S2=[2**e]Phi([mu]P1,[c]sigma(P1)+[d]sigma(Q1)).
	- T1=[2**(e-m)]Phi([a1]P1-[a2]Q1,sigma(P1)).
	- T2=[2**(e-m)]Phi([a1]Q1+[a2]P1,sigma(Q1)).
	- (S1,S2,T1,T2) is induced by the image by Phi of a symplectic basis of E1*E2[2**(m+2)] lying
	above the symplectic basis of E1*E2[4] outputted by gluing_base_change_matrix_dim2.

	INPUT:
	- a1, a2, q: integers.
	- m: integer (number of steps in dimension 2 and 2-adic valuation of a2).
	- M1: matrix of (x_1, ..., x_4, y_1, ..., y_4) in the symplectic basis of E1^2*E2^2[2**(e1+2)]
	given by:

	B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R2,0),(0,0,0,R2)],
	[(Q1,0,0,0),(0,Q1,0,0),(0,0,[1/q]*S2,0),(0,0,0,[1/q]*S2)]]

	OUTPUT:
	- M: symplectic base change matrix of Bp in B.
	"""

	inv_a1=inverse_mod(a1,2**(m+2))
	inv_q=inverse_mod(q,2**(m+2))
	lamb=ZZ(2**(m+1))
	mu=ZZ((1-2**(m+1)*q)*inv_a1)
	a=ZZ(2**(m+1)*a2*inv_q)
	bq=ZZ((-1-2**(m+1)*a1))
	c=ZZ(2**(m+1))
	dq=-ZZ(mu*a2)

	Z4=Integers(4)

	A,B,C,D=bloc_decomposition(M1)

	Ap=matrix(Z4,[[ZZ(-B[0,j]*a1-A[0,j]*a2-B[2,j]) for j in range(4)],
		[ZZ(A[0,j]*a1-B[0,j]*a2+A[2,j]*q) for j in range(4)],
		[ZZ(-B[1,j]*a1-A[1,j]*a2-B[3,j]) for j in range(4)],
		[ZZ(A[1,j]*a1-B[1,j]*a2+A[3,j]*q) for j in range(4)]])

	Bp=matrix(Z4,[[ZZ(2**m*(B[2,j]*a-A[0,j]*lamb-A[2,j]*bq)) for j in range(4)],
		[ZZ(2**m*(B[2,j]*c+B[0,j]*mu-A[2,j]*dq)) for j in range(4)],
		[ZZ(2**m*(B[3,j]*a-A[1,j]*lamb-A[3,j]*bq)) for j in range(4)],
		[ZZ(2**m*(B[3,j]*c+B[1,j]*mu-A[3,j]*dq)) for j in range(4)]])

	Cp=matrix(Z4,[[ZZ(ZZ(-D[0,j]*a1-C[0,j]*a2-D[2,j])//(2**m)) for j in range(4)],
		[ZZ(ZZ(C[0,j]*a1-D[0,j]*a2+C[2,j]*q)//2**m) for j in range(4)],
		[ZZ(ZZ(-D[1,j]*a1-C[1,j]*a2-D[3,j])//2**m) for j in range(4)],
		[ZZ(ZZ(C[1,j]*a1-D[1,j]*a2+C[3,j]*q)//2**m) for j in range(4)]])

	Dp=matrix(Z4,[[ZZ(D[2,j]*a-C[0,j]*lamb-C[2,j]*bq) for j in range(4)],
		[ZZ(D[0,j]*mu+D[2,j]*c-C[2,j]*dq) for j in range(4)],
		[ZZ(D[3,j]*a-C[1,j]*lamb-C[3,j]*bq) for j in range(4)],
		[ZZ(D[1,j]*mu+D[3,j]*c-C[3,j]*dq) for j in range(4)]])

	M=block_matrix(2,2,[[Ap,Cp],[Bp,Dp]])

	assert is_symplectic_matrix_dim4(M)

	return M

def gluing_base_change_matrix_dim2_F2(a1,a2,q):
	r"""Computes the symplectic base change matrix of a symplectic basis (*,B_K4) of E1*E2[4]
	given by the kernel of the dimension 2 gluing isogeny:
	B_K4=2**(f-2)[([a1]P1-[a2]P2,R1),([a1]P2+[a2]P1,R2)]
	in the basis $2**(f-2)*B1$ given by:
	B1:=[[(P1,0),(0,R1)],[(P2,0),(0,[1/q]*R2)]]
	where:
	- lamb is the inverse of q modulo 2**f.
	- (P1,P2) is the canonical basis of E1[2**f].
	- (R1,R2) is the image of (P1,P2) by sigma.

	Input:
	- a1, q: integers.

	Output:
	- M: symplectic base change matrix of (*,B_K4) in 2**(f-2)*B1.
	"""

	Z4=Integers(4)

	mu=inverse_mod(a1,4)

	A=matrix(Z4,[[0,mu],
		[0,0]])
	B=matrix(Z4,[[0,0],
		[1,-ZZ(mu*a2)]])

	C=matrix(Z4,[[ZZ(a1),-ZZ(a2)],
		[-1,0]])
	D=matrix(Z4,[[ZZ(a2),ZZ(a1)],
		[0,-ZZ(q)]])

	M=block_matrix([[A,C],[B,D]])

	assert is_symplectic_matrix_dim2(M)

	return M

def gluing_base_change_matrix_dim2_dim4_F2(a1,a2,q,m,M2):
	r"""Computes the symplectic base change matrix of the symplectic basis Bp of Am*Am[4] induced
	by the image of the symplectic basis (x_1, ..., x_4, y_1, ..., y_4) of E1^2*E2^2[2**(e1+2)]
	adapted to ker(F1)=[4]<y_1, ..., y_4> in the basis associated to the product theta-structure 
	of level 2 of Am*Am:
	
	B:=[(S1,0),(S2,0),(0,S1),(0,S2),(T1,0),(T2,0),(0,T1),(0,T2)]
	
	where:
	- (P1,Q1) is the canonical basis of E1[2**f].
	- (R2,S2) is the image of (P1,P2) by sigma.
	- Phi is the 2**m-isogeny E1*E2-->Am (m first steps of the chain in dimension 2).
	- S1=[2**e]Phi([lamb]Q1,[a]sigma(P1)+[b]sigma(Q1)).
	- S2=[2**e]Phi([mu]P1,[c]sigma(P1)+[d]sigma(Q1)).
	- T1=[2**(e-m)]Phi([a1]P1-[a2]Q1,sigma(P1)).
	- T2=[2**(e-m)]Phi([a1]Q1+[a2]P1,sigma(Q1)).
	- (S1,S2,T1,T2) is induced by the image by Phi of a symplectic basis of E1*E2[2**(m+2)] lying
	above the symplectic basis of E1*E2[4] outputted by gluing_base_change_matrix_dim2.

	INPUT:
	- a1, a2, q: integers.
	- m: integer (number of steps in dimension 2 and 2-adic valuation of a2).
	- M2: matrix of (x_1, ..., x_4, y_1, ..., y_4) in the symplectic basis of E1^2*E2^2[2**(e1+2)]
	given by:

	B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R2,0),(0,0,0,R2)],
	[(Q1,0,0,0),(0,Q1,0,0),(0,0,[1/q]*S2,0),(0,0,0,[1/q]*S2)]]

	OUTPUT:
	- M: symplectic base change matrix of Bp in B.
	"""

	inv_a1=inverse_mod(a1,2**(m+2))
	inv_q=inverse_mod(q,2**(m+2))
	lamb=ZZ(2**(m+1))
	mu=ZZ((1+2**(m+1)*q)*inv_a1)
	a=ZZ(2**(m+1)*a2*inv_q)
	bq=ZZ((1+2**(m+1)*a1))
	c=ZZ(2**(m+1))
	dq=-ZZ(mu*a2)

	Z4=Integers(4)

	A,B,C,D=bloc_decomposition(M2)

	Ap=matrix(Z4,[[ZZ(-B[0,j]*a1+A[0,j]*a2+B[2,j]) for j in range(4)],
		[ZZ(A[0,j]*a1+B[0,j]*a2-A[2,j]*q) for j in range(4)],
		[ZZ(-B[1,j]*a1+A[1,j]*a2+B[3,j]) for j in range(4)],
		[ZZ(A[1,j]*a1+B[1,j]*a2-A[3,j]*q) for j in range(4)]])

	Bp=matrix(Z4,[[ZZ(2**m*(B[2,j]*a-A[0,j]*lamb-A[2,j]*bq)) for j in range(4)],
		[ZZ(2**m*(B[2,j]*c+B[0,j]*mu-A[2,j]*dq)) for j in range(4)],
		[ZZ(2**m*(B[3,j]*a-A[1,j]*lamb-A[3,j]*bq)) for j in range(4)],
		[ZZ(2**m*(B[3,j]*c+B[1,j]*mu-A[3,j]*dq)) for j in range(4)]])

	Cp=matrix(Z4,[[ZZ(ZZ(-D[0,j]*a1+C[0,j]*a2+D[2,j])//(2**m)) for j in range(4)],
		[ZZ(ZZ(C[0,j]*a1+D[0,j]*a2-C[2,j]*q)//2**m) for j in range(4)],
		[ZZ(ZZ(-D[1,j]*a1+C[1,j]*a2+D[3,j])//2**m) for j in range(4)],
		[ZZ(ZZ(C[1,j]*a1+D[1,j]*a2-C[3,j]*q)//2**m) for j in range(4)]])

	Dp=matrix(Z4,[[ZZ(D[2,j]*a-C[0,j]*lamb-C[2,j]*bq) for j in range(4)],
		[ZZ(D[0,j]*mu+D[2,j]*c-C[2,j]*dq) for j in range(4)],
		[ZZ(D[3,j]*a-C[1,j]*lamb-C[3,j]*bq) for j in range(4)],
		[ZZ(D[1,j]*mu+D[3,j]*c-C[3,j]*dq) for j in range(4)]])

	M=block_matrix(2,2,[[Ap,Cp],[Bp,Dp]])

	assert is_symplectic_matrix_dim4(M)

	return M

def point_matrix_product(M,L_P,J=None,modulus=None):
	r"""
	Input:
	- M: matrix with (modular) integer values.
	- L_P: list of elliptic curve points [P1,P2,R1,R2] such that the rows of M correspond to the vectors
	(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1),(P2,0,0,0),(0,P2,0,0),(0,0,R2,0),(0,0,0,R2).
	- J: list of column indices (default, all the columns).
	- modulus: order of points in L_P (default, None).

	Output:
	- L_ret: list of points corresponding to the columns of M with indices in J.
	"""
	if modulus==None:
		M1=M
	else:
		Zmod=Integers(modulus)
		M1=matrix(Zmod,M)

	if J==None:
		J=range(M1.ncols())

	L_ret=[]
	for j in J:
		L_ret.append(TuplePoint(M1[0,j]*L_P[0]+M1[4,j]*L_P[1],M1[1,j]*L_P[0]+M1[5,j]*L_P[1],
			M1[2,j]*L_P[2]+M1[6,j]*L_P[3],M1[3,j]*L_P[2]+M1[7,j]*L_P[3]))

	return L_ret


def kernel_basis(M,ei,mP1,mP2,mR1,mlambR2):
	r"""
	Input:
	- M: matrix of a symplectic basis in the basis
	B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1)],
	[(P2,0,0,0),(0,P2,0,0),(0,0,lamb*R2,0),(0,0,0,lamb*R2)]]
	as defined in complete_kernel_matrix_F1.
	- ei: length of F1 or F2.
	- mP1,mP2: canonical basis (P1,P2) of E1[2**f] multiplied by m:=2**(f-ei-2).
	- mR1,mlambR2: (mR1,mlambR2)=(m*sigma(P1),m*sigma(P2)), where lamb is the
	inverse of q=deg(sigma) modulo 2**f.

	Output:
	- Basis of the second symplectic subgroup basis of E1^2*E2^2[2**(ei+2)] induced by M. 
	"""
	modulus=2**(ei+2)

	return point_matrix_product(M,[mP1,mP2,mR1,mlambR2],[4,5,6,7],modulus)

def base_change_canonical_dim4(P1,P2,R1,R2,q,f,e1,e2):
	lamb=inverse_mod(q,2**f)
	
	lambR2=lamb*R2

	P1_doubles=[P1]
	P2_doubles=[P2]
	R1_doubles=[R1]
	lambR2_doubles=[lambR2]

	for i in range(f-2):
		P1_doubles.append(2*P1_doubles[-1])
		P2_doubles.append(2*P2_doubles[-1])
		R1_doubles.append(2*R1_doubles[-1])
		lambR2_doubles.append(2*lambR2_doubles[-1])

	# Constructing canonical basis of E1[4] and E2[4].
	_,_,T1,T2,MT=make_canonical(P1_doubles[-1],P2_doubles[-1],4,preserve_pairing=True)
	_,_,U1,U2,MU=make_canonical(R1_doubles[-1],lambR2_doubles[-1],4,preserve_pairing=True)

	# Base change matrix of the symplectic basis 2**(f-2)*B1 of E1^2*E2^2[4] in the basis:
	# B0:=[[(T1,0,0,0),(0,T1,0,0),(0,0,U1,0),(0,0,0,U1)],
	#[(T2,0,0,0),(0,T2,0,0),(0,0,U2,0),(0,0,0,U2)]]
	# where B1:=[[(P1,0,0,0),(0,P1,0,0),(0,0,R1,0),(0,0,0,R1)],
	#[(P2,0,0,0),(0,P2,0,0),(0,0,lamb*R2,0),(0,0,0,lamb*R2)]]
	Z4=Integers(4)
	M0=matrix(Z4,[[MT[0,0],0,0,0,MT[1,0],0,0,0],
		[0,MT[0,0],0,0,0,MT[1,0],0,0],
		[0,0,MU[0,0],0,0,0,MU[1,0],0],
		[0,0,0,MU[0,0],0,0,0,MU[1,0]],
		[MT[0,1],0,0,0,MT[1,1],0,0,0],
		[0,MT[0,1],0,0,0,MT[1,1],0,0],
		[0,0,MU[0,1],0,0,0,MU[1,1],0],
		[0,0,0,MU[0,1],0,0,0,MU[1,1]]])

	return P1_doubles,P2_doubles,R1_doubles,lambR2_doubles,T1,T2,U1,U2,MT,MU,M0
