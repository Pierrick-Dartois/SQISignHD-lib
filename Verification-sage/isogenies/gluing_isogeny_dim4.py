from sage.all import *

from theta_structures.Theta_dim4 import ThetaStructureDim4, ThetaPointDim4
from theta_structures.theta_helpers_dim4 import hadamard, squared, batch_inversion, multindex_to_index
from isogenies.tree import Tree
from isogenies.isogeny_dim4 import IsogenyDim4, DualIsogenyDim4

def proj_equal(P1, P2):
    if len(P1) != len(P2):
        return False
    for i in range(0, len(P1)):
        if P1[i]==0:
            if P2[i] != 0:
                return False
        else:
            break
    r=P1[i]
    s=P2[i]
    for i in range(0, len(P1)):
        if P1[i]*s != P2[i]*r:
            return False
    return True

class GluingIsogenyDim4(IsogenyDim4):
	def __init__(self,domain,L_K_8,L_K_8_ind):
		r"""
		Input:
		- domain: a ThetaStructureDim4.
		- L_K_8: list of points of 8-torsion in the kernel.
		- L_K_8_ind: list of corresponding multindices (i0,i1,i2,i3) of points in L_K_8
		(L_K_8[i]=i0*P0+i1*P1+i2*P2+i3*P3, where (P0,..,P3) is a basis of K_8 (compatible with
		the canonical basis of K_2) and L_K_8_ind[i]=(i0,i1,i2,i3)).
		"""

		if not isinstance(domain, ThetaStructureDim4):
			raise ValueError("Argument domain should be a ThetaStructureDim4 object.")
		self._domain = domain
		self._precomputation=None
		self._special_compute_codomain(L_K_8,L_K_8_ind)

		#a_i2=squared(self._domain.zero())
		#HB_i2=hadamard(squared(hadamard(self._codomain.zero())))
		#for i in range(16):
			#print(HB_i2[i]/a_i2[i])

	def _special_compute_codomain(self,L_K_8,L_K_8_ind):
		r"""
		Input:
		- L_K_8: list of points of 8-torsion in the kernel.
		- L_K_8_ind: list of corresponding multindices (i0,i1,i2,i3) of points in L_K_8
		(L_K_8[i]=i0*P0+i1*P1+i2*P2+i3*P3, where (P0,..,P3) is a basis of K_8 (compatible with
		the canonical basis of K_2) and L_K_8_ind[i]=(i0,i1,i2,i3)).
		
		Output:
		- codomain of the isogeny.
		Also initializes self._precomputation, containing the inverse of theta-constants.
		"""
		HSK_8=[hadamard(squared(P.coords())) for P in L_K_8]
		
		# Choice of reference index j_0<->chi_0 corresponding to a non-vanishing theta-constant.
		found_tree=False
		j_0=0
		while not found_tree:
			found_k0=False
			for k in range(len(L_K_8)):
				if HSK_8[k][j_0]!=0:
					k_0=k
					found_k0=True
					break
			if not found_k0:
				j_0+=1
			else:
				j0pk0=j_0^multindex_to_index(L_K_8_ind[k_0])
				# List of tuples of indices (index chi of the denominator: HS(f(P_k))_chi, 
				#index chi.chi_k of the numerator: HS(f(P_k))_chi.chi_k, index k).
				L_ratios_ind=[(j_0,j0pk0,k_0)]
				L_covered_ind=[j_0,j0pk0]

				# Tree containing the the theta-null points indices as nodes and the L_ratios_ind reference indices as edges.
				tree_ratios=Tree(j_0)
				tree_ratios.add_child(Tree(j0pk0),0)

				# Filling in the tree
				tree_filled=False
				while not tree_filled:
					found_j=False
					for j in L_covered_ind:
						for k in range(len(L_K_8)):
							jpk=j^multindex_to_index(L_K_8_ind[k])
							if jpk not in L_covered_ind and HSK_8[k][j]!=0:
								L_covered_ind.append(jpk)
								L_ratios_ind.append((j,jpk,k))
								tree_j=tree_ratios.look_node(j)
								tree_j.add_child(Tree(jpk),len(L_ratios_ind)-1)
								found_j=True
								#break
							#if found_j:
								#break
					if not found_j or len(L_covered_ind)==16:
						tree_filled=True
				if len(L_covered_ind)!=16:
					j_0+=1
				else:
					found_tree=True

		L_denom=[HSK_8[t[2]][t[0]] for t in L_ratios_ind]
		L_denom_inv=batch_inversion(L_denom)
		L_num=[HSK_8[t[2]][t[1]] for t in L_ratios_ind]
		L_ratios=[L_num[i]*L_denom_inv[i] for i in range(15)]

		L_coords_ind=tree_ratios.edge_product(L_ratios)

		O_coords=[ZZ(0) for i in range(16)]
		for t in L_coords_ind:
			O_coords[t[1]]=t[0]

		# Precomputation
		# TODO: optimize inversions and give precomputation to the codomain _arithmetic_precomputation
		L_prec=[]
		L_prec_ind=[]
		for i in range(16):
			if O_coords[i]!=0:
				L_prec.append(O_coords[i])
				L_prec_ind.append(i)
		L_prec_inv=batch_inversion(L_prec)
		precomputation=[None for i in range(16)]
		for i in range(len(L_prec)):
			precomputation[L_prec_ind[i]]=L_prec_inv[i]

		self._precomputation=precomputation

		for k in range(len(L_K_8)):
			for j in range(16):
				jpk=j^multindex_to_index(L_K_8_ind[k])
				assert HSK_8[k][j]*O_coords[jpk]==HSK_8[k][jpk]*O_coords[j]
		
		assert proj_equal(squared(self._domain._null_point.coords()), hadamard(squared(O_coords)))

		self._codomain=ThetaStructureDim4(hadamard(O_coords),null_point_dual=O_coords)

	def special_image(self,P,L_trans,L_trans_ind):
		r"""Used when we cannot evaluate the isogeny self because the codomain has zero 
		dual theta constants.

		Input:
		- P: ThetaPointDim4 of the domain.
		- L_trans: list of translates of P+T of P by points of 4-torsion T above the kernel.
		- L_trans_ind: list of indices of the translation 4-torsion points T.
		If L_trans[i]=\sum i_j*B_K4[j] then L_trans_ind[j]=\sum 2**j*i_j.

		Output:
		- the image of P by the isogeny self.
		"""
		HS_P=hadamard(squared(P.coords()))
		HSL_trans=[hadamard(squared(Q.coords())) for Q in L_trans]
		O_coords=self._codomain.null_point_dual()
		
		# L_lambda_inv: List of inverses of lambda_i such that: 
		# HS(P+Ti)=(lambda_i*U_{chi.chi_i,0}(f(P))*U_{chi,0}(0))_chi.
		L_lambda_inv_num=[]
		L_lambda_inv_denom=[]

		for k in range(len(L_trans)):
			for j in range(16):
				jpk=j^L_trans_ind[k]
				if HSL_trans[k][j]!=0 and O_coords[jpk]!=0:
					L_lambda_inv_num.append(HS_P[jpk]*O_coords[j])
					L_lambda_inv_denom.append(HSL_trans[k][j]*O_coords[jpk])
					break
		L_lambda_inv_denom=batch_inversion(L_lambda_inv_denom)
		L_lambda_inv=[L_lambda_inv_num[i]*L_lambda_inv_denom[i] for i in range(len(L_trans))]

		for k in range(len(L_trans)):
			for j in range(16):
				jpk=j^L_trans_ind[k]
				assert HS_P[jpk]*O_coords[j]==L_lambda_inv[k]*HSL_trans[k][j]*O_coords[jpk]

		U_fP=[]
		for i in range(16):
			if self._precomputation[i]!=None:
				U_fP.append(self._precomputation[i]*HS_P[i])
			else:
				for k in range(len(L_trans)):
					ipk=i^L_trans_ind[k]
					if self._precomputation[ipk]!=None:
						U_fP.append(self._precomputation[ipk]*HSL_trans[k][ipk]*L_lambda_inv[k])
						break

		fP=hadamard(U_fP)

		return self._codomain(fP)

	def dual(self):
		return DualIsogenyDim4(self._codomain,self._domain, hadamard=False)
