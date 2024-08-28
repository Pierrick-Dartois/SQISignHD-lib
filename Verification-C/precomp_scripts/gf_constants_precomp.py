from sage.all import *

# Returns a table of non square elements of Fp2
def compute_nqr_table(Fp2,num=20):
	nqr_table=[]
	while(len(nqr_table)<num):
		x=Fp2.random_element()
		if not x.is_square():
			nqr_table.append(x)
	return nqr_table

# Returns a table of square elements z of Fp2 such that z-1 is non-square
def compute_z_nqr_table(Fp2,num=20):
	z_nqr_table=[]
	while(len(z_nqr_table)<num):
		z=Fp2.random_element()
		if z.is_square() and not (z-1).is_square():
			z_nqr_table.append(z)
	return z_nqr_table

def compute_tables(p,num=20):
	Fp2=GF(p**2,'i',modulus=[1,0,1],proof=False)

	return compute_nqr_table(Fp2,num), compute_z_nqr_table(Fp2,num)

def create_gf_constants(p,num=20,primename=None,path="../src/precomp"):
	
	nqr_table, z_nqr_table=compute_tables(p,num)
	str_nqr_table='const fp2_t NQR_TABLE['+str(num)+'] = {'
	for x in nqr_table:
		pass

	if primename is None:
		primename=str(p)
	with open('gf_constants_'+primename+'.c','w') as cfile:
		cfile.writerow('#include <stddef.h>')
		cfile.writerow('#include <stdint.h>')
		cfile.writerow('#include <tutil.h>')
		cfile.writerow('#include <fp2.h>')

		





