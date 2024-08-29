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

# Returns cofactor (p+1)/2^f with f maximal, the bit length of this cofactor and f
def compute_cofactor_for_2f(p):
	N=p+1
	f=0
	while N%2==0:
		f+=1
		N=N//2

	bit_length=ceil(log(N*1.)/log(2.))

	return N, bit_length, f

def create_gf_constants(p,num=20,primename=None,path="../src/precomp"):
	
	nqr_table, z_nqr_table=compute_tables(p,num)
	str_nqr_table='const fp2_t NQR_TABLE['+str(num)+'] = {'
	str_z_nqr_table='const fp2_t Z_NQR_TABLE['+str(num)+'] = {'
	NWORDS_FIELD=ceil(log(p*1.)/(64*log(2.)))
	for x in nqr_table:
		str_nqr_table+=' {{ '
		str_xre=hex(x[0])
		for i in range(NWORDS_FIELD):
			str_nqr_table+='0x'+str_xre[2+16*(NWORDS_FIELD-1-i):2+16*(NWORDS_FIELD-i)]
			if i<NWORDS_FIELD-1:
				str_nqr_table+=', '
			else:
				str_nqr_table+='}, {'
		str_xim=hex(x[1])
		for i in range(NWORDS_FIELD):
			str_nqr_table+='0x'+str_xim[2+16*(NWORDS_FIELD-1-i):2+16*(NWORDS_FIELD-i)]
			if i<NWORDS_FIELD-1:
				str_nqr_table+=', '
			else:
				str_nqr_table+='}},'
	str_nqr_table=str_nqr_table[:-1]+'}'

	for x in z_nqr_table:
		str_z_nqr_table+=' {{ '
		str_xre=hex(x[0])
		for i in range(NWORDS_FIELD):
			str_z_nqr_table+='0x'+str_xre[2+16*(NWORDS_FIELD-1-i):2+16*(NWORDS_FIELD-i)]
			if i<NWORDS_FIELD-1:
				str_z_nqr_table+=', '
			else:
				str_z_nqr_table+='}, {'
		str_xim=hex(x[1])
		for i in range(NWORDS_FIELD):
			str_z_nqr_table+='0x'+str_xim[2+16*(NWORDS_FIELD-1-i):2+16*(NWORDS_FIELD-i)]
			if i<NWORDS_FIELD-1:
				str_z_nqr_table+=', '
			else:
				str_z_nqr_table+='}},'
	str_z_nqr_table=str_z_nqr_table[:-1]+'}'

	p_cofactor_for_2f, P_COFACTOR_FOR_2F_BITLENGTH, POWER_OF_2 = compute_cofactor_for_2f(p)
	str_p_cofactor_for_2f='const digit_t p_cofactor_for_2f[NWORDS_FIELD] = {'
	str_hex_cofactor=hex(p_cofactor_for_2f)
	for i in range(NWORDS_FIELD):
		chunk=str_hex_cofactor[2+16*(NWORDS_FIELD-1-i):2+16*(NWORDS_FIELD-i)]
		if chunk=='':
			str_p_cofactor_for_2f+='0x0'
		else:
			str_p_cofactor_for_2f+='0x'+chunk
		if i<NWORDS_FIELD-1:
			str_p_cofactor_for_2f+=', '
		else:
			str_p_cofactor_for_2f+='}'



	if primename is None:
		primename=str(p)
	with open(path+'/gf_constants_'+primename+'.c','w') as cfile:
		cfile.write('#include <stddef.h>\n')
		cfile.write('#include <stdint.h>\n')
		cfile.write('#include <tutil.h>\n')
		cfile.write('#include <fp2.h>\n')
		cfile.write('\n')
		cfile.write(str_nqr_table+';\n')
		cfile.write(str_z_nqr_table+';\n')
		cfile.write('\n')
		cfile.write(str_p_cofactor_for_2f+';\n')
		cfile.write('const uint16_t P_COFACTOR_FOR_2F_BITLENGTH = '+str(P_COFACTOR_FOR_2F_BITLENGTH)+';\n')
		cfile.write('const uint16_t POWER_OF_2 = '+str(POWER_OF_2)+';\n')

if __name__=='__main__':
	d_primes={'pHD256':13*2**126*3**78-1,'pHD384':2**191*3**118-1,'pHD512':31*2**256*3**158-1}
	for x in d_primes:
		create_gf_constants(d_primes[x],num=20,primename=x,path="../src/precomp")








		





