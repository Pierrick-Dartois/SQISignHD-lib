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

def create_gf_constants_cfile(p,num=20,primename=None,path="../src/precomp"):
	
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
	NWORDS_P_COFACTOR_FOR_2F = ceil(P_COFACTOR_FOR_2F_BITLENGTH/64)
	str_p_cofactor_for_2f='const digit_t p_cofactor_for_2f[NWORDS_P_COFACTOR_FOR_2F] = {'
	str_hex_cofactor=hex(p_cofactor_for_2f)
	for i in range(NWORDS_P_COFACTOR_FOR_2F):
		str_p_cofactor_for_2f+='0x'+str_hex_cofactor[2+16*(NWORDS_P_COFACTOR_FOR_2F-1-i):2+16*(NWORDS_P_COFACTOR_FOR_2F-i)]
		if i<NWORDS_P_COFACTOR_FOR_2F-1:
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
		cfile.write('#include "gf_constants.h"\n')
		cfile.write('\n')
		cfile.write(str_nqr_table+';\n')
		cfile.write(str_z_nqr_table+';\n')
		cfile.write('\n')
		cfile.write(str_p_cofactor_for_2f+';\n')
		cfile.write('const uint16_t P_COFACTOR_FOR_2F_BITLENGTH = '+str(P_COFACTOR_FOR_2F_BITLENGTH)+';\n')
		cfile.write('const uint16_t POWER_OF_2 = '+str(POWER_OF_2)+';\n')

	return NWORDS_P_COFACTOR_FOR_2F

def create_gf_constants_header(d_nwords,d_prime_codes,num=20,path="../src/precomp"):
	with open(path+'/gf_constants.h','w') as hfile:
		hfile.write('#include <fp2.h>\n')

		hfile.write('\n')

		hfile.write('extern const fp2_t NQR_TABLE['+str(num)+'];\n')
		hfile.write('extern const fp2_t Z_NQR_TABLE['+str(num)+'];\n')

		hfile.write('\n')

		count=0
		for x in d_nwords:
			if count==0:
				hfile.write('#if PRIME_CODE == '+str(d_prime_codes[x])+'\n')
			else:
				hfile.write('#elif PRIME_CODE == '+str(d_prime_codes[x])+'\n')
			hfile.write('\t#define NWORDS_P_COFACTOR_FOR_2F '+str(d_nwords[x])+'\n')
			count+=1
		hfile.write('#endif\n')

		hfile.write('\n')

		hfile.write('extern const digit_t p_cofactor_for_2f[NWORDS_P_COFACTOR_FOR_2F];\n')
		hfile.write('extern const uint16_t P_COFACTOR_FOR_2F_BITLENGTH;\n')
		hfile.write('extern const uint16_t POWER_OF_2;\n')


if __name__=='__main__':
	d_primes = {'pHD256':13*2**126*3**78-1,'pHD384':2**191*3**118-1,'pHD512':31*2**256*3**158-1}
	d_prime_codes = {'pHD256':1,'pHD384':3,'pHD512':5}
	d_nwords = {}
	for x in d_primes:
		nwords=create_gf_constants_cfile(d_primes[x],num=20,primename=x,path="../src/precomp")
		d_nwords[x]=nwords
	create_gf_constants_header(d_nwords,d_prime_codes,num=20,path="../src/precomp")









		





