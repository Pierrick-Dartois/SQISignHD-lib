#importing sys module
import sys
# importing os module
import os
# # importing click
# import click
# importing math
import math
# importing sage
# from sage.all import *

# The following benchmarks are required for computing optimal strategies.
# The costs can be specified in arbitrary units and can be obtained by
# runing src/ec/ref/lvl1/test/mont.test or left as "None" to use a
# generic estimate.
# ++++++++++ 
#p2 =  None   # cost of xdbl
#q4 = None    # cost of xeval4

##SqiSign lvl 1 costs
p2 = 1766
q4 = 2452

##SqiSign lvl 3 costs
#p2 = 4071
#q4 = 5224

##SqiSign lvl 5 costs
#p2 = 7755
#q4 = 9847
# ++++++++++

def fp2str(x : int, p : int):
    nwords = (p.bit_length()-1)//64 + 1
    X = hex(x%p)[2:]
    while len(X) < 16*nwords:
        X='0'+X
    output = '{ '
    for i in range(nwords):
        output += '0x'+X[16*(nwords-1-i):16*(nwords-i)]+', '
    output = output[:-2]+' }'
    return output

def list2str(x):
    out = '{ '
    for a in x:
        out += str(a)+', '
    out = out[:-2] + ' }'
    return out

# Optimal strategy
def strategy(n, p, q):
    S = { 1: [] }
    C = { 1: 0 }
    for i in range(2, n+1):
        b, cost = min(((b, C[i-b] + C[b] + b*p + (i-b)*q) for b in range(1,i)),
                      key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost
    return S[n]

# Function to read sqisign_parameters.txt and initial config
def config():

    # Read file parameters
    f = open('sqisign_parameters.txt')
    f.readline()
    p_input_line = f.readline()
    B_input_line = f.readline()
    f.close()
    try:
        p = int(p_input_line.split('0x')[-1],16)
        B = int(B_input_line.split(' = ')[-1])
        assert(B_input_line[:4+len(str(B))] == 'B = '+str(B) and p_input_line[:4+len(hex(p))] == 'p = '+hex(p))
    except:
        print("Error reading sqisign_parameters.txt. Ensure file contests have the following form:\nlvl = x\np = 0x[hexstring]\nB = [decimalstring]")
        exit(-1)

    # Factorization
    factorization = factor(p+1)
    if not factorization[0][0]==2 and factorization[1][0]==3:
        print("Error: prime does not have the form p = f*2^a*3^b-1")
        exit(-1)
    POWER_OF_2 = factorization[0][1]
    POWER_OF_3 = factorization[1][1]
    Pfactors = [factor[0] for factor in factorization[1:] if factor[0] <= B]
    factorization = factor(p-1)
    Mfactors = [factor[0] for factor in factorization[1:] if factor[0] <= B]
    return p,POWER_OF_2,POWER_OF_3,Pfactors,Mfactors

if __name__ == '__main__':
    p,POWER_OF_2,POWER_OF_3,Pfactors,Mfactors = config()
    PMfactors = Pfactors + Mfactors

    if not p2 and q4:
        p2 = 1.0
        q4 = 1.3

    sizeI = [ceil(sqrt((l-1)/4)) for l in PMfactors]
    sizeJ = [floor((l-1)/4/ceil(sqrt((l-1)/4))) for l in PMfactors]
    sizeK = [int((l-1)/2 - 2*ceil(sqrt((l-1)/4))*floor((l-1)/4/ceil(sqrt((l-1)/4)))) for l in PMfactors]

    f = open('ec_params.h', 'w')
    f.write('#ifndef EC_PARAMS_H\n')
    f.write('#define EC_PARAMS_H\n')
    f.write('\n')
    f.write('#include <fp_constants.h>\n')
    f.write('\n')
    f.write(f'#define POWER_OF_2 {POWER_OF_2}\n')
    f.write(f'#define POWER_OF_3 {POWER_OF_3}\n')
    f.write('\n')
    f.write(f'static digit_t TWOpF[NWORDS_ORDER] = {fp2str(2**POWER_OF_2, p)}; // Fp representation for the power of 2\n')
    f.write(f'static digit_t TWOpFm1[NWORDS_ORDER] = {fp2str(2**(POWER_OF_2-1), p)}; // Fp representation for half the power of 2\n')
    f.write(f'static digit_t THREEpE[NWORDS_ORDER] = {fp2str(3**(POWER_OF_3//2), p)}; // Approximate squareroot of the power of 3\n')
    f.write(f'static digit_t THREEpF[NWORDS_ORDER] = {fp2str(3**(POWER_OF_3), p)}; // Fp representation for the power of 3\n')
    f.write(f'static digit_t THREEpFdiv2[NWORDS_ORDER] = {fp2str(3**(POWER_OF_3)//2, p)}; // Floor of half the power of 3\n')
    f.write('\n')
    f.write('#define scaled 1 // unscaled (0) or scaled (1) remainder tree approach for squareroot velu\n')
    f.write('#define gap 83 // Degree above which we use squareroot velu reather than traditional\n')
    f.write('\n')
    f.write(f'#define P_LEN {len(Pfactors)} // Number of odd primes in p+1\n')
    f.write(f'#define M_LEN {len(Mfactors)} // Number of odd primes in p-1\n')
    f.write('\n')
    f.write('// Bitlength of the odd prime factors\n')
    f.write(f'static digit_t p_plus_minus_bitlength[P_LEN + M_LEN] = \n\t{list2str([int(l).bit_length() for l in Pfactors] + [int(l).bit_length() for l in Mfactors])};\n')
    f.write('\n')
    f.write('// p+1 divided by the power of 2\n')
    f.write(f'static digit_t p_cofactor_for_2f[{(((p+1)//2**POWER_OF_2).bit_length()-1)//64+1}] = {fp2str((p+1)//2**POWER_OF_2, (p+1)//2**POWER_OF_2+1)};\n')
    f.write(f'#define P_COFACTOR_FOR_2F_BITLENGTH {((p+1)//2**POWER_OF_2).bit_length()}\n')
    f.write('\n')
    f.write('// p+1 divided by the power of 3\n')
    f.write(f'static digit_t p_cofactor_for_3g[{(((p+1)//3**POWER_OF_3).bit_length()-1)//64+1}] = {fp2str((p+1)//3**POWER_OF_3, (p+1)//3**POWER_OF_3+1)};\n')
    f.write(f'#define P_COFACTOR_FOR_3G_BITLENGTH {((p+1)//3**POWER_OF_3).bit_length()}\n')
    f.write('\n')
    f.write('// p+1 divided by the powers of 2 and 3\n')
    f.write(f'static digit_t p_cofactor_for_6fg[{(((p+1)//3**POWER_OF_3//2**POWER_OF_2).bit_length()-1)//64+1}] = {fp2str((p+1)//3**POWER_OF_3//2**POWER_OF_2, (p+1)//3**POWER_OF_3//2**POWER_OF_2+1)};\n')
    f.write(f'#define P_COFACTOR_FOR_6FG_BITLENGTH {((p+1)//3**POWER_OF_3//2**POWER_OF_2).bit_length()}\n')
    f.write('\n')
    f.write('// Strategy for 4-isogenies\n')
    f.write(f'static int STRATEGY4[{POWER_OF_2//2-1}] = {list2str(strategy(POWER_OF_2//2-1, 2*p2, q4))};\n')
    f.write('\n')
    f.write('// Optimal sizes for I,J,K in squareroot Velu\n')
    f.write(f'static int sizeI[{len(sizeI)}] =\n\t{list2str(sizeI)};\n')
    f.write(f'static int sizeJ[{len(sizeJ)}] =\n\t{list2str(sizeJ)};\n')
    f.write(f'static int sizeK[{len(sizeK)}] =\n\t{list2str(sizeK)};\n')
    f.write('\n')
    f.write(f'#define sI_max {max(sizeI)}\n')
    f.write(f'#define sJ_max {max(sizeJ)}\n')
    f.write(f'#define sK_max {max(max(sizeK), 83)}\n')
    f.write('\n')
    f.write(f'#define ceil_log_sI_max {ceil(log(max(sizeI),2))}\n')
    f.write(f'#define ceil_log_sJ_max {ceil(log(max(sizeJ),2))}\n')
    f.write('\n')
    f.write('#endif\n')

    f.close()


