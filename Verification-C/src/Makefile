
CC = gcc

## Finite fields (gf)

ifeq (${PRIME},pHD256)
	# SQISignHD NIST level 1
	V_DEFS_PRIME = -DNWORDS_FIELD=4 -DNWORDS_ORDER=4 -DPRIME_CODE=1 -DBITS=256 -DLOG2P=8 -DFP_ENCODED_BYTES=32 -DFP2_ENCODED_BYTES=64
else ifeq (${PRIME},pHD384)
	# SQISignHD NIST level 3
	V_DEFS_PRIME = -DNWORDS_FIELD=6 -DNWORDS_ORDER=6 -DPRIME_CODE=3 -DBITS=384 -DLOG2P=9 -DFP_ENCODED_BYTES=48 -DFP2_ENCODED_BYTES=96
else ifeq (${PRIME},pHD512)
	# SQISignHD NIST level 5
	V_DEFS_PRIME = -DNWORDS_FIELD=8 -DNWORDS_ORDER=8 -DPRIME_CODE=5 -DBITS=512 -DLOG2P=9 -DFP_ENCODED_BYTES=64 -DFP2_ENCODED_BYTES=128
endif

CFLAGS_GF = -Igf/include -Icommon

TARGETS_FP = gf/fp.c gf/mp.c gf/test/sha3.c gf/primes/${PRIME}.c gf/test/test_fp.c gf/test/test_utils.c
TARGETS_FP2 = gf/fp2.c gf/fp.c gf/mp.c gf/test/sha3.c gf/primes/${PRIME}.c gf/test/test_fp2.c gf/test/test_utils.c

all: test_fp.exe #test_fp2.exe

test_fp.exe: ${TARGETS_FP}
	${CC} ${CFLAGS_GF} ${V_DEFS_PRIME} -o test_fp.exe ${TARGETS_FP}

#test_fp2.exe: ${TARGETS_FP2}
	#${CC} ${CFLAGS_GF} ${V_DEFS_PRIME} -o test_fp2.exe ${TARGETS_FP2}

clean:
	rm test_fp.exe #test_fp2.exe
