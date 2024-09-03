
CC = gcc

ifeq (${PRIME},pHD256)
	# SQISignHD NIST level 1
	V_DEFS_PRIME = -DNWORDS_FIELD=4 -DNWORDS_ORDER=4 -DBITS=256 -DLOG2P=8 -DFP_ENCODED_BYTES=32 -DFP2_ENCODED_BYTES=64
else ifeq (${PRIME},pHD384)
	# SQISignHD NIST level 3
	V_DEFS_PRIME = -DNWORDS_FIELD=6 -DNWORDS_ORDER=6 -DBITS=384 -DLOG2P=9 -DFP_ENCODED_BYTES=48 -DFP2_ENCODED_BYTES=96
else ifeq (${PRIME},pHD512)
	# SQISignHD NIST level 5
	V_DEFS_PRIME = -DNWORDS_FIELD=8 -DNWORDS_ORDER=8 -DBITS=512 -DLOG2P=9 -DFP_ENCODED_BYTES=64 -DFP2_ENCODED_BYTES=128
endif

V_DEFS_FILE_GF_CST = -DFILENAME=\"../src/precomp/gf_constants_${PRIME}.c\"

CFLAGS_GF_CST = -I../src/gf/include -I../src/gf/test -I../src/common

TARGETS_GF_CST = gf_constants_precomp.c ../src/gf/primes/${PRIME}.c ../src/gf/fp.c ../src/gf/fp2.c ../src/gf/mp.c ../src/gf/test/test_utils.c ../src/common/sha3.c

all: gf_constants_precomp

gf_constants_precomp: ${TARGETS_GF_CST}
	${CC} ${CFLAGS_GF_CST} ${V_DEFS_PRIME} ${V_DEFS_FILE_GF_CST} -o gf_constants_precomp.exe ${TARGETS_GF_CST}

clean: clean_gf_constants

clean_gf_constants:
	rm gf_constants_precomp.exe
