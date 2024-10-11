
CC = gcc

## Finite fields (gf)

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

CFLAGS_GF = -Igf/include -Icommon

TARGETS_FP = gf/fp.c gf/mp.c common/sha3.c gf/primes/${PRIME}.c gf/test/test_fp.c gf/test/test_utils.c
TARGETS_FP2 = gf/fp2.c gf/fp.c gf/mp.c common/sha3.c gf/primes/${PRIME}.c gf/test/test_fp2.c gf/test/test_utils.c

test_gf: test_fp.exe test_fp2.exe

test_fp.exe: ${TARGETS_FP}
	${CC} ${CFLAGS_GF} ${V_DEFS_PRIME} -o test_fp.exe ${TARGETS_FP}

test_fp2.exe: ${TARGETS_FP2}
	${CC} ${CFLAGS_GF} ${V_DEFS_PRIME} -o test_fp2.exe ${TARGETS_FP2}

clean_gf:
	rm test_fp.exe test_fp2.exe

## Elliptic curves (EC)

TARGETS_ARITH = gf/fp.c gf/fp2.c gf/mp.c common/sha3.c gf/primes/${PRIME}.c ec/ec.c ec/test/curve-arith-test.c ec/test/test_extras.c
TARGETS_BASIS = gf/fp.c gf/fp2.c gf/mp.c common/sha3.c precomp/gf_constants_${PRIME}.c gf/primes/${PRIME}.c ec/ec.c ec/basis.c ec/test/basis-gen-test.c ec/test/test_extras.c
TARGETS_PAIRING = gf/fp.c gf/fp2.c gf/mp.c gf/test/test_utils.c common/sha3.c precomp/gf_constants_${PRIME}.c gf/primes/${PRIME}.c ec/ec.c ec/basis.c ec/biextension.c ec/test/pairing-test.c common/tools.c
TARGETS_ISOG = gf/fp.c gf/fp2.c gf/mp.c precomp/gf_constants_${PRIME}.c gf/primes/${PRIME}.c ec/ec.c ec/basis.c ec/xisog.c ec/xeval.c ec/isog_chains.c ec/test/isog-test.c ec/strategies.c common/tools.c
CFLAGS_EC = ${CFLAGS_GF} -Iec/include -Iprecomp -Igf/test

test_ec: test_arith.exe test_basis.exe test_pairing.exe test_isog.exe

test_arith.exe: ${TARGETS_ARITH}
	${CC} ${CFLAGS_EC} ${V_DEFS_PRIME} -o test_arith.exe ${TARGETS_ARITH}

test_basis.exe: ${TARGETS_BASIS}
	${CC} ${CFLAGS_EC} ${V_DEFS_PRIME} -o test_basis.exe ${TARGETS_BASIS}

test_pairing.exe: ${TARGETS_PAIRING}
	${CC} ${CFLAGS_EC} ${V_DEFS_PRIME} -o test_pairing.exe ${TARGETS_PAIRING}

test_isog.exe: ${TARGETS_ISOG}
	${CC} ${CFLAGS_EC} ${V_DEFS_PRIME} -o test_isog.exe ${TARGETS_ISOG}


clean_ec:
	rm test_arith.exe test_basis.exe test_pairing.exe test_isog.exe

## Basic mp arithmetic

TARGETS_MP = gf/mp.c gf/test/test_mp.c

test_mp: test_mp.exe

test_mp.exe: ${TARGETS_MP}
	${CC} ${CFLAGS_GF} -o test_mp.exe ${TARGETS_MP}

clean_mp:
	rm test_mp.exe



