COMMON_LIB_PATH = ../common/include
COMMON_LIB_FILES = bench.h

CC = gcc
CFLAGS = -I${COMMON_LIB_PATH}
ifeq (${PRIME},hd256)
	# SQISignHD NIST level 1
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=1
else ifeq (${PRIME},p_2_68_3_31)
	V_DEFS = -DNWORDS_FIELD=2 -DPRIME_CODE=3640
else ifeq (${PRIME},p_2_66_3_36)
	V_DEFS = -DNWORDS_FIELD=2 -DPRIME_CODE=3641
else ifeq (${PRIME},p_2_131_3_45)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31280
else ifeq (${PRIME},p_2_132_3_53)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31281
else ifeq (${PRIME},p_2_130_3_64)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31282
else ifeq (${PRIME},p_2_134_3_69)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31283
else ifeq (${PRIME},p_2_131_3_78)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31284
else ifeq (${PRIME},p_2_131_3_79)
	V_DEFS = -DNWORDS_FIELD=5 -DPRIME_CODE=31285
else ifeq (${PRIME},p_2_194_3_61)
	V_DEFS = -DNWORDS_FIELD=5 -DPRIME_CODE=31920
else ifeq (${PRIME},p_2_67_7_13)
	V_DEFS = -DNWORDS_FIELD=2 -DPRIME_CODE=7640
else ifeq (${PRIME},p_2_70_7_23)
	V_DEFS = -DNWORDS_FIELD=3 -DPRIME_CODE=7641
else ifeq (${PRIME},p_2_131_7_35)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=71280
else ifeq (${PRIME},p_2_131_7_39)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=71281
else ifeq (${PRIME},p_2_132_7_43)
	V_DEFS = -DNWORDS_FIELD=5 -DPRIME_CODE=71282
else ifeq (${PRIME},p_2_197_7_35)
	V_DEFS = -DNWORDS_FIELD=5 -DPRIME_CODE=71920
endif
TARGETS_FP = fp.c fiat_arithmetic/fp_${PRIME}.c test/test_fp.c test/test_extras.c
TARGETS_FP2 = fp2.c fp.c fiat_arithmetic/fp_${PRIME}.c test/test_fp2.c test/test_extras.c



all: test_fp.exe test_fp2.exe

test_fp.exe: ${TARGETS_FP}
	${CC} ${CFLAGS} ${V_DEFS} -o test_fp.exe ${TARGETS_FP}

test_fp2.exe: ${TARGETS_FP2}
	${CC} ${CFLAGS} ${V_DEFS} -o test_fp2.exe ${TARGETS_FP2}

clean:
	rm test_fp.exe test_fp2.exe