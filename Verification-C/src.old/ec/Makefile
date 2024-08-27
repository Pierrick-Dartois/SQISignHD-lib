COMMON_LIB_PATH = ../common/include
GF_LIB_PATH1 = ../gf/include
GF_LIB_PATH2 = ../gf/test

CC = gcc
CFLAGS = -I${COMMON_LIB_PATH} -I${GF_LIB_PATH1} -I${GF_LIB_PATH2}
ifeq (${PRIME},hd256)
	# SQISignHD NIST level 1
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=1
	V_DATA = -n_words=4 -f=13 -l_A=2 -f_A=126 -l_B=3 -f_B=78
else ifeq (${PRIME},p_2_68_3_31)
	V_DEFS = -DNWORDS_FIELD=2 -DPRIME_CODE=3640
	V_DATA = -n_words=2 -f=11 -l_A=2 -f_A=68 -l_B=3 -f_B=31
else ifeq (${PRIME},p_2_66_3_36)
	V_DEFS = -DNWORDS_FIELD=2 -DPRIME_CODE=3641
	V_DATA = -n_words=2 -f=5 -l_A=2 -f_A=66 -l_B=3 -f_B=36
else ifeq (${PRIME},p_2_131_3_45)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31280
	V_DATA = -n_words=4 -f=25 -l_A=2 -f_A=131 -l_B=3 -f_B=45
else ifeq (${PRIME},p_2_132_3_53)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31281
	V_DATA = -n_words=4 -f=11 -l_A=2 -f_A=132 -l_B=3 -f_B=53
else ifeq (${PRIME},p_2_130_3_64)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31282
	V_DATA = -n_words=4 -f=11 -l_A=2 -f_A=130 -l_B=3 -f_B=64
else ifeq (${PRIME},p_2_134_3_69)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31283
	V_DATA = -n_words=4 -f=1 -l_A=2 -f_A=134 -l_B=3 -f_B=69
else ifeq (${PRIME},p_2_131_3_78)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=31284
	V_DATA = -n_words=4 -f=1 -l_A=2 -f_A=131 -l_B=3 -f_B=78
else ifeq (${PRIME},p_2_131_3_79)
	V_DEFS = -DNWORDS_FIELD=5 -DPRIME_CODE=31285
	V_DATA = -n_words=5 -f=25 -l_A=2 -f_A=131 -l_B=3 -f_B=79
else ifeq (${PRIME},p_2_194_3_61)
	V_DEFS = -DNWORDS_FIELD=5 -DPRIME_CODE=31920
	V_DATA = -n_words=5 -f=67 -l_A=2 -f_A=194 -l_B=3 -f_B=61
else ifeq (${PRIME},p_2_67_7_13)
	V_DEFS = -DNWORDS_FIELD=2 -DPRIME_CODE=7640
	V_DATA = -n_words=2 -f=109 -l_A=2 -f_A=67 -l_B=7 -f_B=13
else ifeq (${PRIME},p_2_70_7_23)
	V_DEFS = -DNWORDS_FIELD=3 -DPRIME_CODE=7641
	V_DATA = -n_words=3 -f=5 -l_A=2 -f_A=70 -l_B=7 -f_B=23
else ifeq (${PRIME},p_2_131_7_35)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=71280
	V_DATA = -n_words=4 -f=3 -l_A=2 -f_A=131 -l_B=7 -f_B=35
else ifeq (${PRIME},p_2_131_7_39)
	V_DEFS = -DNWORDS_FIELD=4 -DPRIME_CODE=71281
	V_DATA = -n_words=4 -f=261 -l_A=2 -f_A=131 -l_B=7 -f_B=39
else ifeq (${PRIME},p_2_132_7_43)
	V_DEFS = -DNWORDS_FIELD=5 -DPRIME_CODE=71282
	V_DATA = -n_words=5 -f=15 -l_A=2 -f_A=132 -l_B=7 -f_B=43
else ifeq (${PRIME},p_2_197_7_35)
	V_DEFS = -DNWORDS_FIELD=5 -DPRIME_CODE=71920
	V_DATA = -n_words=5 -f=75 -l_A=2 -f_A=197 -l_B=7 -f_B=35
endif
TARGETS = ec_void.c ../gf/fp2.c ../gf/fp.c ../gf/fiat_arithmetic/fp_${PRIME}.c test/ec_void-test.c #Change ec_void.c into ec.c when it is correct.

all: test_arith_data.h ec-test.exe

test_arith_data.h: 
	sage test/ec_test_data.py -arith ${V_DATA} test/test_arith_data.h

ec-test.exe: ${TARGETS}
	${CC} ${CFLAGS} ${V_DEFS} -o ec-test.exe ${TARGETS}

clean:
	rm test/test_arith_data.h ec-test.exe