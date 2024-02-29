#ifndef EC_PARAMS_H
#define EC_PARAMS_H

#include <fp_constants.h>

#define POWER_OF_2 126
#define POWER_OF_3 78

static digit_t TWOpF[NWORDS_ORDER] = { 0x0000000000000000, 0x4000000000000000, 0x0000000000000000, 0x0000000000000000 }; // Fp representation for the power of 2
static digit_t TWOpFm1[NWORDS_ORDER] = { 0x0000000000000000, 0x2000000000000000, 0x0000000000000000, 0x0000000000000000 }; // Fp representation for half the power of 2
static digit_t THREEpE[NWORDS_ORDER] = { 0x383d9170b85ff80b, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 }; // Approximate squareroot of the power of 3
static digit_t THREEpF[NWORDS_ORDER] = { 0x94fd9829d87f5079, 0x0c5afe6ff302bcbf, 0x0000000000000000, 0x0000000000000000 }; // Fp representation for the power of 3
static digit_t THREEpFdiv2[NWORDS_ORDER] = { 0xca7ecc14ec3fa83c, 0x062d7f37f9815e5f, 0x0000000000000000, 0x0000000000000000 }; // Floor of half the power of 3

#define scaled 1 // unscaled (0) or scaled (1) remainder tree approach for squareroot velu
#define gap 83 // Degree above which we use squareroot velu reather than traditional

#define P_LEN 2 // Number of odd primes in p+1
#define M_LEN 1 // Number of odd primes in p-1

// Bitlength of the odd prime factors
static digit_t p_plus_minus_bitlength[P_LEN + M_LEN] = 
	{ 2, 4, 4 };

// p+1 divided by the power of 2
static digit_t p_cofactor_for_2f[2] = { 0x90e0ba1ffe771625, 0xa09eebaf572395ba };
#define P_COFACTOR_FOR_2F_BITLENGTH 128

// p+1 divided by the power of 3
static digit_t p_cofactor_for_3g[3] = { 0x0000000000000000, 0x4000000000000000, 0x0000000000000003 };
#define P_COFACTOR_FOR_3G_BITLENGTH 130

// p+1 divided by the powers of 2 and 3
static digit_t p_cofactor_for_6fg[1] = { 0x000000000000000d };
#define P_COFACTOR_FOR_6FG_BITLENGTH 4

// Strategy for 4-isogenies
static int STRATEGY4[62] = { 27, 15, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1 };

// Optimal sizes for I,J,K in squareroot Velu
static int sizeI[3] =
	{ 1, 2, 2 };
static int sizeJ[3] =
	{ 0, 1, 1 };
static int sizeK[3] =
	{ 1, 2, 1 };

#define sI_max 2
#define sJ_max 1
#define sK_max 83

#define ceil_log_sI_max 1
#define ceil_log_sJ_max 0

#endif
