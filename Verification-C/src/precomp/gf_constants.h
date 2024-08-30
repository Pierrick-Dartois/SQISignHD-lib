#include <fp2.h>

extern const fp2_t NQR_TABLE[20];
extern const fp2_t Z_NQR_TABLE[20];

#if PRIME_CODE == 1
	#define NWORDS_P_COFACTOR_FOR_2F 2
#elif PRIME_CODE == 3
	#define NWORDS_P_COFACTOR_FOR_2F 3
#elif PRIME_CODE == 5
	#define NWORDS_P_COFACTOR_FOR_2F 4
#endif

extern const digit_t p_cofactor_for_2f[NWORDS_P_COFACTOR_FOR_2F];
extern const uint16_t P_COFACTOR_FOR_2F_BITLENGTH;
extern const uint16_t POWER_OF_2;
