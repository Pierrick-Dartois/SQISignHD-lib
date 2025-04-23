#ifndef FP_H
#define FP_H

// Include statements
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <tutil.h>
#include <fp_constants.h>

#include "gf5248.h"

// Type for elements of GF(p)
#define fp_t gf5248

// Constants (Assumed to be in Montgomery form)
#define ZERO gf5248_ZERO
#define ONE gf5248_ONE

// Operations in fp
#define fp_neg gf5248_neg
#define fp_add gf5248_add
#define fp_sub gf5248_sub
#define fp_mul gf5248_mul
#define fp_sqr gf5248_square
#define fp_half gf5248_half

// Constant time selection and swapping
#define fp_select gf5248_select
#define fp_cswap gf5248_cswap

// Comparisons for fp elements
#define fp_is_zero gf5248_iszero
#define fp_is_equal gf5248_equals

// Set a uint32 to an Fp value
#define fp_set_small gf5248_set_small

// Encoding and decoding of bytes
#define fp_encode gf5248_encode
#define fp_decode gf5248_decode
#define fp_decode_reduce gf5248_decode_reduce

// These functions are essentially useless because we can just
// use = for the shallow copies we need, but they're here for
// now until we do a larger refactoring
static inline void
fp_copy(fp_t *out, const fp_t *a)
{
    memcpy(out, a, sizeof(fp_t));
}

static inline void
fp_set_zero(fp_t *a)
{
    memcpy(a, &ZERO, sizeof(fp_t));
}

static inline void
fp_set_one(fp_t *a)
{
    memcpy(a, &ONE, sizeof(fp_t));
}

// Functions defined in low level code but with different API
void fp_inv(fp_t *a);
void fp_sqrt(fp_t *a);
uint32_t fp_is_square(const fp_t *a);

#endif
