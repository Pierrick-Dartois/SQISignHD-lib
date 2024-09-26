#ifndef MP_H
#define MP_H

//#include <sqisign_namespace.h>
#include <stdbool.h>
#include <tutil.h>

// Functions taken from the GF module

void mp_add(digit_t *c, const digit_t *a, const digit_t *b, const unsigned int nwords);
digit_t mp_shiftr(digit_t *x, const unsigned int shift, const unsigned int nwords);
void multiple_mp_shiftl(digit_t *x, const unsigned int shift, const unsigned int nwords);
void mp_shiftl(digit_t *x, const unsigned int shift, const unsigned int nwords);
void MUL(digit_t *out, const digit_t a, const digit_t b);
void mp_mul(digit_t *c, const digit_t *a, const digit_t *b, unsigned int nwords);

// Functions taken from the EC module

void mp_sub(digit_t *c, const digit_t *a, const digit_t *b, const unsigned int nwords);
void select_ct(digit_t *c,
               const digit_t *a,
               const digit_t *b,
               const digit_t mask,
               const int nwords);
void swap_ct(digit_t *a, digit_t *b, const digit_t option, const int nwords);
int mp_compare(const digit_t *a, const digit_t *b, unsigned int nwords);
bool mp_is_zero(const digit_t *a, unsigned int nwords);
void mp_mul2(digit_t *c, const digit_t *a, const digit_t *b);
void mp_set_zero(digit_t *a, unsigned int nwords);
void mp_set_small(digit_t *a, const digit_t b, unsigned int nwords);
void mp_set_bit(digit_t *a, int pos);
void mp_copy(digit_t *a, const digit_t *b,  unsigned int nwords);
void mp_div(digit_t *q, const digit_t *a, const digit_t *b, const unsigned int nwords);
void mp_div_with_remainder(digit_t *q, digit_t *r, const digit_t *a, const digit_t *b, const unsigned int nwords);
uint16_t mp_nbits(const uint64_t * a, const unsigned int nwords);
void mp_decode(digit_t *d, const void *src, unsigned int nwords);

/********************** Constant-time unsigned comparisons ***********************/

// The following functions return 1 (TRUE) if condition is true, 0 (FALSE) otherwise
// TODO: these should match the API of fp and use 0xFF..FF and 0 as uint32_t
static inline unsigned int
is_digit_nonzero_ct(digit_t x)
{ // Is x != 0?
    return (unsigned int)((x | (0 - x)) >> (RADIX - 1));
}

static inline unsigned int
is_digit_zero_ct(digit_t x)
{ // Is x = 0?
    return (unsigned int)(1 ^ is_digit_nonzero_ct(x));
}

static inline unsigned int
is_digit_lessthan_ct(digit_t x, digit_t y)
{ // Is x < y?
    return (unsigned int)((x ^ ((x ^ y) | ((x - y) ^ y))) >> (RADIX - 1));
}

/********************** Platform-independent macros for digit-size operations
 * **********************/

// Digit addition with carry
#define ADDC(sumOut, carryOut, addend1, addend2, carryIn)                                          \
    {                                                                                              \
        digit_t tempReg = (addend1) + (digit_t)(carryIn);                                          \
        (sumOut) = (addend2) + tempReg;                                                            \
        (carryOut) = (is_digit_lessthan_ct(tempReg, (digit_t)(carryIn)) |                          \
                      is_digit_lessthan_ct((sumOut), tempReg));                                    \
    }

// Digit subtraction with borrow
#define SUBC(differenceOut, borrowOut, minuend, subtrahend, borrowIn)                              \
    {                                                                                              \
        digit_t tempReg = (minuend) - (subtrahend);                                                \
        unsigned int borrowReg = (is_digit_lessthan_ct((minuend), (subtrahend)) |                  \
                                  ((borrowIn) & is_digit_zero_ct(tempReg)));                       \
        (differenceOut) = tempReg - (digit_t)(borrowIn);                                           \
        (borrowOut) = borrowReg;                                                                   \
    }

// Shift right with flexible datatype
#define SHIFTR(highIn, lowIn, shift, shiftOut, DigitSize)                                          \
    (shiftOut) = ((lowIn) >> (shift)) ^ ((highIn) << (DigitSize - (shift)));

// Digit shift left
#define SHIFTL(highIn, lowIn, shift, shiftOut, DigitSize)                                          \
    (shiftOut) = ((highIn) << (shift)) ^ ((lowIn) >> (RADIX - (shift)));

#endif
