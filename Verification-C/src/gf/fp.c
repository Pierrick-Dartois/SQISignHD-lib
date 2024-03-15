#include "include/fp.h"

// const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0x252C9E49355147FF, 0x33A6A86587407437, 0x34E29E286B95D98C };
#if PRIME_CODE==1
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0x3fffffffffffffff, 0xa4382e87ff9dc589, 0x2827baebd5c8e56e };
#elif PRIME_CODE==3640
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0x18237a3d198fd8f };
#elif PRIME_CODE==3641
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0x29a8dd82e7638853 };
#elif PRIME_CODE==31280
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0xb8f8cdaf414c38d7, 0x7d1e };
#elif PRIME_CODE==31281
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0x26e4a79072d4430f, 0xb05e2e6 };
#elif PRIME_CODE==31282
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0xd49ecc24e1c87c2b, 0x772ec8e2ca8 };
#elif PRIME_CODE==31283
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0xa6e295d4d359fcbf, 0xa48daf8a94c02 };
#elif PRIME_CODE==31284
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0xa7ecc14ec3fa83c7, 0x62d7f37f9815e5fc };
#elif PRIME_CODE==31285
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0x325ca2136a649b97, 0xf54456618e6a6105, 0x1c };
#elif PRIME_CODE==31920
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0xe2639b854da79fe3, 0x1ae2e77fe01 };
#elif PRIME_CODE==7640
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0x4cd7368fd2d7 };
#elif PRIME_CODE==7641
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xc5a75ffd2cd154bf, 0x1da };
#elif PRIME_CODE==71280
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0x0fd313a36e115827, 0x72c0b3784f };
#elif PRIME_CODE==71281
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0x402d3e186a636697, 0x16dc233492cc7a9 };
#elif PRIME_CODE==71282
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0x5b7092bc11d2fd8f, 0x8a4cf6c83dd08772, 0x1 };
#elif PRIME_CODE==71920
    const uint64_t p[NWORDS_FIELD] =  { 0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x2e73abd6fec66f9f, 0x2cd3461afee2 };
#endif
//const uint64_t R2[NWORDS_FIELD] = { 0x00, 0x00, 0x00, 0x00 }; // not used?
//const uint64_t pp[NWORDS_FIELD] = { 0x00, 0x00, 0x00, 0x00 }; // not used?


void fp_set(digit_t* x, const digit_t val)
{ // Set field element x = val, where val has wordsize

    x[0] = val;
    for (unsigned int i = 1; i < NWORDS_FIELD; i++) {
        x[i] = 0;
    }
}

bool fp_is_equal(const digit_t* a, const digit_t* b)
{ // Compare two field elements in constant time
  // Returns 1 (true) if a=b, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= a[i] ^ b[i];

    return (bool)is_digit_zero_ct(r);
}

bool fp_is_zero(const digit_t* a)
{ // Is a field element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < NWORDS_FIELD; i++)
        r |= a[i] ^ 0;

    return (bool)is_digit_zero_ct(r);
}

void fp_copy(digit_t* out, const digit_t* a)
{
    memcpy(out, a, NWORDS_FIELD*RADIX/8);
}

void fp_neg(digit_t* out, const digit_t* a)
{ // Modular negation, out = -a mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
    unsigned int i, borrow = 0;

    for (i = 0; i < NWORDS_FIELD; i++) {
        SUBC(out[i], borrow, ((digit_t*)p)[i], a[i], borrow);
    }
    fp_sub(out, out, (digit_t*)p);
}

void MUL(digit_t* out, const digit_t a, const digit_t b)
{ // Digit multiplication, digit*digit -> 2-digit result 
  // Inputs: a, b in [0, 2^w-1], where w is the computer wordsize 
  // Output: 0 < out < 2^(2w)-1    
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t)*4), mask_high = (digit_t)(-1) << (sizeof(digit_t)*4);

    al = a & mask_low;                        // Low part
    ah = a >> (sizeof(digit_t)*4);            // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t)*4);

    albl = al * bl;
    albh = al * bh;
    ahbl = ah * bl;
    ahbh = ah * bh;
    out[0] = albl & mask_low;                 // out00

    res1 = albl >> (sizeof(digit_t)*4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t)*4);
    out[0] ^= temp << (sizeof(digit_t)*4);    // out01   

    res1 = ahbl >> (sizeof(digit_t)*4);
    res2 = albh >> (sizeof(digit_t)*4);
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    out[1] = temp & mask_low;                 // out10 
    carry = temp & mask_high;
    out[1] ^= (ahbh & mask_high) + carry;     // out11
}

digit_t mp_shiftr(digit_t* x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision right shift
    digit_t bit_out = x[0] & 1;

    for (unsigned int i = 0; i < nwords-1; i++) {
        SHIFTR(x[i+1], x[i], shift, x[i], RADIX);
    }
    x[nwords-1] >>= shift;
    return bit_out;
}

void mp_shiftl(digit_t* x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision left shift

    for (int i = nwords-1; i > 0; i--) {
        SHIFTL(x[i], x[i-1], shift, x[i], RADIX);
    }
    x[0] <<= shift;
}

static void fp_exp3div4(digit_t* out, const digit_t* a)
{ // Fixed exponentiation out = a^((p-3)/4) mod p
  // Input: a in [0, p-1] 
  // Output: out in [0, p-1] 
  // Requirement: p = 3(mod 4)
    fp_t p_t, acc;
    digit_t bit;

    memcpy((digit_t*)p_t, (digit_t*)p, NWORDS_FIELD*RADIX/8);
    memcpy((digit_t*)acc, (digit_t*)a, NWORDS_FIELD*RADIX/8);
    mp_shiftr(p_t, 1, NWORDS_FIELD);
    mp_shiftr(p_t, 1, NWORDS_FIELD);
    fp_set(out, 1);
    fp_tomont(out, out);

    for (int i = 0; i < NWORDS_FIELD*RADIX-2; i++) {
        bit = p_t[0] & 1;
        mp_shiftr(p_t, 1, NWORDS_FIELD);
        if (bit == 1) {
            fp_mul(out, out, acc);
        }
        fp_sqr(acc, acc);
    }
}

void fp_inv(digit_t* a)
{ // Modular inversion, out = x^-1*R mod p, where R = 2^(w*nwords), w is the computer wordsize and nwords is the number of words to represent p
  // Input: a=xR in [0, p-1] 
  // Output: out in [0, p-1]. It outputs 0 if the input does not have an inverse  
  // Requirement: Ceiling(Log(p)) < w*nwords
    fp_t t;

    fp_exp3div4(t, a);
    fp_sqr(t, t);
    fp_sqr(t, t);
    fp_mul(a, t, a);    // a^(p-2)
}

bool fp_is_square(const digit_t* a)
{ // Is field element a square?
  // Output: out = 0 (false), 1 (true)
    fp_t t, one;

    fp_exp3div4(t, a);
    fp_sqr(t, t);
    fp_mul(t, t, a);    // a^((p-1)/2)
    fp_frommont(t, t);
    fp_set(one, 1);

    return fp_is_equal(t, one);
}

void fp_sqrt(digit_t* a)
{ // Square root computation, out = a^((p+1)/4) mod p
    fp_t t;

    fp_exp3div4(t, a);
    fp_mul(a, t, a);    // a^((p+1)/4)
}