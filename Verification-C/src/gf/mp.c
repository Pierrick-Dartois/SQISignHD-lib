// TODO this should be its own module
#include <mp.h>
#include <stdio.h>

// The below functions were taken from the GF module

void
MUL(digit_t *out, const digit_t a, const digit_t b)
{ // Digit multiplication, digit*digit -> 2-digit result
  // Inputs: a, b in [0, 2^w-1], where w is the computer wordsize
  // Output: 0 < out < 2^(2w)-1
    register digit_t al, ah, bl, bh, temp;
    digit_t albl, albh, ahbl, ahbh, res1, res2, res3, carry;
    digit_t mask_low = (digit_t)(-1) >> (sizeof(digit_t) * 4), mask_high = (digit_t)(-1)
                                                                           << (sizeof(digit_t) * 4);

    al = a & mask_low;               // Low part
    ah = a >> (sizeof(digit_t) * 4); // High part
    bl = b & mask_low;
    bh = b >> (sizeof(digit_t) * 4);

    albl = al * bl;
    albh = al * bh;
    ahbl = ah * bl;
    ahbh = ah * bh;
    out[0] = albl & mask_low; // out00

    res1 = albl >> (sizeof(digit_t) * 4);
    res2 = ahbl & mask_low;
    res3 = albh & mask_low;
    temp = res1 + res2 + res3;
    carry = temp >> (sizeof(digit_t) * 4);
    out[0] ^= temp << (sizeof(digit_t) * 4); // out01

    res1 = ahbl >> (sizeof(digit_t) * 4);
    res2 = albh >> (sizeof(digit_t) * 4);
    res3 = ahbh & mask_low;
    temp = res1 + res2 + res3 + carry;
    out[1] = temp & mask_low; // out10
    carry = temp & mask_high;
    out[1] ^= (ahbh & mask_high) + carry; // out11
}

void
mp_add(digit_t *c, const digit_t *a, const digit_t *b, const unsigned int nwords)
{ // Multiprecision addition
    unsigned int i, carry = 0;

    for (i = 0; i < nwords; i++) {
        ADDC(c[i], carry, a[i], b[i], carry);
    }
}

digit_t
mp_shiftr(digit_t *x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision right shift
    digit_t bit_out = x[0] & 1;

    for (unsigned int i = 0; i < nwords - 1; i++) {
        SHIFTR(x[i + 1], x[i], shift, x[i], RADIX);
    }
    x[nwords - 1] >>= shift;
    return bit_out;
}

void
mp_shiftl(digit_t *x, const unsigned int shift, const unsigned int nwords)
{ // Multiprecision left shift of at most 64

    for (int i = nwords - 1; i > 0; i--) {
        SHIFTL(x[i], x[i - 1], shift, x[i], RADIX);
    }
    x[0] <<= shift;
}

void
multiple_mp_shiftl(digit_t *x, const unsigned int shift, const unsigned int nwords)
{
    int t = shift;
    while (t > 60) {
        mp_shiftl(x, 60, nwords);
        t = t - 60;
    }
    mp_shiftl(x, t, nwords);
}

// The below functions were taken from the EC module

void
mp_sub(digit_t *c, const digit_t *a, const digit_t *b, const unsigned int nwords)
{ // Multiprecision subtraction, assuming a > b
    unsigned int i, borrow = 0;

    for (i = 0; i < nwords; i++) {
        SUBC(c[i], borrow, a[i], b[i], borrow);
    }
}

void
select_ct(digit_t *c, const digit_t *a, const digit_t *b, const digit_t mask, const int nwords)
{ // Select c <- a if mask = 0, select c <- b if mask = 1...1

    for (int i = 0; i < nwords; i++) {
        c[i] = ((a[i] ^ b[i]) & mask) ^ a[i];
    }
}

void
swap_ct(digit_t *a, digit_t *b, const digit_t option, const int nwords)
{ // Swap entries
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then a <- b and b <- a
    digit_t temp;

    for (int i = 0; i < nwords; i++) {
        temp = option & (a[i] ^ b[i]);
        a[i] = temp ^ a[i];
        b[i] = temp ^ b[i];
    }
}

int
mp_compare(const digit_t *a, const digit_t *b, unsigned int nwords)
{ // Multiprecision comparison, a=b? : (1) a>b, (0) a=b, (-1) a<b

    for (int i = nwords - 1; i >= 0; i--) {
        if (a[i] > b[i])
            return 1;
        else if (a[i] < b[i])
            return -1;
    }
    return 0;
}

bool
mp_is_zero(const digit_t *a, unsigned int nwords)
{ // Is a multiprecision element zero?
  // Returns 1 (true) if a=0, 0 (false) otherwise
    digit_t r = 0;

    for (unsigned int i = 0; i < nwords; i++)
        r |= a[i] ^ 0;

    return (bool)is_digit_zero_ct(r);
}

void
mp_mul2(digit_t *c, const digit_t *a, const digit_t *b)
{ // Multiprecision multiplication fixed to two-digit operands
    unsigned int carry = 0;
    digit_t t0[2], t1[2], t2[2];

    MUL(t0, a[0], b[0]);
    MUL(t1, a[0], b[1]);
    ADDC(t0[1], carry, t0[1], t1[0], carry);
    ADDC(t1[1], carry, 0, t1[1], carry);
    MUL(t2, a[1], b[1]);
    ADDC(t2[0], carry, t2[0], t1[1], carry);
    ADDC(t2[1], carry, 0, t2[1], carry);
    c[0] = t0[0];
    c[1] = t0[1];
    c[2] = t2[0];
    c[3] = t2[1];
}

void 
mp_mul(digit_t *c, const digit_t *a, const digit_t *b, unsigned int nwords){
    // Multiprecision multiplication truncated at nwords
    const unsigned int cnwords=nwords;
    digit_t temp[cnwords];
    uint64_t mask;
    uint8_t is_one=0;

    mp_set_zero(c,nwords);
    mp_copy(temp,a,nwords);

    for(int i=0;i<nwords;i++){
        for(int j=0;j<RADIX;j++){
            mask=1ULL<<j;
            is_one=(mask&b[i])>>j;
            if(is_one){
                mp_add(c,c,temp,nwords);
            }
            mp_shiftl(temp,1,nwords);
        }
    }
}

void
mp_set_zero(digit_t *a, unsigned int nwords){
    for(int i=0;i<nwords;i++){
        a[i]=0;
    }
}

void 
mp_set_small(digit_t *a, const digit_t b, unsigned int nwords){
    a[0]=b;
    for(int i=1;i<nwords;i++){
        a[i]=0;
    }
}

void 
mp_set_bit(digit_t *a, int pos){
    int i_word, pos_in_word;
    digit_t b;

    i_word = pos/RADIX;
    pos_in_word = pos%RADIX;

    b = 1ULL<<pos_in_word;
    a[i_word] = a[i_word] + b;
}

void 
mp_copy(digit_t *a, const digit_t *b,  unsigned int nwords){
    for(int i=0;i<nwords;i++){
        a[i]=b[i];
    }
}

void 
mp_div(digit_t *q, const digit_t *a, const digit_t *b, const unsigned int nwords){
    // Euclidean division a=bq+r
    digit_t r[nwords], b_shifted[nwords], temp[nwords];
    int pos=0;

    mp_copy(b_shifted,b,nwords);
    mp_set_zero(q,nwords);
    mp_set_zero(r,nwords);

    if(mp_compare(a,b,nwords)>=0){
        mp_sub(r,a,b,nwords);
        while(mp_compare(r,b_shifted,nwords)==1){
            pos++;
            mp_shiftl(b_shifted,1,nwords);
            mp_sub(r,a,b_shifted,nwords);
        }
        mp_set_bit(q,pos);

        while(mp_compare(r,b,nwords)>=0){
            while(mp_compare(r,b_shifted,nwords)==-1){
                pos--;
                mp_shiftr(b_shifted,1,nwords);
            }
            mp_sub(r,r,b_shifted,nwords);
            mp_set_bit(q,pos);
        }
    }
}

void 
mp_div_with_remainder(digit_t *q, digit_t *r, const digit_t *a, const digit_t *b, const unsigned int nwords){
    // Euclidean division a=bq+r
    digit_t b_shifted[nwords], temp[nwords];
    int pos=0;

    mp_copy(b_shifted,b,nwords);
    mp_set_zero(q,nwords);
    mp_set_zero(r,nwords);

    if(mp_compare(a,b,nwords)>=0){
        mp_sub(r,a,b,nwords);
        while(mp_compare(r,b_shifted,nwords)==1){
            pos++;
            mp_shiftl(b_shifted,1,nwords);
            mp_sub(r,a,b_shifted,nwords);
        }
        mp_set_bit(q,pos);

        while(mp_compare(r,b,nwords)>=0){
            while(mp_compare(r,b_shifted,nwords)==-1){
                pos--;
                mp_shiftr(b_shifted,1,nwords);
            }
            mp_sub(r,r,b_shifted,nwords);
            mp_set_bit(q,pos);
        }
    }
}

uint16_t mp_nbits(const uint64_t * a, const unsigned int nwords){
    uint16_t nbits=nwords*RADIX;
    uint8_t is_one=0;
    uint64_t mask;
    for(int i=0;i<nwords;i++){
        for(int j=0;j<RADIX;j++){
            // 1ULL otherwise, he cannot shift more than 32 bits... C is so wierd...
            mask=1ULL<<(RADIX-1-j);
            is_one=(mask&a[nwords-1-i])>>(RADIX-1-j);
            if(is_one){
                return nbits;
            }
            else{
                nbits--;
            }
        }
    }
    return nbits;
}

// Little-endian decoding of a 64-bit integer.
static inline uint64_t
mp_dec64le(const void *src)
{
    const uint8_t *buf = src;
    return (uint64_t)buf[0] | ((uint64_t)buf[1] << 8) | ((uint64_t)buf[2] << 16) |
           ((uint64_t)buf[3] << 24) | ((uint64_t)buf[4] << 32) | ((uint64_t)buf[5] << 40) |
           ((uint64_t)buf[6] << 48) | ((uint64_t)buf[7] << 56);
}

void 
mp_decode(digit_t *d, const void *src, unsigned int nwords){
    const uint8_t *buf = src;

    for (int i = 0; i < nwords; i++) {
        d[i] = mp_dec64le(buf + i * 8);
    }
}
