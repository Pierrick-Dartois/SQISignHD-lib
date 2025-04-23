#include "fp2.h"
#include <encoded_sizes.h>
#include <inttypes.h>

/* Arithmetic modulo X^2 + 1 */

void
fp2_encode(void *dst, const fp2_t *a)
{
    uint8_t *buf = dst;
    fp_encode(buf, &(a->re));
    fp_encode(buf + FP_ENCODED_BYTES, &(a->im));
}

uint32_t
fp2_decode(fp2_t *d, const void *src)
{
    const uint8_t *buf = src;
    uint32_t re, im;

    re = fp_decode(&(d->re), buf);
    im = fp_decode(&(d->im), buf + FP_ENCODED_BYTES);
    return re & im;
}

void
fp2_inv(fp2_t *x)
{
    fp_t t0, t1;

    fp_sqr(&t0, &(x->re));
    fp_sqr(&t1, &(x->im));
    fp_add(&t0, &t0, &t1);
    fp_inv(&t0);
    fp_mul(&(x->re), &(x->re), &t0);
    fp_mul(&(x->im), &(x->im), &t0);
    fp_neg(&(x->im), &(x->im));
}

void
fp2_batched_inv(fp2_t *x, int len)
{
    fp2_t t1[len], t2[len];
    fp2_t inverse;

    // x = x0,...,xn
    // t1 = x0, x0*x1, ... ,x0 * x1 * ... * xn
    t1[0] = x[0];
    for (int i = 1; i < len; i++) {
        fp2_mul(&t1[i], &t1[i - 1], &x[i]);
    }

    // inverse = 1/ (x0 * x1 * ... * xn)
    inverse = t1[len - 1];
    fp2_inv(&inverse);
    t2[0] = inverse;

    // t2 = 1/ (x0 * x1 * ... * xn), 1/ (x0 * x1 * ... * x(n-1)) , ... , 1/xO
    for (int i = 1; i < len; i++) {
        fp2_mul(&t2[i], &t2[i - 1], &x[len - i]);
    }

    x[0] = t2[len - 1];
    for (int i = 1; i < len; i++) {
        fp2_mul(&x[i], &t1[i - 1], &t2[len - i - 1]);
    }
}

uint32_t
fp2_is_square(const fp2_t *x)
{
    fp_t t0, t1;

    fp_sqr(&t0, &(x->re));
    fp_sqr(&t1, &(x->im));
    fp_add(&t0, &t0, &t1);

    return fp_is_square(&t0);
}

void
fp2_sqrt(fp2_t *x)
// x^p = (x0 + i*x1)^p = x0 - i*x1  (Frobenius automorphism)
// Thus: x^(p+1) = (x0 + i*x1)*(x0 - i*x1) = x0^2 + x1^2, which
// is an element of GF(p). All elements of GF(p) are squares in
// GF(p^2), but x0^2 + x1^2 is not necessarily a square in GF(p).
//
// Let conj(p) = x^p = x0 - i*x1. Note that conj() is analogous to
// the conjugate in complex numbers. In particular:
//    conj(a + b) = conj(a) + conj(b)
//    conj(a * b) = conj(a) * conj(b)
// This implies that conj(x) is a square if and only if x is a
// square, and conj(sqrt(x)) = sqrt(conj(x)). Thus, if x is a
// square, then:
//    (sqrt(x)*conj(sqrt(x)))^2 = x*conj(x) = x0^2 + x1^2
// But sqrt(x)*conj(sqrt(x)) is in GF(p); therefore, if x is a
// square, then x0^2 + x1^2 must be a square in GF(p).
{
    fp_t sqrt_delta, tmp;
    fp_t y0, y1;

    // sqrt_delta = sqrt(re^2 + im^2)
    fp_sqr(&sqrt_delta, &(x->re));
    fp_sqr(&tmp, &(x->im));
    fp_add(&sqrt_delta, &sqrt_delta, &tmp);
    fp_sqrt(&sqrt_delta);

    // y0^2 = (x0 + sqrt(delta)) / 2
    fp_add(&y0, &(x->re), &sqrt_delta);
    fp_half(&y0, &y0);

    // Now if the imaginary part of x = 0
    // then we instead want to set y0^2 = x0
    uint32_t x1_is_zero = fp_is_zero(&(x->im));
    fp_select(&y0, &y0, &(x->re), x1_is_zero);

    // If y0^2 is a square, nqr = 0, otherwise
    // nqr = 0xFF...FF
    uint32_t nqr = ~fp_is_square(&y0);

    // Now we need to check if y02 is a square!
    // If y0 is not a square:
    //     If x1 = 0 then we want to set y0 = -x0
    //     If x1 != 0, then we want to set y0 = y0 - sqrt(delta)
    fp_neg(&tmp, &y0);
    fp_select(&y0, &y0, &tmp, nqr & x1_is_zero);
    fp_sub(&tmp, &y0, &sqrt_delta);
    fp_select(&y0, &y0, &tmp, nqr & ~x1_is_zero);

    // Now we take the square root
    fp_sqrt(&y0);
    fp_add(&tmp, &y0, &y0);
    fp_inv(&tmp);
    fp_mul(&y1, &(x->im), &tmp);

    // If x1 = 0 then the sqrt worked and y1 = 0, but if
    // x0 was not a square, we must swap y0 and y1
    fp_cswap(&y0, &y1, nqr & x1_is_zero);

    //
    // Sign management
    //
    // To ensure deterministic square-root we conditionally negate the output

    // Check if y0 is zero
    uint32_t y0_is_zero = fp_is_zero(&y0);

    // Check whether y0, y1 are odd
    // Note: we encode to ensure canonical representation
    uint8_t tmp_bytes[FP_ENCODED_BYTES];
    fp_encode(tmp_bytes, &y0);
    uint32_t y0_is_odd = -((uint32_t)tmp_bytes[0] & 1);
    fp_encode(tmp_bytes, &y1);
    uint32_t y1_is_odd = -((uint32_t)tmp_bytes[0] & 1);

    // We negate the output if:
    // y0 is odd, or
    // y0 is zero and y1 is odd
    uint32_t negate_output = y0_is_odd | (y0_is_zero & y1_is_odd);
    fp_neg(&tmp, &y0);
    fp_select(&(x->re), &y0, &tmp, negate_output);
    fp_neg(&tmp, &y1);
    fp_select(&(x->im), &y1, &tmp, negate_output);
}

// exponentiation
// TODO could be improved
void
fp2_pow_vartime(fp2_t *out, const fp2_t *x, const digit_t *exp, const int size)
{
    fp2_t acc;
    digit_t bit;

    fp2_copy(&acc, x);
    fp2_set_one(out);

    // Iterate over each word of exp
    for (int j = 0; j < size; j++) {
        // Iterate over each bit of the word
        for (int i = 0; i < RADIX; i++) {
            bit = (exp[j] >> i) & 1;
            if (bit == 1) {
                fp2_mul(out, out, &acc);
            }
            fp2_sqr(&acc, &acc);
        }
    }
}

void
fp2_print(char *name, const fp2_t *a)
{
    printf("%s0x", name);

    uint8_t buf[FP_ENCODED_BYTES];
    fp_encode(&buf, &a->re); // Encoding ensures canonical rep
    for (int i = 0; i < FP_ENCODED_BYTES; i++) {
        printf("%02x", buf[FP_ENCODED_BYTES - i - 1]);
    }

    printf(" + i*0x");

    fp_encode(&buf, &a->im);
    for (int i = 0; i < FP_ENCODED_BYTES; i++) {
        printf("%02x", buf[FP_ENCODED_BYTES - i - 1]);
    }
    printf("\n");
}
