#include "include/fp.h"

// TODO should we modify fp_sqrt to take an out
// variable?
void
fp_sqrt(fp_t *x)
{
    fp_t tmp = *x;
    (void)gf27500_sqrt(x, &tmp);
}

uint32_t
fp_is_square(const fp_t *a)
{
    // ls is (0, 1, -1) and we want fp_is_square
    // to return 0xFF..FF when ls is 1 or 0 and 0x00..00 otherwise
    int32_t ls = gf27500_legendre(a);
    return ~(uint32_t)(ls >> 1);
}

// TODO should we modify fp_inv to take an out
// variable?
void
fp_inv(fp_t *x)
{
    fp_t tmp = *x;
    (void)gf27500_invert(x, &tmp);
}
