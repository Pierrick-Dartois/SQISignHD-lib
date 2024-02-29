#include <stddef.h>
#include <stdint.h>
#include <torsion_constants.h>
#if 0
#elif 8*DIGIT_LEN == 16
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7e;
const uint64_t TORSION_ODD_PRIMES[3] = {0x3, 0xd, 0xb};
const uint64_t TORSION_ODD_POWERS[3] = {0x4e, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[2] = {0x3, 0xd};
const size_t TORSION_PLUS_ODD_POWERS[2] = {0x4e, 0x1};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0xb};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x1};
const size_t DEGREE_COMMITMENT_POWERS[3] = {0x0, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0xffff,0x3fff,0xc589,0xff9d,0x2e87,0xa438,0xe56e,0xd5c8,0xbaeb,0x2827}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 9, ._mp_d = (mp_limb_t[]) {0xf397,0xef1d,0xff5f,0x39a7,0x6f04,0xbe87,0x2088,0xe6d4,0x6}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[3] = {{{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x5079,0xd87f,0x9829,0x94fd,0xbcbf,0xf302,0xfe6f,0xc5a}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xd}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x1625,0xfe77,0xba1f,0x90e0,0x95ba,0x5723,0xebaf,0xa09e}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x4000}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x5079,0xd87f,0x9829,0x94fd,0xbcbf,0xf302,0xfe6f,0xc5a}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x4000,0xd41e,0x761f,0x660a,0xe53f,0xaf2f,0xfcc0,0xbf9b,0x316}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x8f}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xd}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 16, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x4000,0xd41e,0x761f,0x660a,0xe53f,0xaf2f,0xfcc0,0xbf9b,0x316}}};
#elif 8*DIGIT_LEN == 32
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7e;
const uint64_t TORSION_ODD_PRIMES[3] = {0x3, 0xd, 0xb};
const uint64_t TORSION_ODD_POWERS[3] = {0x4e, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[2] = {0x3, 0xd};
const size_t TORSION_PLUS_ODD_POWERS[2] = {0x4e, 0x1};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0xb};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x1};
const size_t DEGREE_COMMITMENT_POWERS[3] = {0x0, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0xffffffff,0xffffffff,0xffffffff,0x3fffffff,0xff9dc589,0xa4382e87,0xd5c8e56e,0x2827baeb}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 5, ._mp_d = (mp_limb_t[]) {0xef1df397,0x39a7ff5f,0xbe876f04,0xe6d42088,0x6}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[3] = {{{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xd87f5079,0x94fd9829,0xf302bcbf,0xc5afe6f}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xd}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xfe771625,0x90e0ba1f,0x572395ba,0xa09eebaf}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x40000000}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xd87f5079,0x94fd9829,0xf302bcbf,0xc5afe6f}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x40000000,0x761fd41e,0xe53f660a,0xfcc0af2f,0x316bf9b}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x8f}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xd}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 8, ._mp_d = (mp_limb_t[]) {0x0,0x0,0x0,0x40000000,0x761fd41e,0xe53f660a,0xfcc0af2f,0x316bf9b}}};
#elif 8*DIGIT_LEN == 64
const uint64_t TORSION_PLUS_EVEN_POWER = 0x7e;
const uint64_t TORSION_ODD_PRIMES[3] = {0x3, 0xd, 0xb};
const uint64_t TORSION_ODD_POWERS[3] = {0x4e, 0x1, 0x1};
const uint64_t TORSION_PLUS_ODD_PRIMES[2] = {0x3, 0xd};
const size_t TORSION_PLUS_ODD_POWERS[2] = {0x4e, 0x1};
const uint64_t TORSION_MINUS_ODD_PRIMES[1] = {0xb};
const size_t TORSION_MINUS_ODD_POWERS[1] = {0x1};
const size_t DEGREE_COMMITMENT_POWERS[3] = {0x0, 0x1, 0x1};
const ibz_t CHARACTERISTIC = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0xffffffffffffffff,0x3fffffffffffffff,0xa4382e87ff9dc589,0x2827baebd5c8e56e}}};
const ibz_t TORSION_ODD = {{._mp_alloc = 0, ._mp_size = 3, ._mp_d = (mp_limb_t[]) {0x39a7ff5fef1df397,0xe6d42088be876f04,0x6}}};
const ibz_t TORSION_ODD_PRIMEPOWERS[3] = {{{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x94fd9829d87f5079,0xc5afe6ff302bcbf}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xd}}}, {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}}};
const ibz_t TORSION_ODD_PLUS = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x90e0ba1ffe771625,0xa09eebaf572395ba}}};
const ibz_t TORSION_ODD_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}};
const ibz_t TORSION_PLUS_2POWER = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x0,0x4000000000000000}}};
const ibz_t TORSION_PLUS_3POWER = {{._mp_alloc = 0, ._mp_size = 2, ._mp_d = (mp_limb_t[]) {0x94fd9829d87f5079,0xc5afe6ff302bcbf}}};
const ibz_t TORSION_PLUS_23POWER = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x4000000000000000,0xe53f660a761fd41e,0x316bf9bfcc0af2f}}};
const ibz_t DEGREE_COMMITMENT = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0x8f}}};
const ibz_t DEGREE_COMMITMENT_PLUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xd}}};
const ibz_t DEGREE_COMMITMENT_MINUS = {{._mp_alloc = 0, ._mp_size = 1, ._mp_d = (mp_limb_t[]) {0xb}}};
const ibz_t DEGREE_CHALLENGE = {{._mp_alloc = 0, ._mp_size = 4, ._mp_d = (mp_limb_t[]) {0x0,0x4000000000000000,0xe53f660a761fd41e,0x316bf9bfcc0af2f}}};
#endif
