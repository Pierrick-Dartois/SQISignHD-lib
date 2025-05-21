#include <rng.h>
#include <stdio.h>
#include <ec.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>
#include <gf_constants.h>
#include <fips202.h>
#include <unistd.h>

#include "test_sqisignhd.h"
#include <tools.h>

static inline int64_t
cpucycles(void)
{
#if (defined(TARGET_AMD64) || defined(TARGET_X86))
    unsigned int hi, lo;

    asm volatile("rdtsc" : "=a"(lo), "=d"(hi));
    return ((int64_t)lo) | (((int64_t)hi) << 32);
#elif (defined(TARGET_S390X))
    uint64_t tod;
    asm volatile("stckf %0\n" : "=Q"(tod) : : "cc");
    return (tod * 1000 / 4096);
#else
    struct timespec time;
    clock_gettime(CLOCK_REALTIME, &time);
    return (int64_t)(time.tv_sec * 1e9 + time.tv_nsec);
#endif
}

static __inline__ uint64_t
rdtsc(void)
{
    return (uint64_t)cpucycles();
}

int
test_sqisign(int repeat, uint64_t bench)
{
    int res = 1;

    int num_sig = 10;

    public_key_t pk;
    secret_key_t sk;
    signature_t sig;
    unsigned char msg[32] = { 0 };

    public_key_init(&pk);
    secret_key_init(&sk);
    secret_sig_init(&(sig));

    // printf("Printing details of first signature\n");
    // int val = protocols_sign(&sig, &pk, &sk, msg, 32, 1);
    setlocale(LC_NUMERIC, "");
    uint64_t t0, t1, dt_kg, dt_sig, dt_ms_kg, dt_ms_sig;
    clock_t t;

    //unsigned char *buf = malloc(FP2_ENCODED_BYTES + FP2_ENCODED_BYTES + 32);
    //{
        //fp2_encode(buf, &NQR_TABLE[0]);
        //fp2_encode(buf + FP2_ENCODED_BYTES, &NQR_TABLE[1]); // TODO use defined constant
        //memcpy(buf + FP2_ENCODED_BYTES + FP2_ENCODED_BYTES,
               //msg,
               //32); // TODO use defined constant
    //}
    //for(int i=0;i<FP2_ENCODED_BYTES + FP2_ENCODED_BYTES + 32;i++){
        //printf("%x",buf[i]);
    //}
    //printf("\n");

    //digit_t digits[NWORDS_FIELD];

    // FIXME should use SHAKE128 for smaller parameter sets?
    //  TODO we want to use a bit differently (first we hash first half and then derive the
    //  second half)
    //SHAKE256((void *)digits, sizeof(digits), buf, FP2_ENCODED_BYTES + FP2_ENCODED_BYTES + 32);

    //ibz_t scalars;
    //ibz_init(&scalars);
    //ibz_copy_digit_array(&scalars, digits);
    //ibz_printf("%Zx\n",scalars);

    //printf("\nNQR_TABLE\n");
    //for (int i = 0; i < 20; ++i) {
        //fp2_print("",&NQR_TABLE[i]);
    //}
    //printf("\nZ_NQR_TABLE\n");
    //for (int i = 0; i < 20; ++i) {
        //fp2_print("",&Z_NQR_TABLE[i]);
    //}

    printf("\n\nTesting signatures\n");
    for (int i = 0; i < repeat; ++i) {
        printf("#%d \n", i);
        t = tic();
        t0 = rdtsc();
        protocols_keygen(&pk, &sk);
        t1 = rdtsc();
        TOC(t, "Keygen");
        printf("\x1b[34mkeygen  took %'" PRIu64 " cycles\x1b[0m\n", t1 - t0);

        t = tic();
        t0 = rdtsc();
        int val = protocols_sign(&sig, &pk, &sk, msg, 32, 0);
        t1 = rdtsc();
        TOC(t, "Signing");
        printf("\x1b[34msigning took %'" PRIu64 " cycles\x1b[0m\n", t1 - t0);
    }

    float ms;

    char buffer[]="";
    getcwd(buffer, 100);
    printf("%s\n", buffer);

    FILE *p_file_pk;
    FILE *p_file_sign;
    p_file_pk = fopen("Public_keys.txt", "w");
    p_file_sign = fopen("Signatures.txt", "w");

    printf("\n\nBenchmarking signatures\n");
    dt_kg = 0;
    dt_sig = 0;
    dt_ms_kg = 0;
    dt_ms_sig = 0;
    for (int i = 0; i < bench; ++i) {
        t = tic();
        t0 = rdtsc();
        protocols_keygen(&pk, &sk);
        t1 = rdtsc();
        ms = (1000. * (float)(clock()-t) / CLOCKS_PER_SEC);
        dt_kg = dt_kg+t1-t0;
        dt_ms_kg = dt_ms_kg+ms;
        fprint_public_key(p_file_pk, &pk);

        t = tic();
        t0 = rdtsc();
        int val = protocols_sign(&sig, &pk, &sk, msg, 32, 0);
        t1 = rdtsc();
        ms = (1000. * (float)(clock()-t) / CLOCKS_PER_SEC);
        dt_sig = dt_sig+t1-t0;
        dt_ms_sig = dt_ms_sig+ms;
        fprint_signature(p_file_sign, &sig);

    }
    
    printf("Average keygen time [%.2f ms]\n", (float)(dt_ms_kg / bench));
    printf("\x1b[34mAvg keygen: %'" PRIu64 " cycles\x1b[0m\n", dt_kg / bench);

    printf("Average signature time [%.2f ms]\n", (float)(dt_ms_sig / bench));
    printf("\x1b[34mAvg signature: %'" PRIu64 " cycles\x1b[0m\n", dt_sig / bench);

    public_key_finalize(&pk);
    secret_key_finalize(&sk);
    secret_sig_finalize(&sig);

    return res;
}

// run all tests in module
int
main()
{
    int res = 1;

    randombytes_init((unsigned char *)"some", (unsigned char *)"string", 128);

    // printf("\nRunning encoding tests\n");
    // res &= test_encode();

    printf("\nRunning sqisignhd tests\n \n");

    res &= test_sqisign(3, 100);

    if (!res) {
        printf("\nSome tests failed!\n");
    } else {
        printf("All tests passed!\n");
    }
    return (!res);
}
