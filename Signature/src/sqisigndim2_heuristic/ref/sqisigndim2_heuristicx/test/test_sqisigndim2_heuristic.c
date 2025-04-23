#include <rng.h>
#include <stdio.h>
#include <ec.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>

#include "test_sqisigndim2_heuristic.h"
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
    uint64_t t0, t1;
    clock_t t;

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

        t = tic();
        t0 = rdtsc();
        int check = protocols_verif(&sig, &pk, msg, 32);
        if (!check) {
            printf("verif failed ! \n");
        }
        t1 = rdtsc();
        TOC(t, "Verification");
        printf("\x1b[34mverif   took %'" PRIu64 " cycles\x1b[0m\n", t1 - t0);

        printf(" \x1b[35mfull\x1b[0m signature was: %s\n\n",
               check ? "\x1b[32mvalid\x1b[0m" : "\x1b[31minvalid\x1b[0m");
    }

    // TOC(t, "protocols_sign");
    // float ms = (1000. * (float) (clock() - t) / CLOCKS_PER_SEC);
    // printf("average signing time [%.2f ms]\n", (float) (ms/repeat));

    // clock_t t_verif = tic();
    // for (int i = 0; i < repeat; ++i)
    // {

    // }
    // TOC(t_verif,"protocols_verif");
    // ms = (1000. * (float) (clock() - t_verif) / CLOCKS_PER_SEC);
    // printf("average verif time [%.2f ms]\n", (float) (ms/repeat));

    if (bench) {
        float ms;

        printf("\n\nBenchmarking signatures\n");
        t = tic();
        t0 = rdtsc();
        for (int i = 0; i < bench; ++i) {
            protocols_keygen(&pk, &sk);
        }
        t1 = rdtsc();
        // ms = tac();
        ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
        printf("Average keygen time [%.2f ms]\n", (float)(ms / bench));
        printf("\x1b[34mAvg keygen: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

        t = tic();
        t0 = rdtsc();
        for (int i = 0; i < bench; ++i) {
            int val = protocols_sign(&sig, &pk, &sk, msg, 32, 0);
        }
        t1 = rdtsc();
        // ms = tac();
        ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
        printf("Average signature time [%.2f ms]\n", (float)(ms / bench));
        printf("\x1b[34mAvg signature: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

        t = tic();
        t0 = rdtsc();
        for (int i = 0; i < bench; ++i) {
            int check = protocols_verif(&sig, &pk, msg, 32);
            if (!check) {
                printf("verif failed ! \n");
            }
        }
        t1 = rdtsc();
        // ms = tac();
        ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
        printf("Average verification time [%.2f ms]\n", (float)(ms / bench));
        printf("\x1b[34mAvg verification: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);
    }

    public_key_finalize(&pk);
    secret_key_finalize(&sk);
    secret_sig_finalize(&sig);

    return res;
}

// run all tests in module
int
main(int argc, char **argv)
{
    int res = 1;

    randombytes_init((unsigned char *)"some", (unsigned char *)"string", 128);

    // printf("\nRunning encoding tests\n");
    // res &= test_encode();

    printf("\nRunning sqisigndim2_heuristic tests\n \n");

    int runs = 0;
    if (argc > 1)
        runs = atoi(argv[1]);

    res &= test_sqisign(3, runs);

    if (!res) {
        printf("\nSome tests failed!\n");
    } else {
        printf("All tests passed!\n");
    }
    return (!res);
}
