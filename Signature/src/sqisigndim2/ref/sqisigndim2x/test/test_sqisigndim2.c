#include <rng.h>
#include <stdio.h>
#include <ec.h>
#include <inttypes.h>
#include <locale.h>
#include <time.h>

#include "test_sqisigndim2.h"
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

void
bench_fp2_operations(int repeat)
{

    uint64_t t0, t1;

    fp2_t a = BASIS_EVEN.P.x;
    fp2_t b = BASIS_EVEN.PmQ.x;
    t0 = rdtsc();
    for (int i = 0; i < repeat; i++) {
        fp2_add(&a, &b, &a);
    }
    t1 = rdtsc();
    printf("\x1b[34m fp2_add %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / repeat);

    a = BASIS_EVEN.P.x;
    b = BASIS_EVEN.Q.x;

    t0 = rdtsc();
    for (int i = 0; i < repeat; i++) {
        fp2_sqr(&a, &a);
    }
    t1 = rdtsc();
    printf("\x1b[34m fp2_sqr %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / repeat);

    a = BASIS_EVEN.P.x;
    b = BASIS_EVEN.PmQ.x;

    t0 = rdtsc();
    for (int i = 0; i < repeat; i++) {
        fp2_mul(&a, &b, &a);
    }
    t1 = rdtsc();
    printf("\x1b[34m fp2_mul %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / repeat);

    t0 = rdtsc();
    for (int i = 0; i < repeat; i++) {
        fp2_inv(&a);
        fp2_mul(&a, &b, &a);
    }
    t1 = rdtsc();
    printf("\x1b[34m fp2_inv %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / repeat);
}

int
test_sqisign(int repeat, uint64_t bench)
{
    int res = 1;

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

    if (bench) {
        float ms;

        public_key_t pks[bench];
        secret_key_t sks[bench];
        signature_t sigs[bench];
        for (int i = 0; i < bench; i++) {
            public_key_init(&(pks[i]));
            secret_key_init(&(sks[i]));
            secret_sig_init(&(sigs[i]));
        }

        printf("\n\nBenchmarking signatures\n");
        t = tic();
        t0 = rdtsc();
        for (int i = 0; i < bench; ++i) {
            protocols_keygen(&pks[i], &sks[i]);
        }
        t1 = rdtsc();
        // ms = tac();
        ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
        printf("Average keygen time [%.2f ms]\n", (float)(ms / bench));
        printf("\x1b[34mAvg keygen: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

        t = tic();
        t0 = rdtsc();
        for (int i = 0; i < bench; ++i) {
            int val = protocols_sign(&(sigs[i]), &(pks[i]), &(sks[i]), msg, 32, 0);
        }
        t1 = rdtsc();
        // ms = tac();
        ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
        printf("Average signature time [%.2f ms]\n", (float)(ms / bench));
        printf("\x1b[34mAvg signature: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

        t = tic();
        t0 = rdtsc();
        for (int i = 0; i < bench; ++i) {
            int check = protocols_verif(&(sigs[i]), &(pks[i]), msg, 32);
            if (!check) {
                printf("verif failed ! \n");
            }
        }
        t1 = rdtsc();
        // ms = tac();
        ms = (1000. * (float)(clock() - t) / CLOCKS_PER_SEC);
        printf("Average verification time [%.2f ms]\n", (float)(ms / bench));
        printf("\x1b[34mAvg verification: %'" PRIu64 " cycles\x1b[0m\n", (t1 - t0) / bench);

        for (int i = 0; i < bench; i++) {
            public_key_init(&(pks[i]));
            secret_key_init(&(sks[i]));
            secret_sig_init(&(sigs[i]));
        }
    }

    public_key_finalize(&pk);
    secret_key_finalize(&sk);
    secret_sig_finalize(&sig);

    return res;
}

void
test_LLL()
{

    // generatting a random ideal
    ibz_t n;
    quat_left_ideal_t ideal;
    ibz_mat_4x4_t gram, reduced;

    ibz_init(&n);
    ibz_mat_4x4_init(&gram);
    ibz_mat_4x4_init(&reduced);
    quat_left_ideal_init(&ideal);

    generate_random_prime(&n, 1, ibz_bitsize(&QUATALG_PINFTY.p) / 2);
    sampling_random_ideal_O0(&ideal, &n, 1);

    quat_lideal_reduce_basis(&reduced, &gram, &ideal, &QUATALG_PINFTY);

    printf("for a %d-bit norm input, the bitsize of the reduced basis is %d %d %d %d \n",
           ibz_bitsize(&n),
           ibz_bitsize(&gram[0][0]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[1][1]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[2][2]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[3][3]) - ibz_bitsize(&n));

    generate_random_prime(&n, 1, ibz_bitsize(&QUATALG_PINFTY.p));
    sampling_random_ideal_O0(&ideal, &n, 1);

    quat_lideal_reduce_basis(&reduced, &gram, &ideal, &QUATALG_PINFTY);

    printf("for a %d-bit norm input, the bitsize of the reduced basis is %d %d %d %d \n",
           ibz_bitsize(&n),
           ibz_bitsize(&gram[0][0]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[1][1]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[2][2]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[3][3]) - ibz_bitsize(&n));

    generate_random_prime(&n, 1, 2 * ibz_bitsize(&QUATALG_PINFTY.p));
    sampling_random_ideal_O0(&ideal, &n, 1);

    quat_lideal_reduce_basis(&reduced, &gram, &ideal, &QUATALG_PINFTY);

    printf("for a %d-bit norm input, the bitsize of the reduced basis is %d %d %d %d \n",
           ibz_bitsize(&n),
           ibz_bitsize(&gram[0][0]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[1][1]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[2][2]) - ibz_bitsize(&n),
           ibz_bitsize(&gram[3][3]) - ibz_bitsize(&n));

    ibz_finalize(&n);
    quat_left_ideal_finalize(&ideal);
    ibz_mat_4x4_finalize(&gram);
    ibz_mat_4x4_finalize(&reduced);
}

// run all tests in module
int
main(int argc, char **argv)
{
    int res = 1;

    randombytes_init((unsigned char *)"some", (unsigned char *)"string", 128);

    // printf("\nRunning encoding tests\n");
    // res &= test_encode();
    printf("\nRunning reduction test\n \n");
    test_LLL();

    int runs = 100;
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
