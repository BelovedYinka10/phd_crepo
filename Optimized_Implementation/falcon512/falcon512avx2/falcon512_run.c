#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <x86intrin.h>   // For __rdtsc() CPU cycle counter

#include "api.h"   // NIST Falcon-512 API

#define NTESTS 100  // Number of iterations for stable benchmarking

int main(void) {
    /* ---- Key buffers (fixed-size from api.h) ---- */
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];

    /* ---- Message ---- */
    const char *msg = "Hello Falcon-512";
    unsigned long long mlen = (unsigned long long)strlen(msg);

    /* ---- Signature buffers (NIST format: sig || message) ---- */
    unsigned char *sm = (unsigned char *)malloc(CRYPTO_BYTES + mlen);
    unsigned char *m2 = (unsigned char *)malloc(CRYPTO_BYTES + mlen);

    if (!sm || !m2) {
        fprintf(stderr, "Allocation failed\n");
        free(sm);
        free(m2);
        return 1;
    }

    unsigned long long smlen = 0;
    unsigned long long m2len = 0;

    /* ---- Cycle counters ---- */
    unsigned long long start, end;
    unsigned long long keygen_cycles = 0;
    unsigned long long sign_cycles = 0;
    unsigned long long verify_cycles = 0;

    printf("========================================\n");
    printf(" Falcon-512 CPU Cycle Benchmark (Reference)\n");
    printf(" Iterations: %d\n", NTESTS);
    printf("========================================\n\n");

    for (int i = 0; i < NTESTS; i++) {

        /* ---- Key Generation (Cycle Measured) ---- */
        start = __rdtsc();
        if (crypto_sign_keypair(pk, sk) != 0) {
            fprintf(stderr, "Keypair generation failed\n");
            free(sm);
            free(m2);
            return 1;
        }
        end = __rdtsc();
        keygen_cycles += (end - start);

        /* ---- Sign (Cycle Measured) ---- */
        start = __rdtsc();
        if (crypto_sign(sm, &smlen,
                        (const unsigned char *)msg, mlen,
                        sk) != 0) {
            fprintf(stderr, "Signing failed\n");
            free(sm);
            free(m2);
            return 1;
        }
        end = __rdtsc();
        sign_cycles += (end - start);

        /* ---- Verify (Cycle Measured) ---- */
        start = __rdtsc();
        if (crypto_sign_open(m2, &m2len,
                             sm, smlen,
                             pk) != 0) {
            fprintf(stderr, "Verification failed\n");
            free(sm);
            free(m2);
            return 2;
        }
        end = __rdtsc();
        verify_cycles += (end - start);

        /* ---- Sanity check (message integrity) ---- */
        if (m2len != mlen || memcmp(m2, msg, mlen) != 0) {
            fprintf(stderr, "Message mismatch after verification\n");
            free(sm);
            free(m2);
            return 3;
        }
    }

    /* ---- Compute Averages ---- */
    unsigned long long avg_keygen = keygen_cycles / NTESTS;
    unsigned long long avg_sign   = sign_cycles / NTESTS;
    unsigned long long avg_verify = verify_cycles / NTESTS;

    printf("========== AVERAGE CPU CYCLES ==========\n");
    printf("KeyGen  : %llu cycles\n", avg_keygen);
    printf("Sign    : %llu cycles\n", avg_sign);
    printf("Verify  : %llu cycles\n", avg_verify);

    printf("\n========== ADDITIONAL METRICS ==========\n");
    printf("Signature size (sig + msg): %llu bytes\n", smlen);
    printf("Public key size           : %d bytes\n", CRYPTO_PUBLICKEYBYTES);
    printf("Secret key size           : %d bytes\n", CRYPTO_SECRETKEYBYTES);

    free(sm);
    free(m2);
    return 0;
}
