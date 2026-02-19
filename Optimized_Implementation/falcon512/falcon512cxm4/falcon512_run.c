#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "api.h"   // NIST Falcon-512 API

#define NTESTS 100  // Number of iterations for stable benchmarking

/* High-resolution timer (nanoseconds) */
static inline uint64_t get_time_ns() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ((uint64_t)ts.tv_sec * 1000000000ULL) + ts.tv_nsec;
}

int main(void) {

    /* ---- Key buffers (from NIST API) ---- */
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];

    /* ---- Message ---- */
    const char *msg = "Hello Falcon-512 on Raspberry Pi";
    unsigned long long mlen = (unsigned long long)strlen(msg);

    /* ---- Signature buffers (sig || message) ---- */
    unsigned char *sm = (unsigned char *)malloc(CRYPTO_BYTES + mlen);
    unsigned char *m2 = (unsigned char *)malloc(CRYPTO_BYTES + mlen);

    if (!sm || !m2) {
        fprintf(stderr, "Memory allocation failed\n");
        free(sm);
        free(m2);
        return 1;
    }

    unsigned long long smlen = 0;
    unsigned long long m2len = 0;

    /* ---- Timing accumulators ---- */
    uint64_t start, end;
    uint64_t keygen_time = 0;
    uint64_t sign_time = 0;
    uint64_t verify_time = 0;

    printf("========================================\n");
    printf(" Falcon-512 Benchmark (Raspberry Pi ARM)\n");
    printf(" Iterations: %d\n", NTESTS);
    printf("========================================\n\n");

    for (int i = 0; i < NTESTS; i++) {

        /* ---- Key Generation ---- */
        start = get_time_ns();
        if (crypto_sign_keypair(pk, sk) != 0) {
            fprintf(stderr, "Keypair generation failed\n");
            free(sm);
            free(m2);
            return 1;
        }
        end = get_time_ns();
        keygen_time += (end - start);

        /* ---- Signing ---- */
        start = get_time_ns();
        if (crypto_sign(sm, &smlen,
                        (const unsigned char *)msg, mlen,
                        sk) != 0) {
            fprintf(stderr, "Signing failed\n");
            free(sm);
            free(m2);
            return 1;
        }
        end = get_time_ns();
        sign_time += (end - start);

        /* ---- Verification ---- */
        start = get_time_ns();
        if (crypto_sign_open(m2, &m2len,
                             sm, smlen,
                             pk) != 0) {
            fprintf(stderr, "Verification failed\n");
            free(sm);
            free(m2);
            return 2;
        }
        end = get_time_ns();
        verify_time += (end - start);

        /* ---- Integrity Check ---- */
        if (m2len != mlen || memcmp(m2, msg, mlen) != 0) {
            fprintf(stderr, "Message mismatch after verification\n");
            free(sm);
            free(m2);
            return 3;
        }
    }

    /* ---- Compute averages (nanoseconds) ---- */
    uint64_t avg_keygen = keygen_time / NTESTS;
    uint64_t avg_sign   = sign_time / NTESTS;
    uint64_t avg_verify = verify_time / NTESTS;

    printf("========== AVERAGE TIME (nanoseconds) ==========\n");
    printf("KeyGen  : %lu ns\n", avg_keygen);
    printf("Sign    : %lu ns\n", avg_sign);
    printf("Verify  : %lu ns\n", avg_verify);

    printf("\n========== ADDITIONAL METRICS ==========\n");
    printf("Signature size (sig + msg): %llu bytes\n", smlen);
    printf("Public key size           : %d bytes\n", CRYPTO_PUBLICKEYBYTES);
    printf("Secret key size           : %d bytes\n", CRYPTO_SECRETKEYBYTES);

    free(sm);
    free(m2);
    return 0;
}
