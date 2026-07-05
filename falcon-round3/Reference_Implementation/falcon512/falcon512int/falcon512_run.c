#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "api.h"   // NIST Falcon-512 API

int main(void) {
    /* ---- Key buffers (fixed-size from api.h) ---- */
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];

    /* ---- Message ---- */
    const char *msg = "Hello Falcon-512";
    unsigned long long mlen = (unsigned long long)strlen(msg);

    /* ---- Signature buffers (NIST format: sig || message) ---- */
    unsigned char *sm = malloc(CRYPTO_BYTES + mlen);
    unsigned char *m2 = malloc(CRYPTO_BYTES + mlen);
    if (!sm || !m2) {
        fprintf(stderr, "Allocation failed\n");
        free(sm);
        free(m2);
        return 1;
    }

    unsigned long long smlen = 0;
    unsigned long long m2len = 0;

    /* ---- Key generation ---- */
    if (crypto_sign_keypair(pk, sk) != 0) {
        fprintf(stderr, "Keypair generation failed\n");
        free(sm);
        free(m2);
        return 1;
    }
    printf("✔ Keypair generated (Falcon-512)\n");

    /* ---- Sign ---- */
    if (crypto_sign(sm, &smlen,
                    (const unsigned char *)msg, mlen,
                    sk) != 0) {
        fprintf(stderr, "Signing failed\n");
        free(sm);
        free(m2);
        return 1;
    }
    printf("✔ Message signed (signature+msg = %llu bytes)\n", smlen);

    /* ---- Verify ---- */
    if (crypto_sign_open(m2, &m2len,
                         sm, smlen,
                         pk) != 0) {
        fprintf(stderr, "✘ Verification failed\n");
        free(sm);
        free(m2);
        return 2;
    }

    /* ---- Check recovered message ---- */
    if (m2len != mlen || memcmp(m2, msg, mlen) != 0) {
        fprintf(stderr, "✘ Message mismatch after verification\n");
        free(sm);
        free(m2);
        return 3;
    }

    printf("✔ Signature verified successfully\n");

    free(sm);
    free(m2);
    return 0;
}
