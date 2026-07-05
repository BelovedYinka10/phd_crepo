/*
Misaligned-buffer test for XKCP.

Feeds deliberately misaligned input/output/key buffers (offsets 1..7 from
an aligned arena) through the public API and checks that the results match
the same computation on aligned buffers.

Two failure modes it detects in unpatched code:
- Built with -fsanitize=alignment (clang), every misaligned lane access in
  the library is reported as a runtime error, deterministically, on any
  architecture -- including x86-64 and AArch64, whose hardware would
  otherwise tolerate the access silently.
- On strict-alignment hardware (e.g. 32-bit ARM), the unpatched library
  crashes with a bus error outright.

Exit code 0 = all aligned/misaligned result pairs match.
*/

#include <stdio.h>
#include <string.h>
#include "config.h"
#include "align.h"
#include "SimpleFIPS202.h"
#include "KeccakHash.h"
#include "SP800-185.h"
#ifdef XKCP_has_Xoodyak
#include "Xoodyak.h"
#endif

#define INLEN  1000
#define OUTLEN 136     /* > one rate block of SHAKE256, exercises ExtractLanes */
#define KEYLEN 32

ALIGN(64) static unsigned char arena_in[INLEN + 64];
ALIGN(64) static unsigned char arena_out[OUTLEN + 64];
ALIGN(64) static unsigned char arena_key[KEYLEN + 64];

static int nFailed = 0;

static void check(const char *what, unsigned int offset,
                  const unsigned char *ref, const unsigned char *got, size_t len)
{
    if (memcmp(ref, got, len) != 0) {
        printf("MISMATCH: %s at offset %u\n", what, offset);
        nFailed++;
    }
}

int main(void)
{
    unsigned char refSHA3[32], refSHAKE[OUTLEN], refKMAC[32], refInc[32];
#ifdef XKCP_has_Xoodyak
    unsigned char refXoodyak[32];
#endif
    unsigned int offset, i;

    for (i = 0; i < INLEN; i++)
        arena_in[i] = (unsigned char)(i * 7 + 1);
    for (i = 0; i < KEYLEN; i++)
        arena_key[i] = (unsigned char)(0x40 + i);

    /* Reference results on aligned buffers (offset 0). */
    SHA3_256(refSHA3, arena_in, INLEN);
    SHAKE256(refSHAKE, OUTLEN, arena_in, INLEN);
    KMAC128(arena_key, KEYLEN*8, arena_in, INLEN*8, refKMAC, 256, NULL, 0);
    {
        Keccak_HashInstance h;
        Keccak_HashInitialize_SHA3_256(&h);
        Keccak_HashUpdate(&h, arena_in, 7*8);
        Keccak_HashUpdate(&h, arena_in + 7, (INLEN - 7)*8);
        Keccak_HashFinal(&h, refInc);
    }
#ifdef XKCP_has_Xoodyak
    {
        Xoodyak_Instance x;
        Xoodyak_Initialize(&x, NULL, 0, NULL, 0, NULL, 0);
        Xoodyak_Absorb(&x, arena_in, INLEN);
        Xoodyak_Squeeze(&x, refXoodyak, 32);
    }
#endif

    for (offset = 1; offset < 8; offset++) {
        unsigned char *in  = arena_in  + offset;   /* misaligned input  */
        unsigned char *out = arena_out + offset;   /* misaligned output */
        unsigned char *key = arena_key + offset;   /* misaligned key    */
        unsigned char aref[64];

        memmove(in, arena_in, INLEN);              /* same content, shifted */
        memmove(key, arena_key, KEYLEN);

        /* one-shot hash, misaligned input and output */
        SHA3_256(out, in, INLEN);
        check("SHA3-256", offset, refSHA3, out, 32);

        /* XOF with long misaligned output (ExtractLanes path) */
        SHAKE256(out, OUTLEN, in, INLEN);
        check("SHAKE256", offset, refSHAKE, out, OUTLEN);

        /* MAC, all three buffers misaligned */
        KMAC128(key, KEYLEN*8, in, INLEN*8, out, 256, NULL, 0);
        check("KMAC128", offset, refKMAC, out, 32);

        /* incremental hashing with a non-lane-multiple first chunk,
           so the second chunk enters AddLanes/FastLoop misaligned */
        {
            Keccak_HashInstance h;
            Keccak_HashInitialize_SHA3_256(&h);
            Keccak_HashUpdate(&h, in, 7*8);
            Keccak_HashUpdate(&h, in + 7, (INLEN - 7)*8);
            Keccak_HashFinal(&h, aref);
            check("SHA3-256 incremental", offset, refInc, aref, 32);
        }

#ifdef XKCP_has_Xoodyak
        {
            Xoodyak_Instance x;
            Xoodyak_Initialize(&x, NULL, 0, NULL, 0, NULL, 0);
            Xoodyak_Absorb(&x, in, INLEN);
            Xoodyak_Squeeze(&x, out, 32);
            check("Xoodyak", offset, refXoodyak, out, 32);
        }
#endif

        memmove(arena_in, in, INLEN);              /* restore for next round */
        memmove(arena_key, key, KEYLEN);
    }

    if (nFailed == 0)
        printf("All misaligned-buffer results match the aligned references.\n");
    return nFailed != 0;
}
