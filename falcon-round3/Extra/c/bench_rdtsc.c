/*
 * RDTSC profiling driver for Falcon sign_dyn.
 * Build via: make bench_rdtsc
 * Run  via:  ./bench_rdtsc [logn]   (default logn=9 => Falcon-512)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "falcon.h"

#define NTRIALS  10

int
main(int argc, char *argv[])
{
	int logn = (argc >= 2) ? atoi(argv[1]) : 9;
	if (logn < 1 || logn > 10) {
		fprintf(stderr, "logn must be 1..10 (9=Falcon-512, 10=Falcon-1024)\n");
		return 1;
	}

	size_t pk_len  = FALCON_PUBKEY_SIZE(logn);
	size_t sk_len  = FALCON_PRIVKEY_SIZE(logn);
	size_t sig_len = FALCON_SIG_COMPRESSED_MAXSIZE(logn);
	size_t tmp_kg  = FALCON_TMPSIZE_KEYGEN(logn);
	size_t tmp_sg  = FALCON_TMPSIZE_SIGNDYN(logn);
	size_t tmp_len = tmp_kg > tmp_sg ? tmp_kg : tmp_sg;

	uint8_t *pk  = malloc(pk_len);
	uint8_t *sk  = malloc(sk_len);
	uint8_t *sig = malloc(sig_len);
	uint8_t *tmp = malloc(tmp_len);
	if (!pk || !sk || !sig || !tmp) {
		fprintf(stderr, "malloc failed\n");
		return 1;
	}

	shake256_context rng;
	shake256_init_prng_from_system(&rng);

	printf("Falcon-%d RDTSC profiling (%d trials)\n",
		1 << logn, NTRIALS);
	printf("Generating keypair...\n");
	if (falcon_keygen_make(&rng, logn,
		sk, sk_len, pk, pk_len, tmp, tmp_len) != 0)
	{
		fprintf(stderr, "keygen failed\n");
		return 1;
	}

	static const uint8_t msg[32] = "rdtsc-profiling-test-message-xx";

	printf("Signing (each trial prints one report):\n\n");
	for (int i = 0; i < NTRIALS; i++) {
		size_t sl = sig_len;
		if (falcon_sign_dyn(&rng,
			sig, &sl, FALCON_SIG_COMPRESSED,
			sk, sk_len,
			msg, sizeof msg,
			tmp, tmp_len) != 0)
		{
			fprintf(stderr, "sign failed\n");
			return 1;
		}
	}

	free(pk); free(sk); free(sig); free(tmp);
	return 0;
}
