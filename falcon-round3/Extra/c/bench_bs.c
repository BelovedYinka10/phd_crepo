/*
 * bench_bs.c — BaseSampler's share of SIGNING time (scalar C reference).
 * =====================================================================
 *
 * Answers the paper's framing ("BaseSampler ~ 30% of signing") against the
 * real C implementation rather than the Python proxy.
 *
 * Method (no perturbation of the timed signing path beyond a 1-cycle counter):
 *   1. keygen + expand a Falcon-512 (and -1024) key.
 *   2. Time falcon_sign_tree() over many iterations  -> t_sign per signature.
 *      A counter inside gaussian0_sampler records BaseSampler calls/signature.
 *   3. Microbenchmark gaussian0_sampler() alone        -> t_g0 per call.
 *   4. share = (calls_per_sig * t_g0) / t_sign.
 *
 * Build (instrumentation flag turns on the counter + microbench helper):
 *   make bench_bs            # see Makefile target added for this
 * or:
 *   clang -O3 -DFALCON_BS_COUNT -o bench_bs bench_bs.c codec.o common.o \
 *         falcon.o fft.o fpr.o keygen.o rng.o shake.o sign.o vrfy.o
 *
 * NOTE: this is the SCALAR path on arm64 (no AVX2 — that's x86 only).  It is
 * an honest baseline, not the paper's AVX2 figure.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "falcon.h"

/* Defined in sign.c under -DFALCON_BS_COUNT. */
extern unsigned long long falcon_bs_count;
extern uint64_t falcon_bench_gaussian0(unsigned long num);

static double
now(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static void
run_degree(unsigned logn, unsigned long sign_iters, unsigned long g0_iters)
{
	size_t n = (size_t)1 << logn;
	shake256_context rng;
	uint8_t *sk, *pk, *esk, *sig, *tmp;
	size_t sig_len, tmp_len;
	unsigned long i;
	double t0, t_sign, t_g0;
	double sign_per, g0_per, calls_per_sig, bs_time_per_sig, share;
	uint64_t sink;

	tmp_len = FALCON_TMPSIZE_SIGNTREE(logn);
	if (FALCON_TMPSIZE_KEYGEN(logn) > tmp_len)
		tmp_len = FALCON_TMPSIZE_KEYGEN(logn);
	if (FALCON_TMPSIZE_EXPANDPRIV(logn) > tmp_len)
		tmp_len = FALCON_TMPSIZE_EXPANDPRIV(logn);

	sk  = malloc(FALCON_PRIVKEY_SIZE(logn));
	pk  = malloc(FALCON_PUBKEY_SIZE(logn));
	esk = malloc(FALCON_EXPANDEDKEY_SIZE(logn));
	sig = malloc(FALCON_SIG_COMPRESSED_MAXSIZE(logn));
	tmp = malloc(tmp_len);
	if (!sk || !pk || !esk || !sig || !tmp) { fprintf(stderr, "OOM\n"); exit(1); }

	if (shake256_init_prng_from_system(&rng) != 0) { fprintf(stderr, "rng\n"); exit(1); }

	if (falcon_keygen_make(&rng, logn, sk, FALCON_PRIVKEY_SIZE(logn),
		pk, FALCON_PUBKEY_SIZE(logn), tmp, tmp_len) != 0) { fprintf(stderr, "keygen\n"); exit(1); }
	if (falcon_expand_privkey(esk, FALCON_EXPANDEDKEY_SIZE(logn),
		sk, FALCON_PRIVKEY_SIZE(logn), tmp, tmp_len) != 0) { fprintf(stderr, "expand\n"); exit(1); }

	/* warm up */
	sig_len = FALCON_SIG_COMPRESSED_MAXSIZE(logn);
	falcon_sign_tree(&rng, sig, &sig_len, FALCON_SIG_COMPRESSED, esk, "data", 4, tmp, tmp_len);

	/* ---- timed signing, with BaseSampler call counting ---- */
	falcon_bs_count = 0;
	t0 = now();
	for (i = 0; i < sign_iters; i++) {
		sig_len = FALCON_SIG_COMPRESSED_MAXSIZE(logn);
		if (falcon_sign_tree(&rng, sig, &sig_len, FALCON_SIG_COMPRESSED,
			esk, "data", 4, tmp, tmp_len) != 0) { fprintf(stderr, "sign\n"); exit(1); }
	}
	t_sign = now() - t0;
	calls_per_sig = (double)falcon_bs_count / (double)sign_iters;

	/* ---- microbenchmark gaussian0_sampler alone ---- */
	sink = falcon_bench_gaussian0(g0_iters / 10);    /* warm up */
	t0 = now();
	sink += falcon_bench_gaussian0(g0_iters);
	t_g0 = now() - t0;

	sign_per = t_sign / (double)sign_iters;
	g0_per   = t_g0 / (double)g0_iters;
	bs_time_per_sig = calls_per_sig * g0_per;
	share    = bs_time_per_sig / sign_per;

	printf("Falcon-%zu  (scalar C, %lu sign iters, %lu g0 iters)\n",
		n, sign_iters, g0_iters);
	printf("  sign_tree per signature   : %9.2f us\n", sign_per * 1e6);
	printf("  BaseSampler calls / sig    : %9.1f\n", calls_per_sig);
	printf("  gaussian0_sampler per call : %9.2f ns\n", g0_per * 1e9);
	printf("  BaseSampler time / sig     : %9.2f us\n", bs_time_per_sig * 1e6);
	printf("  ------------------------------------------\n");
	printf("  BaseSampler share of signing : %6.1f %%   (sink=%llu)\n",
		share * 100.0, (unsigned long long)sink);
	printf("\n");

	free(sk); free(pk); free(esk); free(sig); free(tmp);
}

int
main(void)
{
	printf("BaseSampler share of signing — scalar C reference (arm64, no AVX2)\n");
	printf("==================================================================\n\n");
	run_degree(9,  20000, 20000000UL);   /* Falcon-512 */
	run_degree(10, 10000, 20000000UL);   /* Falcon-1024 */
	printf("Note: scalar path (AVX2 is x86-only). Honest baseline, not the\n");
	printf("paper's AVX2 number. Counter overhead (~1 cycle/call) is included\n");
	printf("in both t_sign and t_g0, so it cancels in the ratio.\n");
	return 0;
}
