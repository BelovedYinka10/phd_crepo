/*
 * bench_domains.c — RDTSC domain-level profiler for Falcon keygen / sign / verify
 *
 * Prints a table showing cycle cost of each cryptographic domain:
 *
 *   KEYGEN : NTRU solve | compute_public (NTT) | f,g sampling overhead
 *   SIGN   : Gram (B·B*) | LDL build | BaseSampler | berexp | FFT ops | hash+enc
 *   VERIFY : NTT poly_mul | norm check | other
 *
 * Build:
 *   make bench_domains
 *
 * Usage:
 *   ./bench_domains [logn]   (default logn=9 → Falcon-512)
 *   ./bench_domains 10       (Falcon-1024)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "falcon.h"
#include "inner.h"

/* ── RDTSC (serialised) ────────────────────────────────────────────────────── */
#if defined(__x86_64__) || defined(__i386__)
# include <x86intrin.h>
static inline uint64_t tsc_start(void) {
	_mm_mfence(); return __rdtsc();
}
static inline uint64_t tsc_stop(void) {
	uint64_t v; unsigned aux;
	v = __rdtscp(&aux); _mm_mfence(); return v;
}
#elif defined(__aarch64__)
static inline uint64_t tsc_start(void) {
	uint64_t v; __asm__ volatile("isb; mrs %0, cntvct_el0" : "=r"(v)); return v;
}
static inline uint64_t tsc_stop(void)  { return tsc_start(); }
#else
static inline uint64_t tsc_start(void) { return 0; }
static inline uint64_t tsc_stop(void)  { return 0; }
#endif

/* ── Extern globals from instrumented translation units ─────────────────────
 *   keygen_domains.o  (keygen.c  -DFALCON_RDTSC)
 *   sign_domains.o    (sign.c    -DFALCON_RDTSC)
 *   vrfy_domains.o    (vrfy.c    -DFALCON_RDTSC)
 */
extern uint64_t kg_cyc_total, kg_cyc_ntru, kg_cyc_pubkey;
extern unsigned long long kg_n_attempts;

extern uint64_t g_cyc_sign, g_cyc_ffsamp, g_cyc_samplerz;
extern uint64_t g_cyc_gram, g_cyc_ldl_build, g_cyc_basesampler;
extern unsigned long long g_n_samplerz, g_n_attempts;

extern uint64_t vf_cyc_total, vf_cyc_ntt, vf_cyc_norm;

/* ── Stats helpers ──────────────────────────────────────────────────────────── */
#define MAX_TRIALS 100

static int cmp_u64(const void *a, const void *b) {
	uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
	return (x > y) - (x < y);
}

typedef struct { uint64_t med, mn, mx; } Stats;

static Stats compute_stats(uint64_t *arr, int n) {
	Stats s;
	qsort(arr, n, sizeof *arr, cmp_u64);
	s.med = arr[n / 2];
	s.mn  = arr[0];
	s.mx  = arr[n - 1];
	return s;
}

static double cpu_ghz = 2.6;   /* updated from sysctl at startup */

static void get_cpu_freq(void) {
#if defined(__APPLE__)
	FILE *f = popen("sysctl -n hw.cpufrequency_max 2>/dev/null", "r");
	if (f) { uint64_t hz = 0; if (fscanf(f, "%llu", &hz) == 1 && hz > 0) cpu_ghz = hz / 1e9; pclose(f); }
#endif
}

static void fmt_cyc(char *buf, uint64_t c) {
	if (c >= (uint64_t)1e9)      sprintf(buf, "%8.3f Gc", c / 1e9);
	else if (c >= (uint64_t)1e6) sprintf(buf, "%8.3f Mc", c / 1e6);
	else if (c >= (uint64_t)1e3) sprintf(buf, "%8.1f Kc", c / 1e3);
	else                          sprintf(buf, "%8llu  c", (unsigned long long)c);
}

static void fmt_time(char *buf, uint64_t c) {
	double us = c / (cpu_ghz * 1e9) * 1e6;
	if      (us >= 1e6) sprintf(buf, "%8.3f s ", us / 1e6);
	else if (us >= 1e3) sprintf(buf, "%8.3f ms", us / 1e3);
	else                 sprintf(buf, "%8.2f us", us);
}

#define ROW_WIDTH 14
static void print_row(const char *name, Stats s, uint64_t total) {
	char mc[32], tc[32], mn[32], mx[32];
	fmt_cyc(mc, s.med); fmt_time(tc, s.med);
	fmt_cyc(mn, s.mn);  fmt_cyc(mx, s.mx);
	double pct = total ? 100.0 * s.med / total : 0.0;
	printf("  %-22s  %s  %s  %s  %s  %5.1f%%\n",
	       name, mc, tc, mn, mx, pct);
}

static void print_header(void) {
	printf("  %-22s  %8s  %10s  %8s  %8s  %6s\n",
	       "Domain", "Median", "Time", "Min", "Max", "%");
	printf("  %s\n", "----------------------------------------------------------------------"
	                  "-------------------");
}

/* ── Keygen benchmark ───────────────────────────────────────────────────────── */
static void bench_keygen(unsigned logn, int trials) {
	size_t n     = (size_t)1 << logn;
	size_t pklen = FALCON_PUBKEY_SIZE(logn);
	size_t sklen = FALCON_PRIVKEY_SIZE(logn);
	size_t tmplen= FALCON_TMPSIZE_KEYGEN(logn);

	uint8_t *pk  = malloc(pklen);
	uint8_t *sk  = malloc(sklen);
	uint8_t *tmp = malloc(tmplen + 63);
	uint8_t *tmpa= (uint8_t *)(((uintptr_t)tmp + 63) & ~(uintptr_t)63);

	shake256_context rng;
	shake256_init_prng_from_system(&rng);

	uint64_t s_total[MAX_TRIALS], s_ntru[MAX_TRIALS],
	         s_pubkey[MAX_TRIALS], s_other[MAX_TRIALS];

	printf("\n======================================================================\n");
	printf("  KEYGEN  Falcon-%u   (%d trials)\n", 1u << logn, trials);
	printf("======================================================================\n");

	for (int i = 0; i < trials; i++) {
		kg_cyc_total = kg_cyc_ntru = kg_cyc_pubkey = 0;
		kg_n_attempts = 0;

		if (falcon_keygen_make(&rng, logn, sk, sklen, pk, pklen, tmpa, tmplen) != 0) {
			fprintf(stderr, "keygen failed\n"); exit(1);
		}
		s_total[i]  = kg_cyc_total;
		s_ntru[i]   = kg_cyc_ntru;
		s_pubkey[i] = kg_cyc_pubkey;
		s_other[i]  = kg_cyc_total - kg_cyc_ntru - kg_cyc_pubkey;
	}

	Stats st = compute_stats(s_total,  trials);
	Stats sn = compute_stats(s_ntru,   trials);
	Stats sp = compute_stats(s_pubkey, trials);
	Stats so = compute_stats(s_other,  trials);
	uint64_t tot = st.med;

	print_header();
	print_row("NTRU (solve_NTRU)",      sn, tot);
	print_row("NTT  (compute_public)",  sp, tot);
	print_row("f,g sample+norm+other",  so, tot);
	printf("  %-22s  ", "TOTAL");
	{ char mc[32], tc[32]; fmt_cyc(mc, tot); fmt_time(tc, tot);
	  printf("%s  %s\n", mc, tc); }

	printf("\n  *** Bottleneck: [NTRU solve_NTRU]  %.1f%% ***\n", 100.0*sn.med/tot);
	printf("  Avg retry attempts: %.1f per keygen\n",
	       (double)kg_n_attempts / trials);

	free(pk); free(sk); free(tmp);
	(void)n;
}

/* ── Sign benchmark ─────────────────────────────────────────────────────────── */
static void bench_sign(unsigned logn, int trials,
                       const uint8_t *sk, size_t sklen,
                       const uint8_t *msg, size_t mlen) {
	size_t sigmax = FALCON_SIG_COMPRESSED_MAXSIZE(logn);
	size_t tmplen = FALCON_TMPSIZE_SIGNDYN(logn);
	uint8_t *sig = malloc(sigmax);
	uint8_t *tmp = malloc(tmplen + 63);
	uint8_t *tmpa= (uint8_t *)(((uintptr_t)tmp + 63) & ~(uintptr_t)63);

	shake256_context rng;
	shake256_init_prng_from_system(&rng);

	uint64_t s_total[MAX_TRIALS], s_gram[MAX_TRIALS], s_ldl[MAX_TRIALS],
	         s_bs[MAX_TRIALS], s_berexp[MAX_TRIALS],
	         s_fftops[MAX_TRIALS], s_other[MAX_TRIALS];

	int completed = 0, attempts = 0;
	while (completed < trials) {
		attempts++;
		g_cyc_sign = g_cyc_ffsamp = g_cyc_samplerz = 0;
		g_cyc_gram = g_cyc_ldl_build = g_cyc_basesampler = 0;
		g_n_samplerz = g_n_attempts = 0;

		size_t sl = sigmax;
		int rc = falcon_sign_dyn(&rng, sig, &sl, FALCON_SIG_COMPRESSED,
		                         sk, sklen, msg, mlen, tmpa, tmplen);
		if (rc != 0) continue;

		s_total[completed]  = g_cyc_sign;
		s_gram[completed]   = g_cyc_gram;
		s_ldl[completed]    = g_cyc_ldl_build;
		s_bs[completed]     = g_cyc_basesampler;
		s_berexp[completed] = g_cyc_samplerz - g_cyc_basesampler;
		/* FFT ops inside ffSampling = ffsamp - samplerz - ldl_build */
		uint64_t fftops = (g_cyc_ffsamp > g_cyc_samplerz + g_cyc_ldl_build)
		                ? g_cyc_ffsamp - g_cyc_samplerz - g_cyc_ldl_build : 0;
		s_fftops[completed] = fftops;
		/* hash + compress + target setup = sign - gram - ffsamp */
		uint64_t other = (g_cyc_sign > g_cyc_gram + g_cyc_ffsamp)
		               ? g_cyc_sign - g_cyc_gram - g_cyc_ffsamp : 0;
		s_other[completed] = other;
		completed++;
	}

	printf("\n======================================================================\n");
	printf("  SIGN  Falcon-%u   (%d valid sigs", 1u << logn, trials);
	if (attempts > trials) printf(", rejection rate %.1f%%", 100.0*(attempts-trials)/attempts);
	printf(")\n");
	printf("======================================================================\n");

	Stats st  = compute_stats(s_total,  trials);
	Stats sg  = compute_stats(s_gram,   trials);
	Stats sl  = compute_stats(s_ldl,    trials);
	Stats sbs = compute_stats(s_bs,     trials);
	Stats sbe = compute_stats(s_berexp, trials);
	Stats sf  = compute_stats(s_fftops, trials);
	Stats so  = compute_stats(s_other,  trials);
	uint64_t tot = st.med;

	print_header();
	print_row("Gram (B·B* in FFT domain)", sg, tot);
	print_row("LDL build (poly_LDL_fft)", sl, tot);
	print_row("BaseSampler (RCDT)",       sbs, tot);
	print_row("berexp+approxexp",         sbe, tot);
	print_row("FFT ops (split/merge/mul)",sf, tot);
	print_row("hash+target+compress",     so, tot);
	printf("  %-22s  ", "TOTAL");
	{ char mc[32], tc[32]; fmt_cyc(mc, tot); fmt_time(tc, tot);
	  printf("%s  %s\n", mc, tc); }

	/* Find dominant domain */
	const char *names[] = {"Gram","LDL build","BaseSampler","berexp","FFT ops","hash/enc"};
	uint64_t meds[] = {sg.med, sl.med, sbs.med, sbe.med, sf.med, so.med};
	int worst = 0;
	for (int i = 1; i < 6; i++) if (meds[i] > meds[worst]) worst = i;
	printf("\n  *** Bottleneck: [%s]  %.1f%% ***\n", names[worst], 100.0*meds[worst]/tot);
	printf("  BaseSampler calls/sign: %llu\n", g_n_samplerz);

	free(sig); free(tmp);
}

/* ── Verify benchmark ───────────────────────────────────────────────────────── */
static void bench_verify(unsigned logn, int trials,
                         const uint8_t *pk, size_t pklen,
                         const uint8_t **sigs, const size_t *siglens,
                         const uint8_t *msg, size_t mlen) {
	size_t tmplen = FALCON_TMPSIZE_VERIFY(logn);
	uint8_t *tmp  = malloc(tmplen + 63);
	uint8_t *tmpa = (uint8_t *)(((uintptr_t)tmp + 63) & ~(uintptr_t)63);

	uint64_t s_total[MAX_TRIALS], s_ntt[MAX_TRIALS],
	         s_norm[MAX_TRIALS], s_other[MAX_TRIALS];

	printf("\n======================================================================\n");
	printf("  VERIFY  Falcon-%u   (%d trials)\n", 1u << logn, trials);
	printf("======================================================================\n");

	for (int i = 0; i < trials; i++) {
		vf_cyc_total = vf_cyc_ntt = vf_cyc_norm = 0;
		uint64_t t0 = tsc_start();
		int ok = falcon_verify(sigs[i % trials], siglens[i % trials],
		                       FALCON_SIG_COMPRESSED,
		                       pk, pklen, msg, mlen, tmpa, tmplen);
		uint64_t wall = tsc_stop() - t0;
		if (ok != 0) { fprintf(stderr, "verify failed! rc=%d\n", ok); exit(1); }

		s_total[i] = wall;
		s_ntt[i]   = vf_cyc_ntt;
		s_norm[i]  = vf_cyc_norm;
		s_other[i] = wall > vf_cyc_ntt + vf_cyc_norm
		           ? wall - vf_cyc_ntt - vf_cyc_norm : 0;
	}

	Stats st = compute_stats(s_total, trials);
	Stats sn = compute_stats(s_ntt,   trials);
	Stats sr = compute_stats(s_norm,  trials);
	Stats so = compute_stats(s_other, trials);
	uint64_t tot = st.med;

	print_header();
	print_row("NTT poly_mul (s1·h mod q)", sn, tot);
	print_row("Norm check (is_short)",      sr, tot);
	print_row("decode+hash+other",          so, tot);
	printf("  %-22s  ", "TOTAL");
	{ char mc[32], tc[32]; fmt_cyc(mc, tot); fmt_time(tc, tot);
	  printf("%s  %s\n", mc, tc); }

	printf("\n  *** Bottleneck: [NTT poly_mul]  %.1f%% ***\n", 100.0*sn.med/tot);

	free(tmp);
}

/* ── Throughput summary ─────────────────────────────────────────────────────── */
static void print_throughput(unsigned logn, uint64_t kg, uint64_t sg, uint64_t vf) {
	printf("\n======================================================================\n");
	printf("  THROUGHPUT SUMMARY  Falcon-%u  @ %.3f GHz\n", 1u << logn, cpu_ghz);
	printf("======================================================================\n");
	printf("  %-10s  %12s  %10s  %12s\n", "Operation", "Cycles(med)", "Time", "ops/sec");
	printf("  %s\n", "--------------------------------------------------------------");
	struct { const char *n; uint64_t c; } ops[] = {
		{"keygen", kg}, {"sign", sg}, {"verify", vf}
	};
	for (int i = 0; i < 3; i++) {
		char mc[32], tc[32]; uint64_t c = ops[i].c;
		fmt_cyc(mc, c); fmt_time(tc, c);
		printf("  %-10s  %12s  %10s  %12.1f\n",
		       ops[i].n, mc, tc, cpu_ghz * 1e9 / c);
	}
	printf("\n");
}

/* ── Main ───────────────────────────────────────────────────────────────────── */
int main(int argc, char *argv[]) {
	unsigned logn   = (argc >= 2) ? (unsigned)atoi(argv[1]) : 9;
	int kg_trials   = (argc >= 3) ? atoi(argv[2]) : 5;
	int sg_trials   = (argc >= 4) ? atoi(argv[3]) : 30;

	if (logn < 1 || logn > 10) {
		fprintf(stderr, "logn must be 1..10\n"); return 1;
	}

	get_cpu_freq();

	printf("\nFalcon-%u domain-level RDTSC profiler\n", 1u << logn);
	printf("  CPU    : %.3f GHz\n", cpu_ghz);
	printf("  keygen : %d trials\n", kg_trials);
	printf("  sign   : %d valid signatures\n", sg_trials);
	printf("  verify : %d trials\n\n", sg_trials);

	size_t pklen  = FALCON_PUBKEY_SIZE(logn);
	size_t sklen  = FALCON_PRIVKEY_SIZE(logn);
	size_t tmpkg  = FALCON_TMPSIZE_KEYGEN(logn);

	uint8_t *pk  = malloc(pklen);
	uint8_t *sk  = malloc(sklen);
	uint8_t *tmp = malloc(tmpkg + 63);
	uint8_t *tmpa= (uint8_t *)(((uintptr_t)tmp + 63) & ~(uintptr_t)63);

	/* ── Keygen ── */
	bench_keygen(logn, kg_trials);

	/* ── Build one key for sign/verify ── */
	shake256_context rng;
	shake256_init_prng_from_system(&rng);
	printf("\nBuilding key for sign/verify benchmarks ...");
	fflush(stdout);
	if (falcon_keygen_make(&rng, logn, sk, sklen, pk, pklen, tmpa, tmpkg) != 0) {
		fprintf(stderr, "keygen failed\n"); return 1;
	}
	printf(" done.\n");

	static const uint8_t msg[32] = "bench-domains-test-message-xxxx";

	/* ── Sign ── */
	bench_sign(logn, sg_trials, sk, sklen, msg, sizeof msg);

	/* ── Collect signatures for verify ── */
	size_t sigmax = FALCON_SIG_COMPRESSED_MAXSIZE(logn);
	size_t tmpsig = FALCON_TMPSIZE_SIGNDYN(logn);
	uint8_t **sigs    = malloc(sg_trials * sizeof *sigs);
	size_t   *siglens = malloc(sg_trials * sizeof *siglens);
	uint8_t  *tmpsg   = malloc(tmpsig + 63);
	uint8_t  *tmpsa   = (uint8_t *)(((uintptr_t)tmpsg + 63) & ~(uintptr_t)63);
	for (int i = 0; i < sg_trials; i++) {
		sigs[i] = malloc(sigmax);
		siglens[i] = sigmax;
		shake256_init_prng_from_system(&rng);
		if (falcon_sign_dyn(&rng, sigs[i], &siglens[i], FALCON_SIG_COMPRESSED,
		                    sk, sklen, msg, sizeof msg, tmpsa, tmpsig) != 0) {
			fprintf(stderr, "sign failed\n"); return 1;
		}
	}

	/* ── Verify ── */
	bench_verify(logn, sg_trials, pk, pklen,
	             (const uint8_t **)sigs, siglens, msg, sizeof msg);

	/* ── Throughput summary ── */
	/* Re-run one keygen/sign/verify to get representative medians in one shot */
	{
		uint64_t kg_med = kg_cyc_total;   /* last keygen from bench_keygen */
		uint64_t sg_med = g_cyc_sign;     /* last sign from bench_sign */
		uint64_t vf_med = vf_cyc_total;   /* last verify from bench_verify */
		print_throughput(logn, kg_med, sg_med, vf_med);
	}

	for (int i = 0; i < sg_trials; i++) free(sigs[i]);
	free(sigs); free(siglens); free(tmpsg);
	free(pk); free(sk); free(tmp);
	return 0;
}
