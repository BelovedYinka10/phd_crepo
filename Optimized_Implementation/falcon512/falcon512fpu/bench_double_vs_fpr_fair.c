#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
  #include <x86intrin.h>
#else
  #include <time.h>
#endif

// IMPORTANT: Use inner.h NOT fpr.h
#include "inner.h"

#ifndef ITER
#define ITER 10000000UL
#endif

#if defined(__GNUC__) || defined(__clang__)
  #define NOINLINE __attribute__((noinline))
#else
  #define NOINLINE __declspec(noinline)
#endif

static volatile double sink_double = 0.0;
static volatile uint64_t sink_u64 = 0;

/* Serialized RDTSC */
static inline uint64_t rdtsc_start(void) {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    _mm_lfence();
    return __rdtsc();
#else
    return (uint64_t)clock();
#endif
}

static inline uint64_t rdtsc_end(void) {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    unsigned aux;
    uint64_t t = __rdtscp(&aux);
    _mm_lfence();
    return t;
#else
    return (uint64_t)clock();
#endif
}

/* Prevent compiler optimization */
static NOINLINE fpr add_fpr_noinline(fpr a, fpr b) {
    return fpr_add(a, b);
}

int main(void) {
    printf("========================================\n");
    printf(" FAIR Benchmark: Hardware Double vs Falcon FPR\n");
    printf(" ITER = %lu\n", (unsigned long)ITER);
    printf("========================================\n\n");

    /* -------- Correctness Check -------- */
    double da = 2.0, db = 2.0;
    double dc = da + db;

    fpr fa = fpr_of(2);
    fpr fb = fpr_of(2);
    fpr fc = add_fpr_noinline(fa, fb);

    printf("Correctness Check:\n");
    printf("  double: 2 + 2 = %.2f\n", dc);
    printf("  fpr   : 2 + 2 = %lld (via fpr_rint)\n\n", (long long)fpr_rint(fc));

    /* -------- Hardware Double Benchmark -------- */
    volatile double x = 0.0;
    volatile double a = 1.0;
    volatile double b = 2.0;

    for (int i = 0; i < 100000; i++) {
        x = x + a + b;
    }

    uint64_t t0 = rdtsc_start();
    for (unsigned long i = 0; i < ITER; i++) {
        x = x + a + b;
    }
    uint64_t t1 = rdtsc_end();

    double cycles_double = (double)(t1 - t0) / (double)ITER;

    /* -------- Falcon FPR Benchmark -------- */
    fpr y = fpr_of(0);
    fpr one = fpr_of(1);
    fpr two = fpr_of(2);

    for (int i = 0; i < 100000; i++) {
        y = add_fpr_noinline(y, one);
        y = add_fpr_noinline(y, two);
    }

    uint64_t t2 = rdtsc_start();
    for (unsigned long i = 0; i < ITER; i++) {
        y = add_fpr_noinline(y, one);
        y = add_fpr_noinline(y, two);
    }
    uint64_t t3 = rdtsc_end();

    double cycles_fpr = (double)(t3 - t2) / (double)ITER;

    /* Prevent optimization removal */
    sink_double += x;
    sink_u64 ^= (uint64_t)fpr_rint(y);

    printf("Final Values:\n");
    printf("  double result: %.2f\n", x);
    printf("  fpr result   : %lld\n\n", (long long)fpr_rint(y));

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    printf("Average CPU cycles per operation:\n");
    printf("  Hardware double add : %.2f cycles\n", cycles_double);
    printf("  Falcon fpr_add      : %.2f cycles\n", cycles_fpr);
    printf("  Slowdown (fpr/double): %.2fx\n", cycles_fpr / cycles_double);
#else
    printf("Non-x86 platform: used clock() timing\n");
#endif

    printf("\n[Research Note]\n");
    printf("- inner.h resolves Zf() macro correctly\n");
    printf("- This is the proper way Falcon benchmarks internal fpr\n");

    return 0;
}
