#include <stdio.h>
#include <stdint.h>

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
  #include <x86intrin.h>   // __rdtsc()
#else
  #include <time.h>
#endif

// Use Falcon's real internal API (emulated fpr backend in reference dir)
#include "inner.h"

#ifndef ITER
#define ITER 10000000UL
#endif

// --- prevent inlining + prevent optimizer removing the work ---
#if defined(__GNUC__) || defined(__clang__)
  #define NOINLINE __attribute__((noinline))
#else
  #define NOINLINE __declspec(noinline)
#endif

static volatile uint64_t sink_u64;
static volatile int64_t  sink_i64;

// A tiny noinline wrapper so the compiler can't fold everything away
static NOINLINE fpr add_fpr_noinline(fpr a, fpr b) {
    return fpr_add(a, b);
}

int main(void) {
    printf("========================================\n");
    printf(" Emulated Falcon fpr vs Hardware double\n");
    printf(" ITER = %lu\n", (unsigned long)ITER);
    printf("========================================\n\n");

    // --------- correctness demo: 2 + 2 ----------
    double da = 2.0, db = 2.0;
    double dc = da + db;

    fpr fa = fpr_of(2);
    fpr fb = fpr_of(2);
    fpr fc = add_fpr_noinline(fa, fb);

    // decode fpr -> int (for easy viewing)
    int64_t fc_int = fpr_rint(fc);

    printf("Values:\n");
    printf("  double: 2.0 + 2.0 = %.2f\n", dc);
    printf("  fpr   : 2 + 2     = %lld (decoded via fpr_rint)\n", (long long)fc_int);
    printf("  fpr raw bits      = %llu\n\n", (unsigned long long)fc);

    // --------- benchmark hardware double add ----------
    volatile double vda = 2.0;
    volatile double vdb = 2.0;
    volatile double vdc = 0.0;

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    unsigned long long t0 = __rdtsc();
    for (unsigned long i = 0; i < ITER; i++) {
        // dependency chain: prevents full unrolling + constant folding
        vdc = vdc + vda + vdb;
    }
    unsigned long long t1 = __rdtsc();
    double cycles_double = (double)(t1 - t0) / (double)ITER;
#else
    clock_t t0 = clock();
    for (unsigned long i = 0; i < ITER; i++) {
        vdc = vdc + vda + vdb;
    }
    clock_t t1 = clock();
    double cycles_double = 0.0; // no rdtsc on this platform
#endif

    // --------- benchmark Falcon emulated fpr_add ----------
    fpr x = fpr_of(0);
    fpr one = fpr_of(1);
    fpr two = fpr_of(2);

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    unsigned long long t2 = __rdtsc();
    for (unsigned long i = 0; i < ITER; i++) {
        // dependency chain using noinline function
        x = add_fpr_noinline(x, one);
        x = add_fpr_noinline(x, two);
    }
    unsigned long long t3 = __rdtsc();
    double cycles_fpr = (double)(t3 - t2) / (double)ITER;
#else
    clock_t t2 = clock();
    for (unsigned long i = 0; i < ITER; i++) {
        x = add_fpr_noinline(x, one);
        x = add_fpr_noinline(x, two);
    }
    clock_t t3 = clock();
    double cycles_fpr = 0.0;
#endif

    // sink results (so compiler can't throw them away)
    sink_u64 ^= (uint64_t)x;
    sink_i64 ^= fpr_rint(x);

    printf("Final Values:\n");
    printf("  double result: %.2f\n", (double)vdc);
    printf("  fpr result   : %lld\n\n", (long long)fpr_rint(x));

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    printf("Average CPU cycles per operation:\n");
    printf("  Hardware double add : %.2f cycles\n", cycles_double);
    printf("  Falcon fpr_add      : %.2f cycles\n", cycles_fpr);
    printf("  Slowdown (fpr/double): %.2fx\n", cycles_fpr / cycles_double);
#else
    printf("Timing measured with clock() (no rdtsc available)\n");
#endif

    return 0;
}
