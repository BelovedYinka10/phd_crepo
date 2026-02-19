#include <stdio.h>
#include <stdint.h>
#include <x86intrin.h>
#include "inner.h"

#define ITER 10000000

/* Prevent inlining (CRITICAL) */
__attribute__((noinline)) fpr bench_fpr_add(fpr a, fpr b) {
    return fpr_add(a, b);
}

__attribute__((noinline)) double bench_double_add(double a, double b) {
    return a + b;
}

/* Serialized RDTSC */
static inline uint64_t rdtsc_start() {
    _mm_lfence();
    return __rdtsc();
}

static inline uint64_t rdtsc_end() {
    uint64_t t = __rdtsc();
    _mm_lfence();
    return t;
}

int main(void)
{
    /* Use volatile + runtime values to defeat constant folding */
    volatile double a = 2.0;
    volatile double b = 2.0;
    volatile double c = 0.0;

    volatile fpr fa = fpr_of(2);
    volatile fpr fb = fpr_of(2);
    volatile fpr fc = fpr_of(0);

    /* ================= HARDWARE DOUBLE ================= */
    uint64_t start_hw = rdtsc_start();
    for (int i = 0; i < ITER; i++) {
        c = bench_double_add(c, a + b);
    }
    uint64_t end_hw = rdtsc_end();

    /* ================= FALCON FPU BACKEND ================= */
    uint64_t start_fpu = rdtsc_start();
    for (int i = 0; i < ITER; i++) {
        fc = bench_fpr_add(fc, bench_fpr_add(fa, fb));
    }
    uint64_t end_fpu = rdtsc_end();

    int64_t decoded = fpr_rint(fc);

    double cycles_hw = (double)(end_hw - start_hw) / ITER;
    double cycles_fpu = (double)(end_fpu - start_fpu) / ITER;

    printf("========================================\n");
    printf(" TRUE Benchmark (No Optimization Collapse)\n");
    printf(" ITER = %d\n", ITER);
    printf("========================================\n\n");

    printf("Final Values:\n");
    printf("  double result: %.2f\n", c);
    printf("  fpr result   : %lld\n\n", decoded);

    printf("Average CPU cycles per operation:\n");
    printf("  Hardware double add : %.2f cycles\n", cycles_hw);
    printf("  Falcon FPU fpr_add  : %.2f cycles\n", cycles_fpu);
    printf("  Slowdown (FPU/double): %.2fx\n", cycles_fpu / cycles_hw);

    printf("\n[Valid Benchmark] noinline + volatile + dependency chain\n");

    return 0;
}
