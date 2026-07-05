#include <stdint.h>

/* LFENCE + RDTSC + LFENCE: serialized read of the timestamp counter.
   Prevents out-of-order execution from leaking instructions across
   the measurement boundary. */
uint64_t rdtsc_start() {
    uint32_t lo, hi;
    __asm__ __volatile__ (
        "lfence\n\t"
        "rdtsc\n\t"
        "lfence"
        : "=a"(lo), "=d"(hi)
    );
    return ((uint64_t)hi << 32) | lo;
}

/* RDTSCP: self-serializing on the way out — waits for all prior
   instructions to retire before reading the counter. */
uint64_t rdtsc_stop() {
    uint32_t lo, hi, aux;
    __asm__ __volatile__ (
        "rdtscp"
        : "=a"(lo), "=d"(hi), "=c"(aux)
    );
    return ((uint64_t)hi << 32) | lo;
}
