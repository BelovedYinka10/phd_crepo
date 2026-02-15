#include <stdint.h>

#if defined(__APPLE__) || defined(__FreeBSD__)
  #include <sys/random.h>
#elif defined(__linux__)
  #include <unistd.h>
#else
  #error "randombytes(): unsupported platform"
#endif

/*
 * NIST-required RNG hook.
 * Falcon calls randombytes() internally.
 */
void randombytes(unsigned char *out, unsigned long long outlen) {
    if (getentropy(out, outlen) != 0) {
        __builtin_trap();   // fail hard if entropy is unavailable
    }
}
