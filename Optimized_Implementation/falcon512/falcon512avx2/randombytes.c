#include <windows.h>
#include <bcrypt.h>
#include <stdint.h>

#pragma comment(lib, "bcrypt.lib")

void randombytes(uint8_t *out, size_t outlen) {
    BCryptGenRandom(NULL, out, outlen, BCRYPT_USE_SYSTEM_PREFERRED_RNG);
}
