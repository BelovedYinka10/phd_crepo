"""
test_sign.py — Tests for Falcon sign and verify
================================================

Tests:
  1. hash_to_point: output length, range, determinism
  2. compress / decompress: round-trip, rejection cases
  3. sign / verify: valid signature accepted
  4. verify: rejects tampered message
  5. verify: rejects tampered signature
  6. sign: different salts produce different signatures
  7. sign + verify: dynamic-mode key
  8. sign + verify: n=1024  (slow — skipped by default)
"""

import os
import time

from sign import sign, verify, hash_to_point, compress, decompress
from keygen import keygen
from params import Q, SALT_LEN, FALCON_PARAMS

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"

results = []

def check(name, cond):
    tag = PASS if cond else FAIL
    print(f"  [{tag}] {name}")
    results.append(cond)


# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Test 1: hash_to_point ===")

n = 16
salt = os.urandom(SALT_LEN)
msg  = b"hello falcon"
c    = hash_to_point(msg, salt, n)

check("output length == n",       len(c) == n)
check("all coefficients in [0,q)", all(0 <= x < Q for x in c))

c2 = hash_to_point(msg, salt, n)
check("deterministic (same salt → same c)", c == c2)

c3 = hash_to_point(msg, os.urandom(SALT_LEN), n)
check("different salt → different c (overwhelmingly)", c != c3)


# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Test 2: compress / decompress ===")

v    = [3, -7, 0, 127, -128, 1, -1, 256, -256]
slen = 40
enc  = compress(v, slen)
check("compress returns bytes",          isinstance(enc, bytes))
check("compressed length == slen",       len(enc) == slen)

v2 = decompress(enc, slen, len(v))
check("round-trip v → compress → decompress == v", v2 == v)

check("compress rejects slen=1 (too small)", compress(v, 1) is False)

# Typical Falcon-512 payload: n=512 coefficients, ~614 bytes
import random
random.seed(42)
v512 = [random.randint(-30, 30) for _ in range(512)]
p512 = FALCON_PARAMS[512]
slen512 = p512["sig_bytelen"] - 1 - SALT_LEN    # 625 bytes
enc512 = compress(v512, slen512)
v512b = decompress(enc512, slen512, 512)
check("compress/decompress 512 Falcon-range coefficients", v512 == v512b)


# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Test 3: sign / verify (n=512, tree mode) ===")
print("  [keygen n=512 …]", end=" ", flush=True)
t0 = time.perf_counter()
sk, pk = keygen(512, mode="tree")
print(f"done ({time.perf_counter()-t0:.1f}s)")

msg = b"Hello, Falcon-512!"
t0  = time.perf_counter()
sig = sign(sk, msg)
print(f"  [sign done in {(time.perf_counter()-t0)*1000:.1f} ms]")

check("signature is bytes",                      isinstance(sig, bytes))
check("signature length == sig_bytelen",         len(sig) == FALCON_PARAMS[512]["sig_bytelen"])
check("header byte correct (0x39)",              sig[0] == 0x39)
check("verify accepts valid signature",          verify(pk, msg, sig))


# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Test 4: verify rejects tampered message ===")

check("tampered message rejected",
      not verify(pk, b"tampered message", sig))


# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Test 5: verify rejects tampered signature ===")

sig_bad = bytearray(sig)
sig_bad[41] ^= 0xFF                 # flip bits in the salt
check("tampered salt rejected",     not verify(pk, msg, bytes(sig_bad)))

sig_bad2 = bytearray(sig)
sig_bad2[-1] ^= 0x01               # flip last byte of encoded s1
check("tampered s1 rejected",       not verify(pk, msg, bytes(sig_bad2)))


# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Test 6: different salts per signature ===")

sigs = [sign(sk, msg) for _ in range(5)]
salts = [s[1:41] for s in sigs]
check("5 signatures have 5 distinct salts", len(set(salts)) == 5)
check("all 5 signatures verify",            all(verify(pk, msg, s) for s in sigs))


# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Test 7: sign / verify (n=512, dynamic mode) ===")

sk_dyn, pk_dyn = keygen(512, mode="dynamic")
msg_dyn = b"dynamic mode test"
sig_dyn = sign(sk_dyn, msg_dyn)
check("dynamic mode: verify accepts signature", verify(pk_dyn, msg_dyn, sig_dyn))


# ─────────────────────────────────────────────────────────────────────────────
print("\n=== Summary ===")
total  = len(results)
passed = sum(results)
print(f"  {passed}/{total} tests passed\n")
if passed < total:
    raise SystemExit(1)
