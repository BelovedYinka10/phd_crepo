#!/usr/bin/env python3
"""
gen_shake_vectors.py — emit SHAKE256 test vectors for keccak_shake256_wrapper.sv.

Golden model: hashlib.shake_256 (Python stdlib, NIST FIPS 202 compliant).
Message lengths are chosen to straddle the 136-byte (1088-bit) SHAKE256 rate
boundary: sub-block, exact block, block+1, exact double block, etc. Output
lengths are always a multiple of 4 (wrapper requirement) and also straddle
the rate boundary so the "multiple squeeze / re-permute" path is exercised.

Output: shake_vectors.txt
  line 1:            <num_vectors>
  per vector, 4 lines:
    <msg_len_bytes> <out_len_bytes>
    <msg bytes as space-separated hex, "" if msg_len==0>
    <expected output bytes as space-separated hex>
    (blank line separator)
"""
import hashlib
import os
import random

MSG_LENS = [1, 4, 32, 135, 136, 137, 200, 271, 272, 273, 300, 400]
OUT_LENS = [4, 32, 64, 136, 200, 400]


def main():
    rng = random.Random(0xC5411A)
    out_path = os.path.join(os.path.dirname(__file__), "shake_vectors.txt")

    vectors = []
    for mlen in MSG_LENS:
        for olen in OUT_LENS:
            msg = bytes(rng.getrandbits(8) for _ in range(mlen))
            exp = hashlib.shake_256(msg).digest(olen)
            vectors.append((msg, exp))

    # extra fully-random coverage
    for _ in range(30):
        mlen = rng.randint(1, 500)
        olen = rng.randint(1, 100) * 4
        msg = bytes(rng.getrandbits(8) for _ in range(mlen))
        exp = hashlib.shake_256(msg).digest(olen)
        vectors.append((msg, exp))

    with open(out_path, "w") as f:
        f.write(f"{len(vectors)}\n")
        for msg, exp in vectors:
            f.write(f"{len(msg)} {len(exp)}\n")
            f.write(" ".join(f"{b:02x}" for b in msg) + "\n")
            f.write(" ".join(f"{b:02x}" for b in exp) + "\n")

    print(f"wrote {len(vectors)} vectors -> {out_path}")


if __name__ == "__main__":
    main()
