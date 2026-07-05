#!/usr/bin/env python3
"""
gen_vectors.py — emit (u, z0) test vectors for the Verilog BaseSampler.

Golden model is the Python reference in ../samplerz.py.  Each vector is
cross-checked two ways before being written:

  1. direct RCDT count:   z0 = sum(int(u < elt) for elt in RCDT)
  2. reference call:       basesampler() fed the exact 9 bytes of u

so the testbench is checking the RTL against the same algorithm the rest of
the Falcon build uses, not a re-derivation.

Output: vectors.txt — one "<u_hex(18 digits)> <z0_dec>" pair per line.

Usage:  python3 gen_vectors.py [num_random_vectors]   (default 2000)
"""
import os
import sys
import random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from samplerz import RCDT, basesampler   # noqa: E402

U_BITS = 72
U_MAX  = (1 << U_BITS) - 1


def z0_from_u(u):
    """The BaseSampler core, in one line — the spec algorithm."""
    return sum(int(u < elt) for elt in RCDT)


def z0_reference(u):
    """Drive the actual library basesampler() with the bytes of u."""
    b = u.to_bytes(9, "little")            # 72 bits, little-endian (matches lib)
    return basesampler(randombytes=lambda k, _b=b: _b[:k])


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    out = os.path.join(os.path.dirname(__file__), "vectors.txt")
    rng = random.Random(0xFA1C0)

    vecs = []

    # Boundary coverage: straddle every RCDT threshold + the two extremes.
    for elt in RCDT:
        for d in (-1, 0, 1):
            v = elt + d
            if 0 <= v <= U_MAX:
                vecs.append(v)
    vecs.append(0)          # below every entry  -> z0 = 18
    vecs.append(U_MAX)      # above every entry  -> z0 = 0

    # Random coverage.
    for _ in range(n):
        vecs.append(rng.getrandbits(U_BITS))

    with open(out, "w") as f:
        for u in vecs:
            z = z0_from_u(u)
            assert z == z0_reference(u), f"ref mismatch at u={u}: {z} vs ref"
            f.write(f"{u:018x} {z}\n")

    print(f"wrote {len(vecs)} vectors -> {out}")
    print(f"  ({3 * len(RCDT)} boundary + 2 extreme + {n} random)")


if __name__ == "__main__":
    main()
