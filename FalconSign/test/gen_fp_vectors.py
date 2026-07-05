"""
gen_fp_vectors.py — Generate test vectors for the 9 FP replacement modules.

For each operation, Python computes the IEEE 754 reference result using its
native double (which is identical to IEEE 754 binary64).  The RTL must produce
bit-exact results for the normal-number cases; special cases (NaN/Inf) must
follow the saturation / pass-through rules documented in each module.

Output: fp_vectors.txt — one vector per line:
  OP A_hex [B_hex] EXPECTED_hex

All hex fields are 16 digits (64 bits), zero-padded.
For int32 results/inputs the value is sign-extended to 64 bits.
"""

import struct, random, math, sys

NVEC = int(sys.argv[1]) if len(sys.argv) > 1 else 500
rng  = random.Random(42)

def to_bits(f):
    return struct.unpack('Q', struct.pack('d', f))[0]

def from_bits(b):
    return struct.unpack('d', struct.pack('Q', b))[0]

def rand_double():
    """Random normal double in (-1e15, 1e15)."""
    return rng.uniform(-1e15, 1e15)

def rand_normal_bits():
    return to_bits(rand_double())

def trunc_to_int(f):
    """Truncate toward zero, matching C (int) cast."""
    return int(math.trunc(f))

lines = []

# ── ADD ──────────────────────────────────────────────────────────────────────
for _ in range(NVEC):
    a = rand_double(); b = rand_double()
    r = a + b
    lines.append(f"ADD {to_bits(a):016X} {to_bits(b):016X} {to_bits(r):016X}")

# ── SUB ──────────────────────────────────────────────────────────────────────
for _ in range(NVEC):
    a = rand_double(); b = rand_double()
    r = a - b
    lines.append(f"SUB {to_bits(a):016X} {to_bits(b):016X} {to_bits(r):016X}")

# ── MUL ──────────────────────────────────────────────────────────────────────
for _ in range(NVEC):
    a = rand_double(); b = rand_double()
    r = a * b
    lines.append(f"MUL {to_bits(a):016X} {to_bits(b):016X} {to_bits(r):016X}")

# ── FLT2I32 (double → int32, truncate toward zero) ───────────────────────────
INT32_MAX =  (1 << 31) - 1
INT32_MIN = -(1 << 31)
for _ in range(NVEC):
    # mix of small and large values
    if rng.random() < 0.5:
        f = rng.uniform(-2**30, 2**30)  # fits in int32
    else:
        f = rng.uniform(-1e15, 1e15)    # may overflow
    i = trunc_to_int(f)
    i = max(INT32_MIN, min(INT32_MAX, i))   # saturate
    # sign-extend to 64 bits
    u64 = i & 0xFFFF_FFFF_FFFF_FFFF
    lines.append(f"FLT2I32 {to_bits(f):016X} {u64:016X}")

# ── I2FLT32 (int32 → double, exact) ──────────────────────────────────────────
for _ in range(NVEC):
    i = rng.randint(INT32_MIN, INT32_MAX)
    r = float(i)
    # input: sign-extended to 64 bits
    u64_in = i & 0xFFFF_FFFF_FFFF_FFFF
    lines.append(f"I2FLT32 {u64_in:016X} {to_bits(r):016X}")

# ── FLT2I64 (double → int64, truncate toward zero) ───────────────────────────
INT64_MAX =  (1 << 63) - 1
INT64_MIN = -(1 << 63)
for _ in range(NVEC):
    if rng.random() < 0.5:
        f = rng.uniform(-2**52, 2**52)   # fits exactly in double mantissa
    else:
        f = rng.uniform(-2**62, 2**62)   # large but fits in int64
    i = trunc_to_int(f)
    i = max(INT64_MIN, min(INT64_MAX, i))
    u64 = i & 0xFFFF_FFFF_FFFF_FFFF
    lines.append(f"FLT2I64 {to_bits(f):016X} {u64:016X}")

# ── I2FLT64 (int64 → double, may round for |i| > 2^52) ──────────────────────
for _ in range(NVEC):
    i = rng.randint(INT64_MIN, INT64_MAX)
    r = float(i)
    u64_in = i & 0xFFFF_FFFF_FFFF_FFFF
    lines.append(f"I2FLT64 {u64_in:016X} {to_bits(r):016X}")

# ── Special cases: zero, inf, NaN (shared across add/sub/mul) ────────────────
ZERO_POS = 0x0000_0000_0000_0000
ZERO_NEG = 0x8000_0000_0000_0000
INF_POS  = 0x7FF0_0000_0000_0000
INF_NEG  = 0xFFF0_0000_0000_0000

for a_bits, b_bits in [
    (to_bits(1.5), ZERO_POS),
    (ZERO_POS, to_bits(2.5)),
    (INF_POS, to_bits(1.0)),
    (to_bits(1.0), INF_POS),
]:
    a = from_bits(a_bits); b = from_bits(b_bits)
    lines.append(f"ADD {a_bits:016X} {b_bits:016X} {to_bits(a+b):016X}")
    lines.append(f"SUB {a_bits:016X} {b_bits:016X} {to_bits(a-b):016X}")
    lines.append(f"MUL {a_bits:016X} {b_bits:016X} {to_bits(a*b):016X}")

out = "fp_vectors.txt"
with open(out, "w") as f:
    f.write("\n".join(lines) + "\n")

counts = {}
for l in lines:
    op = l.split()[0]
    counts[op] = counts.get(op, 0) + 1

print(f"Written {len(lines)} vectors to {out}")
for op, n in sorted(counts.items()):
    print(f"  {op:10s}: {n}")
