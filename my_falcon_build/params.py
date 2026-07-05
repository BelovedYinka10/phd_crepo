"""
params.py — Falcon parameter sets for n = 512 and n = 1024.

Reference: Falcon specification v1.2, Table 3.3 and §3.4, §3.9.3, §3.11.
"""

import math


# ── Shared constants ───────────────────────────────────────────────────────────

Q        = 12289   # Ring modulus: prime, = 1 + 12·1024  (§3.4.1, eq 2.10)
SALT_LEN = 40      # Salt r length in bytes (320 bits)    (§3.9.1)
PHI      = "x^n+1" # Cyclotomic polynomial (symbolic note)


# ── Production parameter sets (Table 3.3) ─────────────────────────────────────
#
# sigma      : std deviation used during signing (Fast Fourier Sampling)
# sigma_min  : minimum leaf std dev in the Falcon tree (used in SamplerZ)
# sigma_max  : maximum leaf std dev in the Falcon tree (used in SamplerZ)
# beta_sq    : max squared norm of a valid signature ⌊β²⌋
# pk_bytelen : public key size in bytes (including 1-byte header)
# sig_bytelen: signature size in bytes  (padded, including header + salt)

FALCON_PARAMS = {
    512: {
        "n"          : 512,
        "sigma"      : 165.736617183,   # §Table 3.3
        "sigma_min"  : 1.277833697,     # §Table 3.3 / §3.9.3
        "sigma_max"  : 1.8205,          # §Table 3.3
        "beta_sq"    : 34_034_726,      # §Table 3.3
        "pk_bytelen" : 897,             # §Table 3.3
        "sig_bytelen": 666,             # §Table 3.3
    },
    1024: {
        "n"          : 1024,
        "sigma"      : 168.388571447,   # §Table 3.3
        "sigma_min"  : 1.298280334,     # §Table 3.3 / §3.9.3
        "sigma_max"  : 1.8205,          # §Table 3.3
        "beta_sq"    : 70_265_242,      # §Table 3.3
        "pk_bytelen" : 1793,            # §Table 3.3
        "sig_bytelen": 1280,            # §Table 3.3
    },
}

# Valid degrees (powers of 2 from 2 to 1024)
# Small values (2..256) are used for testing only — no security
VALID_DEGREES = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]


# ── Derived: key-generation Gaussian std deviation ────────────────────────────
#
# From Algorithm 5 (NTRUGen):
#   sigma_fg = 1.17 * sqrt(q) / sqrt(2n)
# Chosen so that E[||(f,g)||] = 1.17 * sqrt(q),
# which minimises the Gram-Schmidt norm ||B||_GS.  (§2.3.3, eq 2.12)

def sigma_fg(n):
    """
    Std deviation for sampling each coefficient of f and g during key generation.
    Both f and g have n coefficients, each drawn independently from D_{Z, sigma_fg}.
    """
    return 1.17 * math.sqrt(Q) / math.sqrt(2 * n)


# ── Derived: signing sigma for arbitrary n ────────────────────────────────────
#
# From §2.6 (eq 2.13):
#   sigma = (1/pi) * sqrt(ln(4n(1 + 1/eps)) / 2) * 1.17 * sqrt(q)
# For production (n=512, n=1024) use the exact values from Table 3.3 above.
# For testing at smaller n this formula gives approximate values.

def signing_sigma(n, eps=1e-9):
    """Approximate signing std deviation for degree n (use Table 3.3 for n=512/1024)."""
    return (1.0 / math.pi) * math.sqrt(math.log(4 * n * (1 + 1 / eps)) / 2) * 1.17 * math.sqrt(Q)


# ── Private key coefficient bit widths (§3.11.5) ──────────────────────────────
#
# f and g coefficients (signed, two's complement, minimal value forbidden):
#   n <= 32   → 8 bits  (range −127 to +127, −128 forbidden)
#   n  = 64   → 7 bits
#   n  = 128  → 7 bits
#   n  = 256  → 6 bits
#   n  = 512  → 6 bits  (range −31 to +31, −32 forbidden)
#   n  = 1024 → 5 bits  (range −15 to +15, −16 forbidden)
#
# F coefficients always use 8 bits regardless of degree.
# G is NOT encoded — recomputed from f, g, F via the NTRU equation (§3.11.5).

FG_BITWIDTH = {
    2: 8,
    4: 8,
    8: 8,
    16: 8,
    32: 8,
    64: 7,
    128: 7,
    256: 6,
    512: 6,
    1024: 5,
}
F_BITWIDTH = 8   # F always uses 8 bits


def fg_max_coeff(n):
    """Maximum absolute value of f/g coefficients for degree n."""
    bits = FG_BITWIDTH[n]
    return (1 << (bits - 1)) - 1   # e.g. bits=6 → max=31


# ── Header bytes for encoding (§3.11.3, §3.11.4, §3.11.5) ────────────────────
#
# Format: high nibble encodes type, low nibble encodes log2(n)
#
#   Public key : 0000 nnnn  → base = 0x00
#   Private key: 0101 nnnn  → base = 0x50
#   Signature  : 0011 nnnn  → base = 0x30  (cc=01 = compressed)
#
# Example: Falcon-512 public key header  = 0x00 | 9  = 0x09
#          Falcon-512 private key header = 0x50 | 9  = 0x59
#          Falcon-512 signature header   = 0x30 | 9  = 0x39

HEAD_PK  = 0x00   # Public key
HEAD_SK  = 0x50   # Private key
HEAD_SIG = 0x30   # Signature (compressed encoding, cc=01)


def log2n(n):
    """Return log2(n); n must be a power of two."""
    assert n in VALID_DEGREES, f"n={n} is not a supported degree"
    return n.bit_length() - 1


def pk_header(n):
    return HEAD_PK | log2n(n)

def sk_header(n):
    return HEAD_SK | log2n(n)

def sig_header(n):
    return HEAD_SIG | log2n(n)


# ── Signature encoding constants (§3.11.2) ────────────────────────────────────
#
# Each signature coefficient s_i is encoded as:
#   [sign 1 bit] [|s_i| low 7 bits in binary] [|s_i| >> 7 in unary: 0^k 1]
#
# The 7 low bits are sent as-is (close to uniform → not worth compressing).
# The high bits use unary (= Huffman code for the Gaussian high-bit distribution).

SIG_LOWBITS  = 7          # Number of low bits encoded in binary
SIG_LOWMASK  = (1 << SIG_LOWBITS) - 1   # 0x7F


# ── Sanity check ──────────────────────────────────────────────────────────────

def get_params(n):
    """Return the parameter dict for n=512 or n=1024."""
    assert n in FALCON_PARAMS, f"n={n} not a production parameter set (use 512 or 1024)"
    return FALCON_PARAMS[n]


if __name__ == "__main__":
    print("=" * 60)
    print("Falcon Parameter Sets")
    print("=" * 60)

    print(f"\nShared:")
    print(f"  Q (modulus)  = {Q}")
    print(f"  Salt length  = {SALT_LEN} bytes ({SALT_LEN * 8} bits)")

    for n in [512, 1024]:
        p = FALCON_PARAMS[n]
        print(f"\nFalcon-{n}:")
        print(f"  sigma        = {p['sigma']}")
        print(f"  sigma_min    = {p['sigma_min']}")
        print(f"  sigma_max    = {p['sigma_max']}")
        print(f"  beta_sq      = {p['beta_sq']:,}")
        print(f"  pk_bytelen   = {p['pk_bytelen']} bytes")
        print(f"  sig_bytelen  = {p['sig_bytelen']} bytes")
        print(f"  sigma_fg     = {sigma_fg(n):.10f}  (keygen Gaussian)")
        print(f"  f/g bitwidth = {FG_BITWIDTH[n]} bits  (max coeff ±{fg_max_coeff(n)})")
        print(f"  pk header    = 0x{pk_header(n):02X}")
        print(f"  sk header    = 0x{sk_header(n):02X}")
        print(f"  sig header   = 0x{sig_header(n):02X}")
