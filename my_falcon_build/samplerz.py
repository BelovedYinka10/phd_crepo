"""
samplerz.py — Discrete Gaussian sampler for Falcon
====================================================

Samples integers from the discrete Gaussian distribution D_{Z, σ, μ}:
    Pr[z = k]  ∝  exp( -(k - μ)² / (2σ²) )

Algorithm (Falcon spec §3.9, Algorithm 12)
------------------------------------------
  1. basesampler()  — samples z0 ≥ 0 from a half-Gaussian using RCDT
  2. approxexp()    — polynomial approximation of 2^63 · ccs · exp(-x)
  3. berexp()       — Bernoulli bit with probability ≈ ccs · exp(-x)
  4. samplerz()     — full sampler: combines above to produce z ~ D_{Z,σ,μ}

Security note
-------------
  Randomness comes from os.urandom (cryptographically secure).
  All arithmetic is branch-free where possible to resist timing attacks.

Reference
---------
  Falcon spec v1.2 §3.9 and reference Python implementation by Thomas Prest.
"""

__all__ = ["samplerz", "basesampler", "berexp", "approxexp"]

from math import floor
from os   import urandom


# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

# Upper bound on σ — the RCDT table is built for this value
MAX_SIGMA    = 1.8205
INV_2SIGMA2  = 1.0 / (2.0 * MAX_SIGMA ** 2)   # 1 / (2 · σ_max²)

# RCDT table precision: 72 bits → 9 bytes of randomness per base sample
RCDT_PREC = 72

# Natural log constants used in berexp
LN2  = 0.69314718056     # ln(2)
ILN2 = 1.44269504089     # 1 / ln(2)


# ─────────────────────────────────────────────────────────────────────────────
# RCDT — Reverse Cumulative Distribution Table
# ─────────────────────────────────────────────────────────────────────────────
#
# RCDT[k] = Pr[z0 > k] * 2^72   for z0 ~ half-Gaussian(σ_max)
#
# basesampler draws a 72-bit uniform u and counts how many RCDT entries
# exceed u — that count is z0.  This is a branch-free CDF inversion.
#
# Table has 18 entries → z0 ∈ {0, 1, …, 18}.
# (Tail probability beyond 18 is negligible for σ_max = 1.8205.)
#
RCDT = [
    3024686241123004913666,
    1564742784480091954050,
     636254429462080897535,
     199560484645026482916,
      47667343854657281903,
       8595902006365044063,
       1163297957344668388,
        117656387352093658,
          8867391802663976,
           496969357462633,
            20680885154299,
               638331848991,
                14602316184,
                   247426747,
                     3104126,
                       28824,
                         198,
                           1,
]


# ─────────────────────────────────────────────────────────────────────────────
# C — Polynomial coefficients for exp(-x) approximation
# ─────────────────────────────────────────────────────────────────────────────
#
# Approximates  2^63 · exp(-x)  for x ∈ [0, ln2) via a degree-12 polynomial:
#
#   approx = (2^-63) · Σ C[12-i] · x^i   for i = 0..12
#
# Lifted from FACCT: https://doi.org/10.1109/TC.2019.2940949
#
C = [
    0x00000004741183A3,
    0x00000036548CFC06,
    0x0000024FDCBF140A,
    0x0000171D939DE045,
    0x0000D00CF58F6F84,
    0x000680681CF796E3,
    0x002D82D8305B0FEA,
    0x011111110E066FD0,
    0x0555555555070F00,
    0x155555555581FF00,
    0x400000000002B400,
    0x7FFFFFFFFFFF4800,
    0x8000000000000000,
]


# ─────────────────────────────────────────────────────────────────────────────
# basesampler — sample z0 from half-Gaussian(σ_max)
# ─────────────────────────────────────────────────────────────────────────────

def basesampler(randombytes=urandom) -> int:
    """
    Sample z0 ∈ {0, 1, …, 18} from the half-Gaussian D_{Z≥0, σ_max}.

    Method: RCDT (Reverse Cumulative Distribution Table) inversion.
      1. Draw u uniformly at random from [0, 2^72)
      2. z0 = #{k : RCDT[k] > u}   (count how many table entries exceed u)

    This is branch-free: every entry is compared regardless of outcome,
    giving constant-time behaviour (no early exit).
    """
    u  = int.from_bytes(randombytes(RCDT_PREC >> 3), "little")   # 9 random bytes
    z0 = 0
    for elt in RCDT:
        z0 += int(u < elt)
    return z0


# ─────────────────────────────────────────────────────────────────────────────
# approxexp — integer approximation of 2^63 · ccs · exp(-x)
# ─────────────────────────────────────────────────────────────────────────────

def approxexp(x: float, ccs: float) -> int:
    """
    Compute an integer approximation of  2^63 · ccs · exp(-x).

    Both x and ccs must be positive.

    Steps:
      1. Evaluate the degree-12 polynomial in C at z = x · 2^63
         using Horner's method (iterative: y = C[i] - (z·y >> 63))
      2. Scale by ccs: y ← (ccs · 2^64) · y >> 63

    The result y satisfies:   y ≈ 2^63 · ccs · exp(-x)
    """
    y = C[0]
    z = int(x * (1 << 63))           # fixed-point x in [0, 2^63)
    for elt in C[1:]:
        y = elt - ((z * y) >> 63)    # Horner step
    z = int(ccs * (1 << 63)) << 1    # ccs · 2^64
    y = (z * y) >> 63
    return y


# ─────────────────────────────────────────────────────────────────────────────
# berexp — Bernoulli sampler: 1 with probability ≈ ccs · exp(-x)
# ─────────────────────────────────────────────────────────────────────────────

def berexp(x: float, ccs: float, randombytes=urandom) -> bool:
    """
    Return 1 with probability ≈ ccs · exp(-x),  0 otherwise.

    Both x and ccs must be positive.

    Method: bit-by-bit comparison
      Write x = s·ln2 + r  (s integer, r ∈ [0, ln2))
      Then  exp(-x) = 2^{-s} · exp(-r)

      approxexp(r, ccs) ≈ 2^63 · ccs · exp(-r)
      Right-shift by s gives  z ≈ 2^{63-s} · ccs · exp(-x)

      Compare z against a fresh random 64-bit number byte by byte.
      Return 1 iff  random < z.
    """
    # Use floor (not truncation) so that r = x - s*ln2 stays in [0, ln2)
    # even when x < 0.  Python's int() truncates towards zero, which gives
    # r < 0 for negative x — breaking approxexp's assumed input range.
    s = int(floor(x * ILN2))        # s = ⌊x / ln2⌋  (can be negative)
    r = x - s * LN2                 # r ∈ [0, ln2)  always

    # approxexp returns ~ 2^64 · ccs · exp(-r)
    # Right-shift by s divides by 2^s  →  z ≈ 2^{64-s} · ccs · exp(-r)
    #                                       = 2^64 · ccs · exp(-x)
    # For s < 0 (x < 0): left-shift instead; if z overflows 64 bits → accept.
    if s < 0:
        z = (approxexp(r, ccs) - 1) << (-s)
        if z.bit_length() > 64:
            return True             # probability ≥ 1 → always accept
    else:
        s = min(s, 63)              # clamp: exp(-x) < 2^{-63} → z ≈ 0 → reject
        z = (approxexp(r, ccs) - 1) >> s

    # Compare z against 8 random bytes (64 bits) byte by byte
    for i in range(56, -8, -8):
        p = int.from_bytes(randombytes(1), "little")
        w = p - ((z >> i) & 0xFF)
        if w:
            break
    return w < 0


# ─────────────────────────────────────────────────────────────────────────────
# samplerz — main discrete Gaussian sampler
# ─────────────────────────────────────────────────────────────────────────────

def samplerz(mu: float, sigma: float, sigmin: float,
             randombytes=urandom) -> int:
    """
    Sample z ~ D_{Z, σ, μ}  (discrete Gaussian, center μ, std dev σ).

    Inputs
    ------
    mu     : center of the distribution (floating-point)
    sigma  : standard deviation  (must satisfy sigmin < sigma < MAX_SIGMA)
    sigmin : lower bound on sigma  (from params: sigma_min for this n)

    Output
    ------
    Integer z with  Pr[z=k] ∝ exp(-(k-μ)²/(2σ²))

    Algorithm (Falcon spec Algorithm 12)
    -------------------------------------
    1. s ← ⌊μ⌋,   r ← μ - s          (split μ into integer + fractional)
    2. ccs ← σ_min / σ                (compression-correction scaling)
    3. Loop:
         z0 ← basesampler()           (sample non-negative half-Gaussian)
         b  ← random bit
         z  ← b + (2b-1)·z0          (apply random sign: z ∈ Z)
             if b=0: z = -z0  (≤ 0)
             if b=1: z = z0+1 (≥ 1)
         x  ← (z-r)²·dss - z0²·INV_2SIGMA2
         if berexp(x, ccs): return s + z

    The rejection step ensures the output is a true discrete Gaussian
    (not just the base distribution).
    """
    s   = int(floor(mu))             # integer part
    r   = mu - s                     # fractional part  r ∈ [0, 1)
    dss = 1.0 / (2.0 * sigma * sigma)   # 1 / (2σ²)
    ccs = sigmin / sigma             # compression-correction scaling

    while True:
        # Step 1: sample z0 from half-Gaussian(σ_max)
        z0 = basesampler(randombytes=randombytes)

        # Step 2: apply random sign to get z ∈ Z
        b  = int.from_bytes(randombytes(1), "little") & 1
        z  = b + (2 * b - 1) * z0   # b=0 → z=-z0;  b=1 → z=z0+1

        # Step 3: rejection test — accept with prob ∝ exp(-(z-r)²/(2σ²))
        x  = ((z - r) ** 2) * dss - (z0 ** 2) * INV_2SIGMA2
        if berexp(x, ccs, randombytes=randombytes):
            return z + s
