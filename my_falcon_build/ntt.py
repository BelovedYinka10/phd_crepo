"""
ntt.py — Iterative NTT / iNTT for Falcon
=========================================

Negacyclic Number Theoretic Transform over Z_q[x]/(x^n+1)
  q = 12289,  n ∈ {512, 1024}

Forward NTT  : iterative Cooley-Tukey (CT) butterfly
Inverse iNTT : iterative Gentleman-Sande (GS) butterfly

Two arithmetic backends are provided so you can benchmark them:
  1. Plain '%' arithmetic  — clean, readable Python
  2. Montgomery arithmetic — shift-based, no modular division, FPGA-ready

Spec reference : Falcon spec §3.8
C reference    : falcon-round3/.../vrfy.c  (mq_NTT / mq_iNTT)
"""

# ──────────────────────────────────────────────────────────────────────────────
# Exports
# ──────────────────────────────────────────────────────────────────────────────
__all__ = [
    "ntt", "intt",
    "ntt_montgomery",
    "intt_montgomery",
    "ntt_mul",
    "ntt_add",
    "ntt_sub",
    "poly_mul_ntt",
    "benchmark",
]

# ──────────────────────────────────────────────────────────────────────────────
# Domain constants
# ──────────────────────────────────────────────────────────────────────────────
from params import Q  # single source of truth: Q = 12289
G    = 11             # primitive root mod Q  (ord(7)=2048 ≠ Q-1; ord(11)=12288 ✓)
LOGN = {512: 9, 1024: 10}   # log2(n) for each supported n

# ──────────────────────────────────────────────────────────────────────────────
# Montgomery constants   (R = 2^16)
#
#  R      = 65536          "Montgomery radix"
#  Q0I    = 12287          satisfies  Q · Q0I ≡ −1  (mod R)
#  R_MOD_Q  = R mod Q = 4091     = multiplicative identity in Montgomery domain
#  R2_MOD_Q = R² mod Q = 10952   used to convert a value INTO Montgomery domain
#
# Verification:
#   Q · Q0I = 12289 · 12287 = 150 994 943 = 2305 · 65536 − 1  ✓
#   R mod Q = 65536 mod 12289 = 4091                           ✓
#   R² mod Q = 4091² mod 12289 = 16 736 281 mod 12289 = 10952 ✓
# ──────────────────────────────────────────────────────────────────────────────
_R        = 1 << 16          # 65536
_Q0I      = 12287
_R_MOD_Q  = _R % Q           # 4091
_R2_MOD_Q = (_R_MOD_Q * _R_MOD_Q) % Q   # 10952


# ──────────────────────────────────────────────────────────────────────────────
# Bit-reversal helper
# ──────────────────────────────────────────────────────────────────────────────
def _bitrev(k: int, logn: int) -> int:
    """Reverse the bottom `logn` bits of k."""
    r = 0
    for _ in range(logn):
        r = (r << 1) | (k & 1)
        k >>= 1
    return r


# ──────────────────────────────────────────────────────────────────────────────
# Montgomery arithmetic
# ──────────────────────────────────────────────────────────────────────────────
def mq_montymul(x: int, y: int) -> int:
    """
    Montgomery multiplication: returns x·y·R⁻¹ mod Q.

    Step-by-step:
      z = x·y                          (raw product)
      t = (z · Q0I) mod R              (find how many Q's to subtract)
      z = (z + t·Q) >> log2(R)         (subtract and shift — exact division by R)
      conditional subtract if z ≥ Q
    """
    z = x * y
    w = ((z * _Q0I) & 0xFFFF) * Q
    z = (z + w) >> 16
    return z - Q if z >= Q else z


def mq_add(x: int, y: int) -> int:
    """Branchless modular addition: (x + y) mod Q, no Python %."""
    z = x + y - Q
    return z + Q if z < 0 else z


def mq_sub(x: int, y: int) -> int:
    """Branchless modular subtraction: (x − y) mod Q, no Python %."""
    z = x - y
    return z + Q if z < 0 else z


def mq_tomonty(x: int) -> int:
    """Convert x  →  x·R mod Q  (into Montgomery domain)."""
    return mq_montymul(x, _R2_MOD_Q)


def mq_frommonty(x: int) -> int:
    """Convert x·R mod Q  →  x mod Q  (out of Montgomery domain)."""
    return mq_montymul(x, 1)


# ──────────────────────────────────────────────────────────────────────────────
# Twiddle-factor tables
# ──────────────────────────────────────────────────────────────────────────────
#
# For a negacyclic NTT of size n we need ψ, a primitive 2n-th root of unity:
#
#   ψ = g^((q−1)/(2n))  mod q,   g = 11 (true primitive root; ord(7)=2048≠q−1)
#
# The CT and GS butterflies at stage m use the twiddle factor stored at
# position (m + i) in GMB / IGMB:
#
#   GMB[k]  = ψ^bitrev(k, logn)          — forward twiddles
#   IGMB[k] = ψ_inv^bitrev(k, logn)      — inverse twiddles
#
# Why bit-reversed order?  The iterative CT loop increments m as
# 1, 2, 4, …, n/2.  At each level, the required twiddle power doubles,
# which corresponds exactly to a bit-reversal of the index.
#
# Montgomery-domain variants (used by the Montgomery backend):
#   GMB_M[k]  = GMB[k]  · R mod Q
#   IGMB_M[k] = IGMB[k] · R mod Q
# ──────────────────────────────────────────────────────────────────────────────
def _build_tables(n: int):
    logn   = LOGN[n]
    psi     = pow(G, (Q - 1) // (2 * n), Q)   # primitive 2n-th root of unity
    psi_inv = pow(psi, Q - 2, Q)               # ψ⁻¹ mod Q  (Fermat's little theorem)

    # plain powers:  powers[k] = ψ^k
    powers_fwd = [pow(psi,     k, Q) for k in range(n)]
    powers_inv = [pow(psi_inv, k, Q) for k in range(n)]

    # bit-reversed reordering into the twiddle tables
    gmb  = [powers_fwd[_bitrev(k, logn)] for k in range(n)]
    igmb = [powers_inv[_bitrev(k, logn)] for k in range(n)]

    # Montgomery-domain versions
    gmb_m  = [mq_tomonty(v) for v in gmb]
    igmb_m = [mq_tomonty(v) for v in igmb]

    return gmb, igmb, gmb_m, igmb_m


# ── precompute at import time ──────────────────────────────────────────────────
_TABLES: dict = {}
for _n in (512, 1024):
    _gmb, _igmb, _gmb_m, _igmb_m = _build_tables(_n)
    _TABLES[_n] = {
        "gmb":   _gmb,    # plain forward twiddles
        "igmb":  _igmb,   # plain inverse twiddles
        "gmb_m": _gmb_m,  # Montgomery forward twiddles
        "igmb_m":_igmb_m, # Montgomery inverse twiddles
    }

# n⁻¹ mod Q for final scaling in iNTT
N_INV = {n: pow(n, Q - 2, Q) for n in (512, 1024)}


# ──────────────────────────────────────────────────────────────────────────────
# Iterative CT NTT  (plain arithmetic)
# ──────────────────────────────────────────────────────────────────────────────
def ntt(f: list) -> list:
    """
    Forward NTT using iterative Cooley-Tukey butterflies.

    Ring: Z_q[x]/(x^n+1),  q=12289,  n∈{512,1024}

    Butterfly at each stage (one CT step):
        u = a[j]
        v = a[j+t] · s          (s = twiddle factor)
        a[j]   = u + v  mod q
        a[j+t] = u − v  mod q

    Stage loop: m = 1, 2, 4, …, n/2   (half-width t = n, n/2, …, 1)
    Twiddle:    s = GMB[m + i]  for group i in 0..m−1

    Input : list of n integers in Z_q  (natural order)
    Output: list of n integers in Z_q  (NTT domain, natural order)
    """
    n    = len(f)
    gmb  = _TABLES[n]["gmb"]
    a    = list(f)
    t    = n      # current butterfly half-width
    m    = 1      # number of groups at this stage

    while m < n:
        t >>= 1
        for i in range(m):
            j1 = 2 * i * t          # start of this group's first butterfly
            s  = gmb[m + i]         # twiddle factor for this group
            for j in range(j1, j1 + t):
                u        = a[j]
                v        = a[j + t] * s % Q
                a[j]     = (u + v) % Q
                a[j + t] = (u - v) % Q
        m <<= 1

    return a


# ──────────────────────────────────────────────────────────────────────────────
# Iterative GS iNTT  (plain arithmetic)
# ──────────────────────────────────────────────────────────────────────────────
def intt(f_ntt: list) -> list:
    """
    Inverse NTT using iterative Gentleman-Sande butterflies.

    Butterfly at each stage (one GS step):
        u = a[j]
        v = a[j+t]
        a[j]   = u + v       mod q
        a[j+t] = (u − v) · s  mod q   (s = inverse twiddle)

    Stage loop: m = n, n/2, …, 2   (i.e., h = m/2 = n/2, n/4, …, 1)
    Twiddle:    s = IGMB[h + i]  for group i in 0..h−1

    After all stages, multiply every coefficient by n⁻¹ mod q.

    Input : list of n integers in NTT domain
    Output: list of n integers in Z_q (natural order)
    """
    n    = len(f_ntt)
    igmb = _TABLES[n]["igmb"]
    a    = list(f_ntt)
    t    = 1
    m    = n

    while m > 1:
        h  = m >> 1       # half the group count — also the twiddle base index
        j1 = 0
        for i in range(h):
            s = igmb[h + i]
            for j in range(j1, j1 + t):
                u        = a[j]
                v        = a[j + t]
                a[j]     = (u + v) % Q
                a[j + t] = (u - v) * s % Q
            j1 += 2 * t
        t <<= 1
        m >>= 1

    # final scaling by n⁻¹
    n_inv = N_INV[n]
    return [x * n_inv % Q for x in a]


# ──────────────────────────────────────────────────────────────────────────────
# Iterative CT NTT  (Montgomery arithmetic)
# ──────────────────────────────────────────────────────────────────────────────
def ntt_montgomery(f: list) -> list:
    """
    Forward NTT using Montgomery arithmetic throughout.

    Montgomery butterfly:
        v = mq_montymul(a[j+t], s_m)   where s_m = twiddle · R mod Q
        a[j]   = mq_add(u, v)
        a[j+t] = mq_sub(u, v)

    All values live in Montgomery domain during the transform;
    inputs are converted in, outputs are converted out.

    This is the FPGA-friendly version: every modular multiply
    becomes a multiply + shift (no division).

    Input : list of n integers in Z_q  (natural order)
    Output: list of n integers in Z_q  (NTT domain, natural order)
    """
    n     = len(f)
    gmb_m = _TABLES[n]["gmb_m"]

    # bring input into Montgomery domain: a[j] = f[j] · R mod Q
    a = [mq_tomonty(x) for x in f]
    t = n
    m = 1

    while m < n:
        t >>= 1
        for i in range(m):
            j1 = 2 * i * t
            s  = gmb_m[m + i]        # twiddle already in Montgomery domain
            for j in range(j1, j1 + t):
                u        = a[j]
                v        = mq_montymul(a[j + t], s)
                a[j]     = mq_add(u, v)
                a[j + t] = mq_sub(u, v)
        m <<= 1

    # convert back from Montgomery domain
    return [mq_frommonty(x) for x in a]


# ──────────────────────────────────────────────────────────────────────────────
# Iterative GS iNTT  (Montgomery arithmetic)
# ──────────────────────────────────────────────────────────────────────────────
def intt_montgomery(f_ntt: list) -> list:
    """
    Inverse NTT using Montgomery arithmetic throughout.

    Montgomery GS butterfly:
        a[j]   = mq_add(u, v)
        a[j+t] = mq_montymul(mq_sub(u, v), s_m)

    Input : list of n integers in NTT domain
    Output: list of n integers in Z_q (natural order)
    """
    n      = len(f_ntt)
    igmb_m = _TABLES[n]["igmb_m"]
    # n_inv in Montgomery domain: mq_montymul(x, n_inv_m) scales x by n⁻¹
    n_inv_m = mq_tomonty(N_INV[n])

    # bring input into Montgomery domain
    a = [mq_tomonty(x) for x in f_ntt]
    t = 1
    m = n

    while m > 1:
        h  = m >> 1
        j1 = 0
        for i in range(h):
            s = igmb_m[h + i]
            for j in range(j1, j1 + t):
                u        = a[j]
                v        = a[j + t]
                a[j]     = mq_add(u, v)
                a[j + t] = mq_montymul(mq_sub(u, v), s)
            j1 += 2 * t
        t <<= 1
        m >>= 1

    # scale by n⁻¹, then exit Montgomery domain
    return [mq_frommonty(mq_montymul(x, n_inv_m)) for x in a]


# ──────────────────────────────────────────────────────────────────────────────
# Pointwise operations (NTT domain)
# ──────────────────────────────────────────────────────────────────────────────
def ntt_mul(a_ntt: list, b_ntt: list) -> list:
    """Pointwise multiplication in NTT domain."""
    return [x * y % Q for x, y in zip(a_ntt, b_ntt)]


def ntt_add(a_ntt: list, b_ntt: list) -> list:
    """Pointwise addition in NTT domain."""
    return [(x + y) % Q for x, y in zip(a_ntt, b_ntt)]


def ntt_sub(a_ntt: list, b_ntt: list) -> list:
    """Pointwise subtraction in NTT domain."""
    return [(x - y) % Q for x, y in zip(a_ntt, b_ntt)]


def poly_mul_ntt(f: list, g: list) -> list:
    """
    Multiply f · g  mod (x^n + 1, q)  using NTT.

    This is the fundamental operation for Falcon: used in key generation
    (NTRU equation), verification, and the polynomial arithmetic inside
    the Gram-Schmidt / LDL decomposition.

    Cost: 2 × NTT  +  n multiplications  +  1 × iNTT
    """
    return intt(ntt_mul(ntt(f), ntt(g)))


# ──────────────────────────────────────────────────────────────────────────────
# Benchmark harness
# ──────────────────────────────────────────────────────────────────────────────
def benchmark(n: int = 512, iters: int = 200) -> None:
    """
    Time plain vs Montgomery NTT / iNTT for the given n.
    Run as:  python ntt.py  (calls benchmark for both n=512 and n=1024)
    """
    import time, random
    f = [random.randrange(Q) for _ in range(n)]

    results = []
    for label, fn in [
        ("ntt (plain)",       ntt),
        ("intt (plain)",      intt),
        ("ntt (montgomery)",  ntt_montgomery),
        ("intt (montgomery)", intt_montgomery),
    ]:
        t0 = time.perf_counter()
        for _ in range(iters):
            fn(f)
        ms = (time.perf_counter() - t0) / iters * 1000
        results.append((label, ms))

    print(f"\n── n={n},  {iters} iterations ──────────────────────────")
    for label, ms in results:
        print(f"  {label:<26}  {ms:7.4f} ms/call")


# ──────────────────────────────────────────────────────────────────────────────
# Self-test
# ──────────────────────────────────────────────────────────────────────────────
def _test(n: int = 512) -> None:
    import random

    def poly_mul_naive(f, g, n, q):
        """Schoolbook negacyclic multiplication — reference for testing."""
        h = [0] * n
        for i, fi in enumerate(f):
            for j, gj in enumerate(g):
                idx = (i + j) % n
                sign = -1 if (i + j) >= n else 1
                h[idx] = (h[idx] + sign * fi * gj) % q
        return h

    random.seed(42)
    f = [random.randrange(Q) for _ in range(n)]
    g = [random.randrange(Q) for _ in range(n)]

    # ── round-trip tests ────────────────────────────────────────────────────
    assert intt(ntt(f)) == f, \
        f"[FAIL] n={n}: plain round-trip broken"
    assert intt_montgomery(ntt_montgomery(f)) == f, \
        f"[FAIL] n={n}: montgomery round-trip broken"

    # ── cross-backend consistency ────────────────────────────────────────────
    assert ntt(f) == ntt_montgomery(f), \
        f"[FAIL] n={n}: plain NTT ≠ montgomery NTT"

    # ── convolution correctness ──────────────────────────────────────────────
    expected = poly_mul_naive(f, g, n, Q)
    got      = poly_mul_ntt(f, g)
    assert got == expected, \
        f"[FAIL] n={n}: NTT convolution ≠ naive schoolbook"

    print(f"[OK]  n={n:4d} — round-trip, montgomery consistency, convolution all pass.")


if __name__ == "__main__":
    _test(512)
    _test(1024)
    benchmark(512)
    benchmark(1024)
