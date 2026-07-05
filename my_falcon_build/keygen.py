"""
keygen.py — Falcon Key Generation
===================================

Generates Falcon key pairs: secret key (for signing) and public key
(for verification).  Supports two modes:

  "tree"    — precomputes B_fft and LDL tree T during keygen.
              Signing is fast (tree already built).
              Secret key is larger.

  "dynamic" — stores only (f, g, F, G).
              Tree is rebuilt each time you sign.
              Secret key is smaller.

Pipeline (tree mode)
--------------------
  1. ntru_gen(n) → (f, g, F, G)        NTRU key quadruple
  2. B = [[g, -f], [G, -F]]            NTRU lattice basis
  3. gram(B) → G_mat → fft → ffldl → normalize → T
  4. fft(B) → B_fft
  5. h = g · f⁻¹ mod q                 Public key

Pipeline (dynamic mode)
-----------------------
  1. ntru_gen(n) → (f, g, F, G)
  2. h = g · f⁻¹ mod q

Reference
---------
  Falcon spec v1.2, Algorithm 5 (NTRUGen) + key formatting
  Reference Python implementation by Thomas Prest (PQShield)
"""

__all__ = [
    "keygen", "SecretKey", "PublicKey",
    "expand_secret_key",
    "serialize_public_key", "deserialize_public_key",
]

from collections import namedtuple

from params import Q, FALCON_PARAMS
from ntt import ntt, intt
from fft import fft
from ntrugen import ntru_gen
from ffsampling import gram, ffldl_fft, normalize_tree


# ──────────────────────────────────────────────────────────────────────────────
# Key structures
# ──────────────────────────────────────────────────────────────────────────────

SecretKey = namedtuple("SecretKey", ["n", "mode", "f", "g", "F", "G", "B_fft", "T"])
PublicKey = namedtuple("PublicKey", ["n", "h"])


# ──────────────────────────────────────────────────────────────────────────────
# NTT-domain polynomial division mod q
# ──────────────────────────────────────────────────────────────────────────────

def poly_div_ntt(g, f):
    """
    Compute h = g · f⁻¹ mod (q, x^n+1) using NTT.

    Steps:
      1. f_ntt = NTT(f mod q),  g_ntt = NTT(g mod q)
      2. h_ntt[i] = g_ntt[i] · f_ntt[i]⁻¹ mod q   (Fermat inverse)
      3. h = iNTT(h_ntt)

    Raises ZeroDivisionError if f is not invertible mod q.
    """
    n = len(f)
    f_mod = [c % Q for c in f]
    g_mod = [c % Q for c in g]
    f_ntt = ntt(f_mod)
    g_ntt = ntt(g_mod)

    if any(v == 0 for v in f_ntt):
        raise ZeroDivisionError("f is not invertible mod q in NTT domain")

    h_ntt = [g_ntt[i] * pow(f_ntt[i], Q - 2, Q) % Q for i in range(n)]
    return intt(h_ntt)


# ──────────────────────────────────────────────────────────────────────────────
# Key generation
# ──────────────────────────────────────────────────────────────────────────────

def keygen(n, mode="tree"):
    """
    Generate a Falcon key pair for degree n.

    Args:
        n: polynomial degree, must be 512 or 1024
        mode: "tree"    — precompute B_fft and LDL tree (fast signing)
              "dynamic" — store only (f,g,F,G) (small key, slower signing)

    Returns:
        (sk, pk) where:
          sk = SecretKey(n, mode, f, g, F, G, B_fft, T)
          pk = PublicKey(n, h)

    In dynamic mode, sk.B_fft and sk.T are None.
    Call expand_secret_key(sk) to compute them before signing.
    """
    assert n in FALCON_PARAMS, f"Unsupported degree n={n}"
    assert mode in ("tree", "dynamic"), f"Unknown mode: {mode}"

    # Step 1: Generate NTRU key quadruple
    f, g, F, G = ntru_gen(n)

    # Step 2: Public key h = g / f mod q
    h = poly_div_ntt(g, f)

    # Step 3: Precompute tree data (tree mode only)
    if mode == "tree":
        B_fft, T = _build_tree(n, f, g, F, G)
    else:
        B_fft, T = None, None

    sk = SecretKey(n=n, mode=mode, f=f, g=g, F=F, G=G, B_fft=B_fft, T=T)
    pk = PublicKey(n=n, h=h)
    return sk, pk


def _build_tree(n, f, g, F, G):
    """
    Build the FFT basis and normalized LDL tree from (f, g, F, G).

    Returns (B_fft, T).
    """
    sigma = FALCON_PARAMS[n]["sigma"]

    neg_f = [-c for c in f]
    neg_F = [-c for c in F]
    B = [[g, neg_f], [G, neg_F]]

    # Gram matrix → FFT domain → LDL tree → normalize
    G_mat = gram(B)
    G_fft = [[fft(G_mat[i][j]) for j in range(2)] for i in range(2)]
    T = ffldl_fft(G_fft)
    normalize_tree(T, sigma)

    # FFT-domain basis
    B_fft = [[fft([float(c) for c in poly]) for poly in row] for row in B]

    return B_fft, T


def expand_secret_key(sk):
    """
    Expand a dynamic-mode secret key into tree mode.

    Takes a SecretKey with B_fft=None, T=None and returns a new
    SecretKey with B_fft and T computed.

    Useful for sign_dyn: rebuild the tree before each signature.
    """
    B_fft, T = _build_tree(sk.n, sk.f, sk.g, sk.F, sk.G)
    return SecretKey(n=sk.n, mode="tree", f=sk.f, g=sk.g,
                     F=sk.F, G=sk.G, B_fft=B_fft, T=T)


# ──────────────────────────────────────────────────────────────────────────────
# Public key serialization (14 bits per coefficient)
# ──────────────────────────────────────────────────────────────────────────────

LOGN = {512: 9, 1024: 10}

def serialize_public_key(h, n):
    """
    Serialize public key h to bytes.

    Format: 1 header byte (0x00 | logn) + 14 bits per coefficient.
    Coefficients must be in [0, q).

    Returns bytes of length 1 + ceil(14*n / 8).
    """
    logn = LOGN[n]
    header = bytes([0x00 | logn])

    # Pack 14-bit values into a byte stream
    acc = 0
    acc_len = 0
    buf = bytearray()
    for coeff in h:
        acc = (acc << 14) | (coeff % Q)
        acc_len += 14
        while acc_len >= 8:
            acc_len -= 8
            buf.append((acc >> acc_len) & 0xFF)

    # Flush remaining bits (padded with zeros)
    if acc_len > 0:
        buf.append((acc << (8 - acc_len)) & 0xFF)

    return header + bytes(buf)


def deserialize_public_key(data, n):
    """
    Deserialize public key from bytes.

    Returns h as a list of n integers in [0, q).
    """
    # Skip header byte
    payload = data[1:]

    acc = 0
    acc_len = 0
    h = []
    idx = 0
    for _ in range(n):
        while acc_len < 14:
            acc = (acc << 8) | payload[idx]
            idx += 1
            acc_len += 8
        acc_len -= 14
        val = (acc >> acc_len) & 0x3FFF  # 14-bit mask
        if val >= Q:
            raise ValueError(f"Coefficient {val} >= q={Q}")
        h.append(val)

    return h
