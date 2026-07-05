"""
sign.py — Falcon Sign and Verify
==================================

Implements sign() and verify() for Falcon-512 and Falcon-1024.

Sign pipeline
-------------
  1. salt  ← random(40 bytes)
  2. c     = hash_to_point(salt || message)      SHAKE256 → Z_q polynomial
  3. t     = (c, 0) · B^{-1} / q               FFT-domain target vector
  4. z     = ffsampling_fft(t, T, sigmin, rng)  Gaussian lattice sample
  5. v     = z · B                               nearest lattice point
  6. (s0, s1) = (c − v0, −v1)                   short signature vector
  7. check ||s||² ≤ β²; encode s1               retry if norm or length fails

Verify pipeline
---------------
  1. c  = hash_to_point(salt || message)
  2. s0 = c − s1·h  mod q  (centered in (−q/2, q/2])
  3. check ||s0||² + ||s1||² ≤ β²

Reference
---------
  Falcon spec v1.2, Algorithms 10–11
  Reference Python by Thomas Prest (PQShield) — falcon.py-master 3/falcon.py
"""

__all__ = ["sign", "verify", "hash_to_point", "compress", "decompress"]

import hashlib
import os

from params import Q, SALT_LEN, FALCON_PARAMS
from fft import fft, ifft, poly_mul_fft, poly_add_fft
from ffsampling import ffsampling_fft
from keygen import expand_secret_key
from poly_arith import poly_mul, poly_center


# ─────────────────────────────────────────────────────────────────────────────
# Hash to point  (Falcon spec §3.7.1)
# ─────────────────────────────────────────────────────────────────────────────

def hash_to_point(message: bytes, salt: bytes, n: int) -> list:
    """
    SHAKE256(salt || message) → n coefficients in Z_q.

    Rejection sampling: read 16-bit words from the XOF; accept elt if
    elt < k·q (k = 5), then take elt mod q.  This gives a uniform distribution
    over Z_q with ~6.25% rejection rate.
    """
    k = (1 << 16) // Q          # 65536 // 12289 = 5
    nbytes = n * 3               # yields ~1.5n candidates; far more than enough
    data = hashlib.shake_256(salt + message).digest(nbytes)
    c = []
    idx = 0
    while len(c) < n:
        if idx + 2 > len(data):  # expand if (very unlikely) we run dry
            nbytes *= 2
            data = hashlib.shake_256(salt + message).digest(nbytes)
        elt = (data[idx] << 8) | data[idx + 1]
        idx += 2
        if elt < k * Q:
            c.append(elt % Q)
    return c


# ─────────────────────────────────────────────────────────────────────────────
# Signature encoding  (Falcon spec §3.11.2)
# ─────────────────────────────────────────────────────────────────────────────

def compress(v: list, slen: int):
    """
    Encode list of signed integers v into exactly slen bytes.

    Per coefficient: [sign 1 bit] [|coef| low 7 bits, MSB first] [|coef|>>7 zeros + 1]
    Returns bytes on success, False if slen is too small.
    """
    bits = []
    for coef in v:
        bits.append(1 if coef < 0 else 0)
        low = abs(coef) & 0x7F
        for i in range(6, -1, -1):
            bits.append((low >> i) & 1)
        high = abs(coef) >> 7
        bits.extend([0] * high)
        bits.append(1)                   # unary terminator

    if len(bits) > 8 * slen:
        return False

    bits.extend([0] * (8 * slen - len(bits)))
    out = bytearray(slen)
    for i in range(slen):
        b = 0
        for j in range(8):
            b = (b << 1) | bits[8 * i + j]
        out[i] = b
    return bytes(out)


def decompress(data: bytes, slen: int, n: int):
    """
    Decode slen bytes into n signed integers.
    Returns list[int] or False if the encoding is invalid.
    """
    if len(data) > slen:
        return False

    bits = []
    for byte in data:
        for i in range(7, -1, -1):
            bits.append((byte >> i) & 1)

    while bits and bits[-1] == 0:       # strip padding zeros
        bits.pop()

    v = []
    idx = 0
    try:
        while idx < len(bits) and len(v) < n:
            sign = -1 if bits[idx] else 1
            idx += 1
            low = 0
            for _ in range(7):
                low = (low << 1) | bits[idx]
                idx += 1
            high = 0
            while bits[idx] == 0:
                high += 1
                idx += 1
            idx += 1                    # consume unary terminator
            coef = sign * (low + (high << 7))
            if coef == 0 and sign == -1:
                return False            # −0 is not a valid encoding
            v.append(coef)
    except IndexError:
        return False

    return v if len(v) == n else False


# ─────────────────────────────────────────────────────────────────────────────
# Sign  (Falcon spec Algorithm 10)
# ─────────────────────────────────────────────────────────────────────────────

def sign(sk, message: bytes, rng=os.urandom) -> bytes:
    """
    Sign a message using a Falcon secret key.

    Args:
        sk:      SecretKey from keygen() — tree or dynamic mode both work
        message: arbitrary bytes
        rng:     callable(nbytes) → bytes; default os.urandom

    Returns:
        signature: header(1) + salt(40) + compressed_s1
    """
    n   = sk.n
    p   = FALCON_PARAMS[n]
    sigmin      = p["sigma_min"]
    beta_sq     = p["beta_sq"]
    payload_len = p["sig_bytelen"] - 1 - SALT_LEN

    # Expand dynamic-mode key to get B_fft and T
    if sk.mode == "dynamic":
        sk = expand_secret_key(sk)

    [[a, b], [c_mat, d]] = sk.B_fft
    # a = fft(g), b = fft(-f), c_mat = fft(G), d = fft(-F)

    header = bytes([0x30 | (n.bit_length() - 1)])

    while True:
        salt   = rng(SALT_LEN)
        c_poly = hash_to_point(message, salt, n)
        c_fft  = fft([float(x) for x in c_poly])

        # Target t = (c, 0) · B^{-1} in FFT domain.
        # B^{-1} = (1/q)[[-F, f], [-G, g]],  so:
        #   t0 = c·(−F)/q  = c_fft * d / q
        #   t1 = c·f/q     = −c_fft * b / q  (b = fft(−f))
        t0 = [x / Q for x in poly_mul_fft(c_fft, d)]
        t1 = [-x / Q for x in poly_mul_fft(c_fft, b)]

        # Sample z ∼ D_{Z^{2n}, σ} centered at t
        z = ffsampling_fft([t0, t1], sk.T, sigmin, rng)

        # Nearest lattice point v = z · B
        v0_fft = poly_add_fft(poly_mul_fft(z[0], a), poly_mul_fft(z[1], c_mat))
        v1_fft = poly_add_fft(poly_mul_fft(z[0], b), poly_mul_fft(z[1], d))
        v0 = [int(round(x)) for x in ifft(v0_fft)]
        v1 = [int(round(x)) for x in ifft(v1_fft)]

        # Signature vector (s0, s1) = (c − v0, −v1)
        s0 = [c_poly[i] - v0[i] for i in range(n)]
        s1 = [-v1[i] for i in range(n)]

        # Reject if norm too large
        norm_sq = sum(x * x for x in s0) + sum(x * x for x in s1)
        if norm_sq > beta_sq:
            continue

        # Reject if compressed s1 doesn't fit
        enc = compress(s1, payload_len)
        if enc is False:
            continue

        return header + salt + enc


# ─────────────────────────────────────────────────────────────────────────────
# Verify  (Falcon spec Algorithm 12 / §3.10)
# ─────────────────────────────────────────────────────────────────────────────

def verify(pk, message: bytes, signature: bytes) -> bool:
    """
    Verify a Falcon signature.

    Args:
        pk:        PublicKey from keygen()
        message:   the signed message (bytes)
        signature: bytes produced by sign()

    Returns:
        True if the signature is valid, False otherwise.
    """
    n = pk.n
    h = pk.h
    p = FALCON_PARAMS[n]
    beta_sq     = p["beta_sq"]
    payload_len = p["sig_bytelen"] - 1 - SALT_LEN

    if len(signature) < 1 + SALT_LEN:
        return False

    salt   = signature[1 : 1 + SALT_LEN]
    enc_s1 = signature[1 + SALT_LEN :]

    # Decode s1
    s1 = decompress(enc_s1, payload_len, n)
    if s1 is False:
        return False

    # Hash message to the same point
    c = hash_to_point(message, salt, n)

    # Recover s0 = c − s1·h  mod q  (centered)
    s1_mod = [x % Q for x in s1]
    s1h    = poly_mul(s1_mod, h)
    s0     = [(c[i] - s1h[i]) % Q for i in range(n)]
    s0     = poly_center(s0)

    # Accept iff ||s0||² + ||s1||² ≤ β²
    norm_sq = sum(x * x for x in s0) + sum(x * x for x in s1)
    return norm_sq <= beta_sq
