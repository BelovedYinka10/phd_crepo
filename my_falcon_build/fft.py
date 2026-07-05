"""
fft.py — Iterative FFT / iFFT for Falcon (hardware-translatable)
=================================================================

Complex FFT over R[x]/(x^n+1) using floating-point arithmetic.

This is the *signing-side* transform (not the NTT used in verification).
Fast Fourier Sampling, LDL tree construction, and Gram-Schmidt all
operate in this FFT domain.

Design
------
  • Iterative Cooley-Tukey (forward) and Gentleman-Sande (inverse)
    — NO recursion, directly translatable to Verilog FSM
  • Split-halves memory layout (mirrors the C implementation):
      f[0 .. n/2-1]   = real parts
      f[n/2 .. n-1]   = imaginary parts
    One array of n floats stores n/2 complex values.
  • Precomputed twiddle table GM_TAB (bit-reversed roots of unity)

Reference: falcon-round3/Extra/c/fft.c  (Thomas Pornin / NCC Group)

Notation
--------
  n   = polynomial degree (512 or 1024)
  hn  = n // 2
  w   = exp(i·π/n)  — primitive 2n-th root of unity
  GM_TAB[k] = (Re(w^bitrev(k)), Im(w^bitrev(k)))   for k = 0..n-1
"""

__all__ = [
    # Core transforms
    "fft", "ifft",
    # FFT-domain polynomial operations
    "poly_add_fft", "poly_sub_fft", "poly_neg_fft",
    "poly_mul_fft", "poly_adj_fft", "poly_div_fft",
    "poly_mulconst_fft",
    "poly_mulselfadj_fft", "poly_muladj_fft",
    "poly_mul_autoadj_fft", "poly_div_autoadj_fft",
    "poly_invnorm2_fft",
    # Split / merge (needed by ffsampling)
    "poly_split_fft", "poly_merge_fft",
    # LDL operations
    "poly_LDL_fft", "poly_LDLmv_fft",
]

import math

# ──────────────────────────────────────────────────────────────────────────────
# Supported degrees
# ──────────────────────────────────────────────────────────────────────────────
LOGN = {2: 1, 4: 2, 8: 3, 16: 4, 32: 5, 64: 6, 128: 7, 256: 8, 512: 9, 1024: 10}


# ──────────────────────────────────────────────────────────────────────────────
# Bit-reversal helper
# ──────────────────────────────────────────────────────────────────────────────
def _bitrev(k: int, bits: int) -> int:
    """Reverse the bottom `bits` bits of k."""
    r = 0
    for _ in range(bits):
        r = (r << 1) | (k & 1)
        k >>= 1
    return r


# ──────────────────────────────────────────────────────────────────────────────
# Twiddle-factor table  (precomputed at import time)
# ──────────────────────────────────────────────────────────────────────────────
#
# GM_TAB[n] is a flat list of 2*n floats:
#   GM_TAB[n][2*k + 0] = Re( w^bitrev(k, logn) )
#   GM_TAB[n][2*k + 1] = Im( w^bitrev(k, logn) )
#
# where w = exp(i·π/n), a primitive 2n-th root of unity.
#
# The forward FFT indexes this as GM_TAB[n][2*(m + i) + 0/1].
# The inverse FFT uses the conjugate: same real part, negated imag part.
# ──────────────────────────────────────────────────────────────────────────────
GM_TAB: dict = {}

for _n, _logn in LOGN.items():
    _w = math.pi / _n                      # argument: π/n
    _tab = [0.0] * (2 * _n)
    for _k in range(_n):
        _rev = _bitrev(_k, _logn)
        _angle = _rev * _w                  # bitrev(k) · π/n
        _tab[2 * _k + 0] = math.cos(_angle)
        _tab[2 * _k + 1] = math.sin(_angle)
    GM_TAB[_n] = _tab

# Scaling constants:  P2_TAB[logn] = 2^(1-logn) = 2/n
# Used by iFFT for final normalisation (compensates for skipped last level).
P2_TAB = {logn: 2.0 / (1 << logn) for logn in LOGN.values()}


# ──────────────────────────────────────────────────────────────────────────────
# Forward FFT  (iterative Cooley-Tukey)
# ──────────────────────────────────────────────────────────────────────────────
def fft(f: list) -> list:
    """
    In-place iterative forward FFT.

    Input : n real coefficients  (polynomial in time domain)
    Output: n floats in split-halves layout:
              out[0..n/2-1]   = real parts of n/2 complex FFT values
              out[n/2..n-1]   = imaginary parts

    The first butterfly level (m=1, twiddle=i) is implicit in the
    split-halves representation, so the loop starts at m=2.

    Butterfly (CT):
        y = s · (y_re, y_im)       complex multiply by twiddle
        f[j]    = x + y            (top)
        f[j+ht] = x - y            (bottom)
    """
    n    = len(f)
    logn = LOGN[n]
    hn   = n >> 1
    tab  = GM_TAB[n]
    a    = list(f)          # work on a copy

    t = hn                  # half-width starts at n/2
    m = 2                   # skip m=1 (implicit)

    for _ in range(1, logn):        # logn-1 levels
        ht = t >> 1
        hm = m >> 1                 # number of butterfly groups

        for i1 in range(hm):
            j1   = i1 * t
            # twiddle factor:  s = GM[m + i1]
            idx  = (m + i1) << 1
            s_re = tab[idx]
            s_im = tab[idx + 1]

            for j in range(j1, j1 + ht):
                x_re = a[j]
                x_im = a[j + hn]
                y_re = a[j + ht]
                y_im = a[j + ht + hn]

                # complex multiply:  y = s * y
                ty_re = y_re * s_re - y_im * s_im
                ty_im = y_re * s_im + y_im * s_re

                # butterfly
                a[j]          = x_re + ty_re
                a[j + hn]     = x_im + ty_im
                a[j + ht]     = x_re - ty_re
                a[j + ht + hn] = x_im - ty_im

        t = ht
        m <<= 1

    return a


# ──────────────────────────────────────────────────────────────────────────────
# Inverse FFT  (iterative Gentleman-Sande)
# ──────────────────────────────────────────────────────────────────────────────
def ifft(f_fft: list) -> list:
    """
    In-place iterative inverse FFT.

    Input : n floats in split-halves FFT layout
    Output: n real coefficients (imaginary parts should be ~zero)

    The last butterfly level (m=2→1, twiddle=i) is implicit, so
    the loop stops at m=2 and the final scaling uses 2/n instead of 1/n.

    Butterfly (GS):
        f[j]   = x + y              (sum)
        f[j+t] = s · (x - y)        (twiddle the difference)

    where s = conj(GM[hm + i1])  (conjugate of forward twiddle).
    """
    n    = len(f_fft)
    logn = LOGN[n]
    hn   = n >> 1
    tab  = GM_TAB[n]
    a    = list(f_fft)

    t = 1
    m = n

    for _ in range(logn, 1, -1):    # logn-1 levels  (stop before last)
        hm = m >> 1
        dt = t << 1
        j1 = 0

        for i1 in range(hm):
            if j1 >= hn:
                break
            j2 = j1 + t

            # inverse twiddle = conjugate of forward: negate imag part
            idx  = (hm + i1) << 1
            s_re = tab[idx]
            s_im = -tab[idx + 1]       # conjugate

            for j in range(j1, j2):
                x_re = a[j]
                x_im = a[j + hn]
                y_re = a[j + t]
                y_im = a[j + t + hn]

                # sum → top
                a[j]      = x_re + y_re
                a[j + hn] = x_im + y_im

                # difference · twiddle → bottom
                d_re = x_re - y_re
                d_im = x_im - y_im
                a[j + t]      = d_re * s_re - d_im * s_im
                a[j + t + hn] = d_re * s_im + d_im * s_re

            j1 += dt

        t = dt
        m = hm

    # Final scaling: 2/n  (compensates for skipped last level + 1/n norm)
    ni = P2_TAB[logn]
    return [x * ni for x in a]


# ──────────────────────────────────────────────────────────────────────────────
# FFT-domain pointwise operations
# ──────────────────────────────────────────────────────────────────────────────
# All operate on the split-halves layout:
#   a[0..hn-1] = real,  a[hn..n-1] = imag
# ──────────────────────────────────────────────────────────────────────────────

def poly_add_fft(a: list, b: list) -> list:
    """a + b  (pointwise complex addition)."""
    return [x + y for x, y in zip(a, b)]


def poly_sub_fft(a: list, b: list) -> list:
    """a - b  (pointwise complex subtraction)."""
    return [x - y for x, y in zip(a, b)]


def poly_neg_fft(a: list) -> list:
    """-a  (pointwise complex negation)."""
    return [-x for x in a]


def poly_adj_fft(a: list) -> list:
    """
    Complex conjugate:  negate only the imaginary half.
    adj(a)[k] = conj(a[k])
    """
    n  = len(a)
    hn = n >> 1
    out = list(a)
    for u in range(hn, n):
        out[u] = -out[u]
    return out


def poly_mul_fft(a: list, b: list) -> list:
    """Pointwise complex multiplication  a * b."""
    n  = len(a)
    hn = n >> 1
    out = [0.0] * n
    for u in range(hn):
        a_re, a_im = a[u], a[u + hn]
        b_re, b_im = b[u], b[u + hn]
        out[u]      = a_re * b_re - a_im * b_im
        out[u + hn] = a_re * b_im + a_im * b_re
    return out


def poly_muladj_fft(a: list, b: list) -> list:
    """a * conj(b)  — multiply a by the conjugate of b."""
    n  = len(a)
    hn = n >> 1
    out = [0.0] * n
    for u in range(hn):
        a_re, a_im = a[u], a[u + hn]
        b_re, b_im = b[u], -b[u + hn]   # conjugate b
        out[u]      = a_re * b_re - a_im * b_im
        out[u + hn] = a_re * b_im + a_im * b_re
    return out


def poly_mulselfadj_fft(a: list) -> list:
    """
    a * conj(a) = |a|²  — result is real (imag = 0).
    Output: real parts = |a_k|², imag parts = 0.
    """
    n  = len(a)
    hn = n >> 1
    out = [0.0] * n
    for u in range(hn):
        out[u] = a[u] * a[u] + a[u + hn] * a[u + hn]
        # out[u + hn] = 0.0  (already zero)
    return out


def poly_mulconst_fft(a: list, c: float) -> list:
    """Multiply all slots by a real scalar c."""
    return [x * c for x in a]


def poly_mul_autoadj_fft(a: list, b: list) -> list:
    """
    Multiply by auto-adjoint b (b is real-only: b[hn..] = 0).
    Each complex a[k] is scaled by real b[k].
    """
    n  = len(a)
    hn = n >> 1
    out = [0.0] * n
    for u in range(hn):
        out[u]      = a[u]      * b[u]
        out[u + hn] = a[u + hn] * b[u]
    return out


def poly_div_autoadj_fft(a: list, b: list) -> list:
    """
    Divide by auto-adjoint b (b is real-only).
    Each complex a[k] is divided by real b[k].
    """
    n  = len(a)
    hn = n >> 1
    out = [0.0] * n
    for u in range(hn):
        out[u]      = a[u]      / b[u]
        out[u + hn] = a[u + hn] / b[u]
    return out


def poly_div_fft(a: list, b: list) -> list:
    """Pointwise complex division  a / b."""
    n  = len(a)
    hn = n >> 1
    out = [0.0] * n
    for u in range(hn):
        a_re, a_im = a[u], a[u + hn]
        b_re, b_im = b[u], b[u + hn]
        m = b_re * b_re + b_im * b_im     # |b|²
        out[u]      = (a_re * b_re + a_im * b_im) / m
        out[u + hn] = (a_im * b_re - a_re * b_im) / m
    return out


def poly_invnorm2_fft(a: list, b: list) -> list:
    """
    d[k] = 1 / (|a_k|² + |b_k|²)

    Result is real-only (auto-adjoint): d[0..hn-1] filled, d[hn..] = 0.
    Used in LDL and Gram-Schmidt computations.
    """
    n  = len(a)
    hn = n >> 1
    out = [0.0] * n
    for u in range(hn):
        norm2 = (a[u] * a[u] + a[u + hn] * a[u + hn]
               + b[u] * b[u] + b[u + hn] * b[u + hn])
        out[u] = 1.0 / norm2
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Split / Merge  (FFT domain)
# ──────────────────────────────────────────────────────────────────────────────
#
# These decompose/recompose a polynomial f(x) = f0(x²) + x·f1(x²)
# entirely in the FFT domain.  Used by ffsampling's recursive tree.
#
# Despite being called by a recursive tree, these operations themselves
# are non-recursive — they are simple loops over the FFT coefficients.
# ──────────────────────────────────────────────────────────────────────────────

def poly_split_fft(f: list) -> tuple:
    """
    Split f (n coefficients in FFT domain) into f0, f1 (each n/2).

    f(x) = f0(x²) + x·f1(x²)

    In FFT domain:
        f0[k] = (f[2k] + f[2k+1]) / 2
        f1[k] = (f[2k] - f[2k+1]) / 2 · conj(GM[k + n/2])

    Output: (f0, f1) each with n/2 floats in split-halves layout.
    """
    n  = len(f)
    hn = n >> 1
    qn = hn >> 1       # n/4
    tab = GM_TAB[n]

    f0 = [0.0] * hn
    f1 = [0.0] * hn

    # first element (u=0) — the formulas below use 2*u which starts at 0
    # but in the C code there's a special case for f0[0] = f[0], f1[0] = f[hn]
    # Let's handle it uniformly through the loop (u=0 maps to 2*0=0, 2*0+1=1).

    for u in range(qn):
        a_re = f[2 * u]
        a_im = f[2 * u + hn]
        b_re = f[2 * u + 1]
        b_im = f[2 * u + 1 + hn]

        # f0 = (a + b) / 2
        t_re = (a_re + b_re) * 0.5
        t_im = (a_im + b_im) * 0.5
        f0[u]      = t_re
        f0[u + qn] = t_im

        # f1 = (a - b) / 2 · conj(GM[u + hn])
        d_re = (a_re - b_re) * 0.5
        d_im = (a_im - b_im) * 0.5
        idx  = (u + hn) << 1
        s_re =  tab[idx]
        s_im = -tab[idx + 1]       # conjugate
        f1[u]      = d_re * s_re - d_im * s_im
        f1[u + qn] = d_re * s_im + d_im * s_re

    return f0, f1


def poly_merge_fft(f0: list, f1: list) -> list:
    """
    Merge f0, f1 (each n/2 in FFT domain) into f (n).

    Inverse of poly_split_fft.

    In FFT domain:
        b = f1[k] · GM[k + n_out/2]       (forward twiddle, NOT conjugate)
        f[2k]   = f0[k] + b
        f[2k+1] = f0[k] - b
    """
    hn  = len(f0)       # half of output size
    n   = hn << 1       # output size
    qn  = hn >> 1       # n/4 = number of complex values in each half
    tab = GM_TAB[n]

    f = [0.0] * n

    for u in range(qn):
        a_re = f0[u]
        a_im = f0[u + qn]

        # b = f1[u] · GM[u + hn]   (forward twiddle)
        idx  = (u + hn) << 1
        s_re = tab[idx]
        s_im = tab[idx + 1]        # NOT conjugate
        f1_re = f1[u]
        f1_im = f1[u + qn]
        b_re = f1_re * s_re - f1_im * s_im
        b_im = f1_re * s_im + f1_im * s_re

        # even slot: a + b
        f[2 * u]          = a_re + b_re
        f[2 * u + hn]     = a_im + b_im

        # odd slot: a - b
        f[2 * u + 1]      = a_re - b_re
        f[2 * u + 1 + hn] = a_im - b_im

    return f


# ──────────────────────────────────────────────────────────────────────────────
# LDL decomposition in FFT domain
# ──────────────────────────────────────────────────────────────────────────────

def poly_LDL_fft(g00: list, g01: list, g11: list) -> tuple:
    """
    In-place LDL decomposition of a 2×2 Gram matrix in FFT domain.

    Input (Gram matrix entries, all in FFT domain):
        G = [[g00, g01], [adj(g01), g11]]

    Output: (g00, g01, g11) modified so that:
        g00  = g00                        (unchanged — this is D00)
        g01  = g01 / g00                  (= L10 = off-diagonal of L)
        g11  = g11 - g01·adj(g01)·g00     (= D11)

    More precisely:
        L10  = g01 · adj(g00)⁻¹           = g01 / g00  (since g00 is auto-adjoint)
        D11  = g11 - L10 · adj(L10) · D00
    """
    n  = len(g00)
    hn = n >> 1

    d11 = [0.0] * n
    l10 = [0.0] * n

    for u in range(hn):
        g00_re = g00[u]
        g00_im = g00[u + hn]
        g01_re = g01[u]
        g01_im = g01[u + hn]
        g11_re = g11[u]
        g11_im = g11[u + hn]

        # g00 is auto-adjoint (real), so inv(g00) = 1/g00_re
        inv_g00 = 1.0 / g00_re

        # L10 = g01 / g00   (complex / real)
        l10_re = g01_re * inv_g00
        l10_im = g01_im * inv_g00
        l10[u]      = l10_re
        l10[u + hn] = l10_im

        # D11 = g11 - L10 · adj(L10) · D00
        #      = g11 - |L10|² · g00
        # Since g00 and g11 are auto-adjoint (real imag=0):
        l10_norm2 = l10_re * l10_re + l10_im * l10_im
        d11[u]      = g11_re - l10_norm2 * g00_re
        d11[u + hn] = 0.0      # auto-adjoint → imag = 0

    return g00, l10, d11


def poly_LDLmv_fft(d11: list, l10: list,
                    g00: list, g01: list, g11: list) -> None:
    """
    LDL decomposition variant — writes into provided output arrays
    (non-destructive to inputs g00, g01, g11).

    Same math as poly_LDL_fft but stores results in d11 and l10.
    g00 is unchanged (it is D00).
    """
    n  = len(g00)
    hn = n >> 1

    for u in range(hn):
        g00_re = g00[u]
        g01_re = g01[u]
        g01_im = g01[u + hn]
        g11_re = g11[u]

        inv_g00 = 1.0 / g00_re

        l_re = g01_re * inv_g00
        l_im = g01_im * inv_g00
        l10[u]      = l_re
        l10[u + hn] = l_im

        l_norm2 = l_re * l_re + l_im * l_im
        d11[u]      = g11_re - l_norm2 * g00_re
        d11[u + hn] = 0.0


# ──────────────────────────────────────────────────────────────────────────────
# Coefficient-domain convenience wrappers  (fft → operate → ifft)
# ──────────────────────────────────────────────────────────────────────────────

def poly_mul_coeff(f: list, g: list) -> list:
    """Multiply two polynomials in R[x]/(x^n+1) via FFT."""
    return ifft(poly_mul_fft(fft(f), fft(g)))


def poly_adj_coeff(f: list) -> list:
    """Adjoint of a polynomial in R[x]/(x^n+1) via FFT."""
    return ifft(poly_adj_fft(fft(f)))


def poly_div_coeff(f: list, g: list) -> list:
    """Divide two polynomials in R[x]/(x^n+1) via FFT."""
    return ifft(poly_div_fft(fft(f), fft(g)))


def poly_add_coeff(f: list, g: list) -> list:
    """Add two polynomials (coefficient representation)."""
    return [a + b for a, b in zip(f, g)]


# ──────────────────────────────────────────────────────────────────────────────
# Self-test
# ──────────────────────────────────────────────────────────────────────────────
def _test(n: int = 512) -> None:
    """Basic round-trip and operation tests."""
    import random
    random.seed(42)

    # ── Round-trip: ifft(fft(f)) ≈ f ──────────────────────────────────
    f = [random.gauss(0, 10) for _ in range(n)]
    f_fft = fft(f)
    f_back = ifft(f_fft)
    max_err = max(abs(a - b) for a, b in zip(f, f_back))
    assert max_err < 1e-6, f"[FAIL] n={n}: round-trip error {max_err}"

    # ── Split/Merge round-trip: merge(split(f_fft)) ≈ f_fft ──────────
    if n >= 8:  # split/merge need n >= 8 (qn >= 1)
        f0, f1 = poly_split_fft(f_fft)
        f_merged = poly_merge_fft(f0, f1)
        max_err = max(abs(a - b) for a, b in zip(f_fft, f_merged))
        assert max_err < 1e-6, f"[FAIL] n={n}: split/merge error {max_err}"

    # ── Multiplication: ifft(mul_fft(fft(f), fft(g))) ≈ f*g mod x^n+1
    g = [random.gauss(0, 10) for _ in range(n)]
    h_naive = [0.0] * n
    for i in range(n):
        for j in range(n):
            idx = (i + j) % n
            sign = -1 if (i + j) >= n else 1
            h_naive[idx] += sign * f[i] * g[j]

    h_fft  = poly_mul_fft(fft(f), fft(g))
    h_back = ifft(h_fft)
    max_err = max(abs(a - b) for a, b in zip(h_naive, h_back))
    assert max_err < 1e-3, f"[FAIL] n={n}: mul error {max_err}"

    # ── Adjoint: adj(adj(f)) ≈ f ─────────────────────────────────────
    f_adj = poly_adj_fft(f_fft)
    f_adj2 = poly_adj_fft(f_adj)
    max_err = max(abs(a - b) for a, b in zip(f_fft, f_adj2))
    assert max_err < 1e-10, f"[FAIL] n={n}: adj round-trip error {max_err}"

    # ── mulselfadj produces real (imag = 0) ──────────────────────────
    hn = n >> 1
    sa = poly_mulselfadj_fft(f_fft)
    max_imag = max(abs(sa[i]) for i in range(hn, n))
    assert max_imag < 1e-10, f"[FAIL] n={n}: mulselfadj imag {max_imag}"

    print(f"[OK]  n={n:4d} — round-trip, split/merge, mul, adj, mulselfadj all pass.")


if __name__ == "__main__":
    for _n in [8, 16, 64, 512, 1024]:
        _test(_n)