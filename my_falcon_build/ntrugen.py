"""
ntrugen.py — NTRU key generation for Falcon
=============================================

Solves the NTRU equation:  f*G - g*F = q  mod (x^n + 1)
to produce the secret key quadruple (f, g, F, G).

Algorithm (Falcon spec Section 3.8.2, Algorithm 5 NTRUGen)
----------------------------------------------------------
  1. Sample small (f, g) with discrete Gaussian coefficients
  2. Check Gram-Schmidt norm of the NTRU lattice basis
  3. Check f is invertible mod q via NTT
  4. Recursively solve the NTRU equation (NTRUSolve)
  5. Return (f, g, F, G)

NTRUSolve uses field norms to reduce degree n → n/2 → ... → 1,
solves the base case via extended GCD, then lifts back up with
Babai reduction at each level to keep (F, G) small.

Reference
---------
  Falcon spec v1.2 Section 3.8.2
  Reference Python implementation by Thomas Prest (PQShield)
"""

__all__ = [
    "ntru_gen", "ntru_solve", "gen_poly",
    "karatsuba", "karamul",
    "galois_conjugate", "field_norm", "lift",
    "reduce", "gs_norm", "xgcd", "sqnorm",
]

from fft import (fft, ifft,
                 poly_add_fft, poly_mul_fft, poly_adj_fft, poly_div_fft,
                 poly_mul_coeff, poly_adj_coeff, poly_div_coeff,
                 poly_add_coeff)
from ntt import ntt
from samplerz import samplerz
from params import Q


# ──────────────────────────────────────────────────────────────────────────────
# Integer polynomial arithmetic (over Z, NOT Z_q)
# ──────────────────────────────────────────────────────────────────────────────

def karatsuba(a, b, n):
    """
    Karatsuba multiplication of two polynomials over Z.
    a, b have n coefficients each (n must be a power of 2).
    Returns a list of 2n coefficients (product, unreduced).
    """
    if n == 1:
        return [a[0] * b[0], 0]
    n2 = n // 2
    a0, a1 = a[:n2], a[n2:]
    b0, b1 = b[:n2], b[n2:]
    ax = [a0[i] + a1[i] for i in range(n2)]
    bx = [b0[i] + b1[i] for i in range(n2)]
    a0b0 = karatsuba(a0, b0, n2)
    a1b1 = karatsuba(a1, b1, n2)
    axbx = karatsuba(ax, bx, n2)
    for i in range(n):
        axbx[i] -= (a0b0[i] + a1b1[i])
    ab = [0] * (2 * n)
    for i in range(n):
        ab[i] += a0b0[i]
        ab[i + n] += a1b1[i]
        ab[i + n2] += axbx[i]
    return ab


def karamul(a, b):
    """
    Karatsuba multiplication mod (x^n + 1) over Z.
    Uses x^n = -1 for reduction: result[i] = ab[i] - ab[i+n].
    """
    n = len(a)
    ab = karatsuba(a, b, n)
    return [ab[i] - ab[i + n] for i in range(n)]


# ──────────────────────────────────────────────────────────────────────────────
# Ring algebra helpers
# ──────────────────────────────────────────────────────────────────────────────

def galois_conjugate(a):
    """
    Galois conjugate: a(x) -> a(-x) in Q[x]/(x^n+1).
    Negates all odd-index coefficients.
    """
    return [a[i] if i % 2 == 0 else -a[i] for i in range(len(a))]


def field_norm(a):
    """
    Field norm: project from Q[x]/(x^n+1) down to Q[x]/(x^(n/2)+1).
    Separates even/odd coefficients, computes a_e^2 - x * a_o^2.
    """
    n2 = len(a) // 2
    ae = [a[2 * i] for i in range(n2)]
    ao = [a[2 * i + 1] for i in range(n2)]
    ae_sq = karamul(ae, ae)
    ao_sq = karamul(ao, ao)
    res = ae_sq[:]
    # Subtract x * ao_sq:  x * p(x) mod (x^(n/2)+1)
    #   coeff[0] += ao_sq[n2-1]  (wraparound with sign flip)
    #   coeff[i] -= ao_sq[i-1]   for i = 1..n2-1
    for i in range(n2 - 1):
        res[i + 1] -= ao_sq[i]
    res[0] += ao_sq[n2 - 1]
    return res


def lift(a):
    """
    Lift from Q[x]/(x^(n/2)+1) to Q[x]/(x^n+1).
    Maps a(x) -> a(x^2): place coefficients at even indices.
    """
    n = len(a)
    res = [0] * (2 * n)
    for i in range(n):
        res[2 * i] = a[i]
    return res


# ──────────────────────────────────────────────────────────────────────────────
# Utilities
# ──────────────────────────────────────────────────────────────────────────────

def bitsize(a):
    """
    Bit length of integer |a|, rounded up to the next multiple of 8.
    Used by reduce() to decide precision scaling.
    """
    val = abs(a)
    res = 0
    while val:
        res += 8
        val >>= 8
    return res


def sqnorm(v):
    """
    Squared Euclidean norm of a vector of polynomials.
    v is a list of coefficient lists.
    """
    return sum(c ** 2 for poly in v for c in poly)


def xgcd(b, n):
    """
    Extended GCD of integers b and n.
    Returns (d, u, v) such that d = u*b + v*n and d = gcd(b, n).
    """
    x0, x1, y0, y1 = 1, 0, 0, 1
    while n != 0:
        q_div, b, n = b // n, n, b % n
        x0, x1 = x1, x0 - q_div * x1
        y0, y1 = y1, y0 - q_div * y1
    return b, x0, y0


# ──────────────────────────────────────────────────────────────────────────────
# Babai reduction (Algorithm 7)
# ──────────────────────────────────────────────────────────────────────────────

def reduce(f, g, F, G):
    """
    Babai reduction of (F, G) relative to (f, g).
    (F, G) <- (F, G) - k*(f, g)  where
      k = round( (F*adj(f) + G*adj(g)) / (f*adj(f) + g*adj(g)) )
    computed in finite-precision FFT domain.
    """
    n = len(f)
    size = max(53, bitsize(min(f)), bitsize(max(f)),
               bitsize(min(g)), bitsize(max(g)))

    f_adjust = [elt >> (size - 53) for elt in f]
    g_adjust = [elt >> (size - 53) for elt in g]
    fa_fft = fft(f_adjust)
    ga_fft = fft(g_adjust)

    while True:
        Size = max(53, bitsize(min(F)), bitsize(max(F)),
                   bitsize(min(G)), bitsize(max(G)))
        if Size < size:
            break

        F_adjust = [elt >> (Size - 53) for elt in F]
        G_adjust = [elt >> (Size - 53) for elt in G]
        Fa_fft = fft(F_adjust)
        Ga_fft = fft(G_adjust)

        den_fft = poly_add_fft(poly_mul_fft(fa_fft, poly_adj_fft(fa_fft)),
                               poly_mul_fft(ga_fft, poly_adj_fft(ga_fft)))
        num_fft = poly_add_fft(poly_mul_fft(Fa_fft, poly_adj_fft(fa_fft)),
                               poly_mul_fft(Ga_fft, poly_adj_fft(ga_fft)))
        k_fft = poly_div_fft(num_fft, den_fft)
        k = ifft(k_fft)
        k = [int(round(elt)) for elt in k]
        if all(elt == 0 for elt in k):
            break
        fk = karamul(f, k)
        gk = karamul(g, k)
        for i in range(n):
            F[i] -= fk[i] << (Size - size)
            G[i] -= gk[i] << (Size - size)
    return F, G


# ──────────────────────────────────────────────────────────────────────────────
# NTRU solver (recursive)
# ──────────────────────────────────────────────────────────────────────────────

def ntru_solve(f, g):
    """
    Recursively solve the NTRU equation: find F, G in Z[x]/(x^n+1)
    such that f*G - g*F = q.
    Reduces degree n -> n/2 via field norms until base case n=1.
    """
    n = len(f)
    if n == 1:
        f0, g0 = f[0], g[0]
        d, u, v = xgcd(f0, g0)
        if d != 1:
            raise ValueError("gcd(f[0], g[0]) != 1, NTRU equation unsolvable")
        return [-Q * v], [Q * u]
    else:
        fp = field_norm(f)
        gp = field_norm(g)
        Fp, Gp = ntru_solve(fp, gp)
        F = karamul(lift(Fp), galois_conjugate(g))
        G = karamul(lift(Gp), galois_conjugate(f))
        F, G = reduce(f, g, F, G)
        return F, G


# ──────────────────────────────────────────────────────────────────────────────
# Gram-Schmidt norm check
# ──────────────────────────────────────────────────────────────────────────────

def gs_norm(f, g, q):
    """
    Squared Gram-Schmidt norm of the NTRU matrix [[g, -f], [G, -F]].
    Equivalent to line 9 of Algorithm 5 (NTRUGen).
    """
    sqnorm_fg = sqnorm([f, g])
    ffgg = poly_add_coeff(poly_mul_coeff(f, poly_adj_coeff(f)),
                          poly_mul_coeff(g, poly_adj_coeff(g)))
    Ft = poly_div_coeff(poly_adj_coeff(g), ffgg)
    Gt = poly_div_coeff(poly_adj_coeff(f), ffgg)
    sqnorm_FG = (q ** 2) * sqnorm([Ft, Gt])
    return max(sqnorm_fg, sqnorm_FG)


# ──────────────────────────────────────────────────────────────────────────────
# Polynomial sampling and key generation
# ──────────────────────────────────────────────────────────────────────────────

def gen_poly(n):
    """
    Generate a degree-(n-1) polynomial with coefficients from
    D_{Z, 0, sigma_fg}  where sigma_fg = 1.17 * sqrt(q / (2n)).

    Uses the sum-of-Gaussians trick: sample 4096 values from
    D_{Z, 0, sigma_base}, then sum groups of (4096/n) to get the
    correct overall standard deviation.
    """
    sigma = 1.43300980528773   # 1.17 * sqrt(12289 / 8192)
    assert n < 4096
    f0 = [samplerz(0, sigma, sigma - 0.001) for _ in range(4096)]
    f = [0] * n
    k = 4096 // n
    for i in range(n):
        f[i] = sum(f0[i * k + j] for j in range(k))
    return f


def ntru_gen(n):
    """
    NTRU key generation (Algorithm 5, NTRUGen).
    Returns (f, g, F, G) in Z[x]/(x^n+1) satisfying f*G - g*F = q.

    Only supports n in {512, 1024} (NTT invertibility check).
    """
    while True:
        f = gen_poly(n)
        g = gen_poly(n)
        # Check Gram-Schmidt norm
        if gs_norm(f, g, Q) > (1.17 ** 2) * Q:
            continue
        # Check f is invertible in Z_q[x]/(x^n+1)
        f_ntt = ntt(f)
        if any(elem == 0 for elem in f_ntt):
            continue
        # Solve NTRU equation
        try:
            F, G = ntru_solve(f, g)
            F = [int(coef) for coef in F]
            G = [int(coef) for coef in G]
            return f, g, F, G
        except ValueError:
            continue
