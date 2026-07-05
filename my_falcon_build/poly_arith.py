"""
poly_arith.py — Polynomial arithmetic over Z_q[x]/(x^n+1)
==========================================================

All polynomials are represented as plain Python lists of n integers in Z_q.

Operations
----------
  poly_add      : f + g  mod q               (coefficient-wise, no NTT)
  poly_sub      : f - g  mod q               (coefficient-wise, no NTT)
  poly_neg      : -f     mod q               (coefficient-wise, no NTT)
  poly_scale    : c·f    mod q               (scalar multiply, no NTT)
  poly_mul      : f · g  mod (x^n+1, q)     (via NTT — O(n log n))
  poly_inv      : f⁻¹    mod (x^n+1, q)     (via NTT — O(n log n))
  poly_sq_norm  : ||f||² = sum of squared coefficients (centered)
  poly_center   : map coefficients from [0,q) → (-q/2, q/2]

Reference: Falcon spec §2.2, §3.4
"""

__all__ = [
    "poly_add", "poly_sub", "poly_neg", "poly_scale",
    "poly_mul", "poly_inv",
    "poly_sq_norm", "poly_center",
]

from params import Q
from ntt import ntt, intt, ntt_mul, poly_mul_ntt


# ─────────────────────────────────────────────────────────────────────────────
# Coefficient-wise operations  (no NTT needed — already O(n))
# ─────────────────────────────────────────────────────────────────────────────

def poly_add(f: list, g: list) -> list:
    """f + g  mod q  (coefficient-wise)."""
    return [(a + b) % Q for a, b in zip(f, g)]


def poly_sub(f: list, g: list) -> list:
    """f - g  mod q  (coefficient-wise)."""
    return [(a - b) % Q for a, b in zip(f, g)]


def poly_neg(f: list) -> list:
    """-f  mod q  (coefficient-wise)."""
    return [(-a) % Q for a in f]


def poly_scale(c: int, f: list) -> list:
    """c · f  mod q  (scalar multiply, coefficient-wise)."""
    c = c % Q
    return [c * a % Q for a in f]


# ─────────────────────────────────────────────────────────────────────────────
# Polynomial multiplication  (via NTT)
# ─────────────────────────────────────────────────────────────────────────────

def poly_mul(f: list, g: list) -> list:
    """
    f · g  mod (x^n+1, q)  using NTT.

    Cost: 2 NTT + n pointwise multiplies + 1 iNTT = O(n log n)
    Alias for poly_mul_ntt from ntt.py — kept here so callers only
    need to import poly_arith.
    """
    return poly_mul_ntt(f, g)


# ─────────────────────────────────────────────────────────────────────────────
# Polynomial inversion  (via NTT)
# ─────────────────────────────────────────────────────────────────────────────

def poly_inv(f: list) -> list:
    """
    Compute f⁻¹  mod (x^n+1, q)  using NTT-domain inversion.

    Method
    ------
    The NTT of f gives its evaluations at n distinct points ωᵢ:
        F[i] = f(ωᵢ)  mod q

    If every F[i] ≠ 0  then f is invertible, and:
        (f⁻¹)(ωᵢ) = F[i]⁻¹  mod q

    So:
        f⁻¹ = iNTT( [pow(F[i], q-2, q) for each i] )

    This is O(n log n) — far better than the extended-Euclidean O(n²) approach.

    Raises
    ------
    ValueError  if f is not invertible (some NTT coefficient is 0 mod q).
    """
    F     = ntt(f)
    F_inv = []
    for i, v in enumerate(F):
        if v == 0:
            raise ValueError(
                f"poly_inv: polynomial is not invertible mod (x^n+1, {Q}) "
                f"— NTT coefficient {i} is 0"
            )
        F_inv.append(pow(v, Q - 2, Q))   # Fermat: v^(q-2) = v^(-1) mod q
    return intt(F_inv)


# ─────────────────────────────────────────────────────────────────────────────
# Norm and centering utilities
# ─────────────────────────────────────────────────────────────────────────────

def poly_center(f: list) -> list:
    """
    Map each coefficient from [0, q)  →  (-q/2, q/2].

    Falcon norms and Gram-Schmidt computations work with signed coefficients.
    This lifts the canonical representative from [0,q) into the symmetric
    range centred at 0.

    Rule:  if a > q//2  then a -= q
    """
    half = Q >> 1   # q // 2 = 6144
    return [a - Q if a > half else a for a in f]


def poly_sq_norm(f: list) -> int:
    """
    Squared Euclidean norm  ||f||²  with centered coefficients.

    ||f||² = Σ  fᵢ²   where fᵢ ∈ (-q/2, q/2]

    Used in:
      • Signature verification: ||s1||² + ||s2||² ≤ β²
      • Key generation: checking Gram-Schmidt norms
    """
    return sum(a * a for a in poly_center(f))
