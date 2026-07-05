"""
test_ntrugen.py — Tests for NTRU key generation module
=======================================================
"""

import random
import time
from ntrugen import (karatsuba, karamul, galois_conjugate, field_norm, lift,
                     bitsize, xgcd, reduce, ntru_solve, sqnorm, gs_norm,
                     gen_poly, ntru_gen)
from ntt import ntt
from params import Q

random.seed(42)


# ── Helper: schoolbook multiplication mod (x^n + 1) over Z ──────────────────

def schoolbook_mul(a, b):
    """Naive O(n^2) multiplication mod (x^n+1) over Z, for reference."""
    n = len(a)
    c = [0] * n
    for i in range(n):
        for j in range(n):
            idx = (i + j) % n
            sign = -1 if (i + j) >= n else 1
            c[idx] += sign * a[i] * b[j]
    return c


def check_ntru_equation(f, g, F, G, q):
    """Verify f*G - g*F = [q, 0, 0, ..., 0] over Z[x]/(x^n+1)."""
    n = len(f)
    fG = karamul(f, G)
    gF = karamul(g, F)
    lhs = [fG[i] - gF[i] for i in range(n)]
    expected = [q] + [0] * (n - 1)
    return lhs == expected


# ── Tests ────────────────────────────────────────────────────────────────────

def test_karatsuba():
    """Karatsuba matches schoolbook (before mod reduction)."""
    for n in [1, 2, 4, 8, 16]:
        a = [random.randint(-10, 10) for _ in range(n)]
        b = [random.randint(-10, 10) for _ in range(n)]
        # Schoolbook: full product (degree 2n-2)
        full = [0] * (2 * n)
        for i in range(n):
            for j in range(n):
                full[i + j] += a[i] * b[j]
        # Karatsuba returns 2n coefficients
        kara = karatsuba(a, b, n)
        assert kara == full, f"n={n}: karatsuba mismatch"
    print("[OK]  karatsuba")


def test_karamul():
    """karamul matches schoolbook mod (x^n+1)."""
    for n in [2, 4, 8, 16, 32]:
        a = [random.randint(-10, 10) for _ in range(n)]
        b = [random.randint(-10, 10) for _ in range(n)]
        expected = schoolbook_mul(a, b)
        result = karamul(a, b)
        assert result == expected, f"n={n}: karamul mismatch"
    print("[OK]  karamul")


def test_galois_conjugate():
    """conj(conj(a)) == a, and conj negates odd coefficients."""
    for n in [4, 8, 16]:
        a = [random.randint(-100, 100) for _ in range(n)]
        ac = galois_conjugate(a)
        # Odd coefficients negated
        for i in range(n):
            if i % 2 == 0:
                assert ac[i] == a[i]
            else:
                assert ac[i] == -a[i]
        # Involution
        acc = galois_conjugate(ac)
        assert acc == a
    print("[OK]  galois_conjugate")


def test_field_norm_and_lift():
    """field_norm produces n/2 coefficients; algebraic identity holds."""
    for n in [4, 8, 16, 32]:
        a = [random.randint(-5, 5) for _ in range(n)]
        fn = field_norm(a)
        assert len(fn) == n // 2, f"n={n}: field_norm wrong length"

        # Algebraic identity: N(a) = a * conj(a) evaluated at x^2
        # i.e., karamul(a, conj(a)) should have zero odd coefficients
        # and even coefficients should match field_norm
        prod = karamul(a, galois_conjugate(a))
        for i in range(n):
            if i % 2 == 1:
                assert prod[i] == 0, f"n={n}: odd coeff {i} nonzero: {prod[i]}"
        even_coeffs = [prod[2 * i] for i in range(n // 2)]
        assert even_coeffs == fn, f"n={n}: field_norm identity failed"

    # Lift round-trip
    for n in [4, 8, 16]:
        a = [random.randint(-10, 10) for _ in range(n)]
        la = lift(a)
        assert len(la) == 2 * n
        for i in range(n):
            assert la[2 * i] == a[i]
            assert la[2 * i + 1] == 0
    print("[OK]  field_norm / lift")


def test_xgcd():
    """Extended GCD: u*b + v*n = d for various pairs."""
    test_cases = [
        (12289, 1),
        (35, 15),
        (97, 37),
        (12289, 7),
        (1024, 513),
    ]
    for b, n in test_cases:
        d, u, v = xgcd(b, n)
        assert d == u * b + v * n, f"xgcd({b},{n}): identity failed"
        # Verify it's actually the GCD
        from math import gcd
        assert d == gcd(b, n), f"xgcd({b},{n}): d={d} != gcd={gcd(b, n)}"
    print("[OK]  xgcd")


def test_bitsize():
    """bitsize rounds up to multiple of 8."""
    assert bitsize(0) == 0
    assert bitsize(1) == 8
    assert bitsize(127) == 8
    assert bitsize(128) == 8
    assert bitsize(255) == 8
    assert bitsize(256) == 16
    assert bitsize(-256) == 16
    print("[OK]  bitsize")


def test_sqnorm():
    """Squared norm of vector of polynomials."""
    assert sqnorm([[1, 2, 3]]) == 14
    assert sqnorm([[1, 0], [0, 1]]) == 2
    assert sqnorm([[3, 4]]) == 25
    print("[OK]  sqnorm")


def test_ntru_solve_small():
    """Solve NTRU equation at small dimensions, verify f*G - g*F = q."""
    for n in [4, 8, 16]:
        # Try multiple random (f, g) until one works (gcd=1 at base)
        solved = False
        for _ in range(200):
            f = [random.randint(-3, 3) for _ in range(n)]
            g = [random.randint(-3, 3) for _ in range(n)]
            if f[0] == 0 or g[0] == 0:
                continue
            try:
                F, G = ntru_solve(f, g)
                assert check_ntru_equation(f, g, F, G, Q), \
                    f"n={n}: NTRU equation failed"
                solved = True
                break
            except ValueError:
                continue
        assert solved, f"n={n}: could not find solvable (f, g) in 200 tries"
    print("[OK]  ntru_solve (small n)")


def test_gs_norm():
    """gs_norm returns a positive float for valid polynomials."""
    n = 8
    f = [random.randint(-3, 3) for _ in range(n)]
    g = [random.randint(-3, 3) for _ in range(n)]
    # Make sure f, g aren't all zeros
    f[0] = 1
    g[0] = 1
    norm = gs_norm(f, g, Q)
    assert norm > 0, "gs_norm should be positive"
    assert isinstance(norm, float) or isinstance(norm, int)
    print("[OK]  gs_norm")


def test_ntru_gen_512():
    """Full integration test: generate NTRU keys at n=512."""
    print("       ntru_gen(512) running (may take a few seconds)...", end="", flush=True)
    t0 = time.time()
    f, g, F, G = ntru_gen(512)
    dt = time.time() - t0

    # Verify lengths
    assert len(f) == len(g) == len(F) == len(G) == 512

    # Verify NTRU equation
    assert check_ntru_equation(f, g, F, G, Q), "NTRU equation failed for n=512"

    # Verify f is NTT-invertible
    f_ntt = ntt(f)
    assert all(elem != 0 for elem in f_ntt), "f not invertible"

    # Verify Gram-Schmidt norm
    norm = gs_norm(f, g, Q)
    assert norm <= (1.17 ** 2) * Q, f"GS norm too large: {norm}"

    print(f"\r[OK]  ntru_gen(512) — {dt:.1f}s")


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    test_karatsuba()
    test_karamul()
    test_galois_conjugate()
    test_field_norm_and_lift()
    test_xgcd()
    test_bitsize()
    test_sqnorm()
    test_ntru_solve_small()
    test_gs_norm()
    test_ntru_gen_512()
    print(f"\nAll tests passed.")
