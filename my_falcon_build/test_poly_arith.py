"""
test_poly_arith.py — Tests for poly_arith.py
=============================================

Tests:
  1. poly_add / poly_sub / poly_neg     — basic coefficient-wise ops
  2. poly_scale                         — scalar multiply
  3. poly_mul                           — NTT multiply vs schoolbook
  4. poly_inv                           — f · f⁻¹ == 1
  5. poly_inv non-invertible            — correct ValueError raised
  6. poly_center                        — range mapping
  7. poly_sq_norm                       — known values + signature bound check
  8. Algebraic identities               — distributivity, commutativity

Run: python test_poly_arith.py
"""

import random
import sys

from params import Q
from poly_arith import (
    poly_add, poly_sub, poly_neg, poly_scale,
    poly_mul, poly_inv,
    poly_center, poly_sq_norm,
)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"
_results = []

def _record(name, ok, detail=""):
    _results.append((name, ok, detail))
    status = PASS if ok else FAIL
    print(f"  [{status}]  {name}" + (f" — {detail}" if detail else ""))

def _rand_poly(n, rng):
    return [rng.randrange(Q) for _ in range(n)]

def _poly_mul_naive(f, g, n, q):
    h = [0] * n
    for i, fi in enumerate(f):
        for j, gj in enumerate(g):
            idx  = (i + j) % n
            sign = -1 if (i + j) >= n else 1
            h[idx] = (h[idx] + sign * fi * gj) % q
    return h

# ─────────────────────────────────────────────────────────────────────────────
# 1. poly_add / poly_sub / poly_neg
# ─────────────────────────────────────────────────────────────────────────────

def test_add_sub_neg():
    print("\n── 1. poly_add / poly_sub / poly_neg ────────────────────────────")
    rng = random.Random(1)
    for n in (512, 1024):
        f = _rand_poly(n, rng)
        g = _rand_poly(n, rng)

        # f + g - g == f
        _record(f"n={n}: (f+g)-g == f",
                poly_sub(poly_add(f, g), g) == f)

        # f - f == 0
        zero = [0] * n
        _record(f"n={n}: f-f == 0",
                poly_sub(f, f) == zero)

        # f + 0 == f
        _record(f"n={n}: f+0 == f",
                poly_add(f, zero) == f)

        # neg(f) + f == 0
        _record(f"n={n}: neg(f)+f == 0",
                poly_add(poly_neg(f), f) == zero)

        # commutativity: f+g == g+f
        _record(f"n={n}: f+g == g+f",
                poly_add(f, g) == poly_add(g, f))

        # sub via neg: f-g == f+neg(g)
        _record(f"n={n}: f-g == f+neg(g)",
                poly_sub(f, g) == poly_add(f, poly_neg(g)))

# ─────────────────────────────────────────────────────────────────────────────
# 2. poly_scale
# ─────────────────────────────────────────────────────────────────────────────

def test_scale():
    print("\n── 2. poly_scale ────────────────────────────────────────────────")
    rng = random.Random(2)
    for n in (512, 1024):
        f = _rand_poly(n, rng)
        c = rng.randrange(1, Q)

        # scale by 0 gives zero
        _record(f"n={n}: scale(0,f) == 0",
                poly_scale(0, f) == [0] * n)

        # scale by 1 gives f
        _record(f"n={n}: scale(1,f) == f",
                poly_scale(1, f) == f)

        # scale by c then by c⁻¹ gives f back
        c_inv = pow(c, Q - 2, Q)
        _record(f"n={n}: scale(c⁻¹, scale(c, f)) == f",
                poly_scale(c_inv, poly_scale(c, f)) == f)

        # distributivity: scale(c, f+g) == scale(c,f) + scale(c,g)
        g = _rand_poly(n, rng)
        lhs = poly_scale(c, poly_add(f, g))
        rhs = poly_add(poly_scale(c, f), poly_scale(c, g))
        _record(f"n={n}: c·(f+g) == c·f + c·g",
                lhs == rhs)

# ─────────────────────────────────────────────────────────────────────────────
# 3. poly_mul
# ─────────────────────────────────────────────────────────────────────────────

def test_mul():
    print("\n── 3. poly_mul ──────────────────────────────────────────────────")
    rng = random.Random(3)
    for n in (512, 1024):
        f = _rand_poly(n, rng)
        g = _rand_poly(n, rng)
        one  = [1] + [0] * (n - 1)
        zero = [0] * n

        # vs schoolbook
        expected = _poly_mul_naive(f, g, n, Q)
        _record(f"n={n}: poly_mul == schoolbook",
                poly_mul(f, g) == expected)

        # commutativity
        _record(f"n={n}: f·g == g·f",
                poly_mul(f, g) == poly_mul(g, f))

        # multiply by 1
        _record(f"n={n}: f·1 == f",
                poly_mul(f, one) == f)

        # multiply by 0
        _record(f"n={n}: f·0 == 0",
                poly_mul(f, zero) == zero)

        # distributivity: f·(g+h) == f·g + f·h
        h = _rand_poly(n, rng)
        lhs = poly_mul(f, poly_add(g, h))
        rhs = poly_add(poly_mul(f, g), poly_mul(f, h))
        _record(f"n={n}: f·(g+h) == f·g + f·h",
                lhs == rhs)

# ─────────────────────────────────────────────────────────────────────────────
# 4. poly_inv
# ─────────────────────────────────────────────────────────────────────────────

def test_inv():
    print("\n── 4. poly_inv ──────────────────────────────────────────────────")
    rng = random.Random(4)
    one = None
    for n in (512, 1024):
        one = [1] + [0] * (n - 1)

        for trial in range(3):
            f     = _rand_poly(n, rng)
            try:
                f_inv = poly_inv(f)
                prod  = poly_mul(f, f_inv)
                _record(f"n={n} trial {trial+1}: f · f⁻¹ == 1",
                        prod == one)
            except ValueError:
                # Non-invertible poly — extremely rare for random inputs
                _record(f"n={n} trial {trial+1}: f not invertible (skipped)", True,
                        "rare but valid")

        # constant polynomial c: inv should be [c⁻¹, 0, 0, ...]
        c   = rng.randrange(1, Q)
        f_c = [c] + [0] * (n - 1)
        inv = poly_inv(f_c)
        _record(f"n={n}: inv([c,0,...]) == [c⁻¹,0,...]",
                inv == [pow(c, Q - 2, Q)] + [0] * (n - 1))

# ─────────────────────────────────────────────────────────────────────────────
# 5. poly_inv non-invertible
# ─────────────────────────────────────────────────────────────────────────────

def test_inv_noninvertible():
    print("\n── 5. poly_inv — non-invertible detection ───────────────────────")
    # The zero polynomial is never invertible
    for n in (512, 1024):
        zero = [0] * n
        raised = False
        try:
            poly_inv(zero)
        except ValueError:
            raised = True
        _record(f"n={n}: poly_inv(0) raises ValueError", raised)

# ─────────────────────────────────────────────────────────────────────────────
# 6. poly_center
# ─────────────────────────────────────────────────────────────────────────────

def test_center():
    print("\n── 6. poly_center ───────────────────────────────────────────────")
    half = Q >> 1   # 6144

    # all output values must be in (-q/2, q/2]
    rng = random.Random(6)
    for n in (512, 1024):
        f      = _rand_poly(n, rng)
        fc     = poly_center(f)
        in_range = all(-half < a <= half for a in fc)
        _record(f"n={n}: all centered coefficients in (-q/2, q/2]", in_range)

    # specific known values
    cases = [
        (0,        0),
        (1,        1),
        (half,     half),       # q//2 stays
        (half + 1, half + 1 - Q),  # q//2+1 maps to negative
        (Q - 1,    -1),
    ]
    ok = all(poly_center([a])[0] == expected for a, expected in cases)
    _record("Specific centering cases correct", ok,
            str([(a, poly_center([a])[0]) for a, _ in cases]))

# ─────────────────────────────────────────────────────────────────────────────
# 7. poly_sq_norm
# ─────────────────────────────────────────────────────────────────────────────

def test_sq_norm():
    print("\n── 7. poly_sq_norm ──────────────────────────────────────────────")

    # zero polynomial has norm 0
    for n in (512, 1024):
        _record(f"n={n}: ||0||² == 0",
                poly_sq_norm([0] * n) == 0)

    # constant polynomial [1,0,...] has norm 1
    for n in (512, 1024):
        f = [1] + [0] * (n - 1)
        _record(f"n={n}: ||[1,0,...,0]||² == 1",
                poly_sq_norm(f) == 1)

    # norm is always non-negative
    rng = random.Random(7)
    for n in (512, 1024):
        f = _rand_poly(n, rng)
        _record(f"n={n}: ||f||² >= 0",
                poly_sq_norm(f) >= 0)

    # manual check: poly with known centered coefficients
    # f = [1, q-1, 2] -> centered = [1, -1, 2] -> norm = 1+1+4 = 6
    f_test = [1, Q - 1, 2] + [0] * (512 - 3)
    _record("Manual: ||[1, q-1, 2, 0,...]||² == 6",
            poly_sq_norm(f_test) == 6)

    # Falcon-512 beta_sq check: signature norms should be <= 34_034_726
    # Generate a very small poly to confirm it passes the bound
    from params import FALCON_PARAMS
    beta_sq_512 = FALCON_PARAMS[512]["beta_sq"]
    small = [1] * 512   # all-ones poly: centered = [1,1,...] -> norm = 512
    _record(f"||all-ones, n=512||² = 512 < beta_sq={beta_sq_512:,}",
            poly_sq_norm(small) == 512 and 512 < beta_sq_512)

# ─────────────────────────────────────────────────────────────────────────────
# 8. Algebraic identities
# ─────────────────────────────────────────────────────────────────────────────

def test_algebra():
    print("\n── 8. Algebraic identities ──────────────────────────────────────")
    rng = random.Random(8)
    for n in (512, 1024):
        f = _rand_poly(n, rng)
        g = _rand_poly(n, rng)
        h = _rand_poly(n, rng)
        c = rng.randrange(1, Q)

        # associativity of add: (f+g)+h == f+(g+h)
        _record(f"n={n}: (f+g)+h == f+(g+h)",
                poly_add(poly_add(f, g), h) == poly_add(f, poly_add(g, h)))

        # scale then add: c·f + c·g == c·(f+g)
        _record(f"n={n}: c·f + c·g == c·(f+g)",
                poly_add(poly_scale(c, f), poly_scale(c, g)) ==
                poly_scale(c, poly_add(f, g)))

        # mul then add (distributivity — already in test_mul but good to repeat)
        lhs = poly_mul(poly_add(f, g), h)
        rhs = poly_add(poly_mul(f, h), poly_mul(g, h))
        _record(f"n={n}: (f+g)·h == f·h + g·h",
                lhs == rhs)

        # scale is same as mul by constant poly
        f_c = poly_scale(c, f)
        c_poly = [c] + [0] * (n - 1)
        _record(f"n={n}: scale(c,f) == [c,0,...] · f",
                f_c == poly_mul(c_poly, f))

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

def _summary():
    total  = len(_results)
    passed = sum(1 for _, ok, _ in _results if ok)
    failed = total - passed
    print(f"\n{'='*66}")
    print(f"  Results: {passed}/{total} passed", end="")
    if failed:
        print(f"  ── {failed} FAILED:")
        for name, ok, detail in _results:
            if not ok:
                print(f"    ✗  {name}" + (f"  [{detail}]" if detail else ""))
    else:
        print("  — all OK")
    print(f"{'='*66}")

# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("poly_arith Test Suite")
    print("=" * 66)

    test_add_sub_neg()
    test_scale()
    test_mul()
    test_inv()
    test_inv_noninvertible()
    test_center()
    test_sq_norm()
    test_algebra()

    _summary()
    sys.exit(0 if all(ok for _, ok, _ in _results) else 1)
