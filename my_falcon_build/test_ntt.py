"""
test_ntt.py — Comprehensive tests and evaluation for ntt.py
============================================================

Tests:
  1. Domain constants        — q, primitive root order, ψ correctness
  2. Montgomery arithmetic   — montymul, tomonty, frommonty round-trips
  3. Twiddle tables          — all evaluation points are distinct odd powers of ψ
  4. Round-trip (plain)      — intt(ntt(f)) == f
  5. Round-trip (Montgomery) — intt_montgomery(ntt_montgomery(f)) == f
  6. Cross-backend           — ntt(f) == ntt_montgomery(f)
  7. Convolution             — NTT poly_mul == schoolbook naive (random inputs)
  8. Edge cases              — zero poly, constant poly, x, x^(n-1)
  9. Linearity               — ntt(af+g) == a·ntt(f) + ntt(g)
 10. Negacyclic property     — x^n ≡ -1  =>  x*x^(n-1) ≡ -1 ≡ q-1
 11. Benchmark               — timed runs for n=512 and n=1024

Run: python test_ntt.py
"""

import random
import sys
import time

from ntt import (
    Q, G, LOGN, N_INV,
    _TABLES, _bitrev,
    mq_montymul, mq_add, mq_sub, mq_tomonty, mq_frommonty,
    ntt, intt,
    ntt_montgomery, intt_montgomery,
    ntt_mul, ntt_add, ntt_sub,
    poly_mul_ntt,
)

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _order(a: int, q: int) -> int:
    """Compute multiplicative order of a mod q."""
    o, cur = 1, a % q
    while cur != 1:
        cur = cur * a % q
        o += 1
    return o


def _poly_mul_naive(f: list, g: list, n: int, q: int) -> list:
    """Schoolbook negacyclic multiplication mod (x^n+1, q) — O(n²)."""
    h = [0] * n
    for i, fi in enumerate(f):
        for j, gj in enumerate(g):
            idx  = (i + j) % n
            sign = -1 if (i + j) >= n else 1
            h[idx] = (h[idx] + sign * fi * gj) % q
    return h


def _poly_add(f: list, g: list, q: int) -> list:
    return [(a + b) % q for a, b in zip(f, g)]


def _poly_scale(c: int, f: list, q: int) -> list:
    return [c * x % q for x in f]


PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"

_results = []  # list of (name, ok, detail)

def _record(name: str, ok: bool, detail: str = "") -> None:
    _results.append((name, ok, detail))
    status = PASS if ok else FAIL
    print(f"  [{status}]  {name}" + (f" — {detail}" if detail else ""))


# ─────────────────────────────────────────────────────────────────────────────
# 1. Domain constants
# ─────────────────────────────────────────────────────────────────────────────

def test_constants() -> None:
    print("\n── 1. Domain constants ──────────────────────────────────────────")

    # q must be an NTT-friendly prime
    _record("Q is prime (Fermat)",
            pow(2, Q - 1, Q) == 1 and Q == 12289)

    # G must be a primitive root mod Q
    ord_g = _order(G, Q)
    _record(f"ord(G={G}) = Q-1 = {Q-1}",
            ord_g == Q - 1,
            f"ord({G})={ord_g}")

    for n in (512, 1024):
        psi     = pow(G, (Q - 1) // (2 * n), Q)
        ord_psi = _order(psi, Q)
        negone  = pow(psi, n, Q)
        _record(f"n={n}: ord(ψ) = 2n = {2*n}",
                ord_psi == 2 * n,
                f"ψ={psi}, ord={ord_psi}")
        _record(f"n={n}: ψ^n = q-1 (i.e. -1)",
                negone == Q - 1,
                f"ψ^n={negone}")


# ─────────────────────────────────────────────────────────────────────────────
# 2. Montgomery arithmetic unit tests
# ─────────────────────────────────────────────────────────────────────────────

def test_montgomery() -> None:
    print("\n── 2. Montgomery arithmetic ─────────────────────────────────────")

    rng = random.Random(1)

    # round-trip: frommonty(tomonty(x)) == x
    vals = [rng.randrange(Q) for _ in range(1000)]
    ok   = all(mq_frommonty(mq_tomonty(x)) == x for x in vals)
    _record("tomonty → frommonty round-trip (1000 random)", ok)

    # montymul: x·y·R⁻¹ mod Q — verify against plain
    R_inv = pow(1 << 16, Q - 2, Q)
    samples = [(rng.randrange(Q), rng.randrange(Q)) for _ in range(1000)]
    ok = all(mq_montymul(x, y) == x * y * R_inv % Q for x, y in samples)
    _record("mq_montymul(x,y) == x·y·R⁻¹ mod Q (1000 pairs)", ok)

    # mq_add / mq_sub
    ok_add = all((mq_add(x, y) == (x + y) % Q) for x, y in samples)
    ok_sub = all((mq_sub(x, y) == (x - y) % Q) for x, y in samples)
    _record("mq_add correct (1000 pairs)", ok_add)
    _record("mq_sub correct (1000 pairs)", ok_sub)

    # identity in Montgomery domain: montymul(x_m, R mod Q) == x_m
    R_MOD_Q = (1 << 16) % Q
    ok = all(
        mq_montymul(mq_tomonty(x), mq_tomonty(R_MOD_Q)) == mq_tomonty(x * R_MOD_Q % Q)
        for x in vals[:200]
    )
    _record("Montgomery domain: mul by R_MOD_Q is consistent", ok)


# ─────────────────────────────────────────────────────────────────────────────
# 3. Twiddle-factor tables
# ─────────────────────────────────────────────────────────────────────────────

def test_twiddle_tables() -> None:
    print("\n── 3. Twiddle-factor tables ─────────────────────────────────────")

    for n in (512, 1024):
        psi  = pow(G, (Q - 1) // (2 * n), Q)
        gmb  = _TABLES[n]["gmb"]
        igmb = _TABLES[n]["igmb"]

        # gmb[k] = ψ^bitrev(k, logn)
        logn = LOGN[n]
        expected_gmb = [pow(psi, _bitrev(k, logn), Q) for k in range(n)]
        _record(f"n={n}: GMB entries match ψ^bitrev(k)", gmb == expected_gmb)

        # all evaluation points (gmb[n//2 .. n-1] are the leaves) are distinct
        # Actually, we check that all n entries in gmb are distinct
        _record(f"n={n}: GMB has {n} distinct values", len(set(gmb)) == n)

        # igmb[k] = (ψ⁻¹)^bitrev(k, logn)
        psi_inv = pow(psi, Q - 2, Q)
        expected_igmb = [pow(psi_inv, _bitrev(k, logn), Q) for k in range(n)]
        _record(f"n={n}: IGMB entries match ψ⁻¹^bitrev(k)", igmb == expected_igmb)

        # gmb[k] · igmb[k] ≡ 1 (they're inverses of each other)
        ok = all(gmb[k] * igmb[k] % Q == 1 for k in range(n))
        _record(f"n={n}: GMB[k] · IGMB[k] ≡ 1 for all k", ok)

        # Montgomery tables are just tomonty of plain tables
        gmb_m  = _TABLES[n]["gmb_m"]
        igmb_m = _TABLES[n]["igmb_m"]
        ok_gmb_m  = all(gmb_m[k]  == mq_tomonty(gmb[k])  for k in range(n))
        ok_igmb_m = all(igmb_m[k] == mq_tomonty(igmb[k]) for k in range(n))
        _record(f"n={n}: GMB_M = tomonty(GMB)", ok_gmb_m)
        _record(f"n={n}: IGMB_M = tomonty(IGMB)", ok_igmb_m)


# ─────────────────────────────────────────────────────────────────────────────
# 4 & 5. Round-trip tests
# ─────────────────────────────────────────────────────────────────────────────

def test_roundtrip() -> None:
    print("\n── 4–5. Round-trip (plain & Montgomery) ─────────────────────────")

    rng = random.Random(42)
    for n in (512, 1024):
        f = [rng.randrange(Q) for _ in range(n)]

        _record(f"n={n}: intt(ntt(f)) == f  [plain]",
                intt(ntt(f)) == f)
        _record(f"n={n}: intt_m(ntt_m(f)) == f  [Montgomery]",
                intt_montgomery(ntt_montgomery(f)) == f)


# ─────────────────────────────────────────────────────────────────────────────
# 6. Cross-backend consistency
# ─────────────────────────────────────────────────────────────────────────────

def test_cross_backend() -> None:
    print("\n── 6. Cross-backend consistency ─────────────────────────────────")

    rng = random.Random(7)
    for n in (512, 1024):
        f = [rng.randrange(Q) for _ in range(n)]
        _record(f"n={n}: ntt(f) == ntt_montgomery(f)",
                ntt(f) == ntt_montgomery(f))
        f_ntt = ntt(f)
        _record(f"n={n}: intt(f_ntt) == intt_montgomery(f_ntt)",
                intt(f_ntt) == intt_montgomery(f_ntt))


# ─────────────────────────────────────────────────────────────────────────────
# 7. Convolution correctness (random inputs)
# ─────────────────────────────────────────────────────────────────────────────

def test_convolution() -> None:
    print("\n── 7. Convolution correctness ───────────────────────────────────")

    rng = random.Random(99)
    for n in (512, 1024):
        for trial in range(5):
            f = [rng.randrange(Q) for _ in range(n)]
            g = [rng.randrange(Q) for _ in range(n)]
            expected = _poly_mul_naive(f, g, n, Q)
            got      = poly_mul_ntt(f, g)
            _record(f"n={n} trial {trial+1}: poly_mul_ntt == schoolbook",
                    got == expected)


# ─────────────────────────────────────────────────────────────────────────────
# 8. Edge cases
# ─────────────────────────────────────────────────────────────────────────────

def test_edge_cases() -> None:
    print("\n── 8. Edge cases ────────────────────────────────────────────────")

    for n in (512, 1024):
        zero = [0] * n
        one  = [1] + [0] * (n - 1)        # polynomial 1
        x    = [0, 1] + [0] * (n - 2)     # polynomial x
        xnm1 = [0] * (n - 1) + [1]        # polynomial x^(n-1)

        # zero polynomial
        _record(f"n={n}: ntt(0) == 0", ntt(zero) == zero)
        _record(f"n={n}: intt(0) == 0", intt(zero) == zero)

        # constant 1: NTT of 1 should equal the evaluation of the constant
        # polynomial 1 at all points = [1, 1, ..., 1]
        _record(f"n={n}: ntt([1,0,...,0]) == [1,1,...,1]",
                ntt(one) == [1] * n)

        # round-trip on specific polys
        _record(f"n={n}: round-trip of x",
                intt(ntt(x)) == x)
        _record(f"n={n}: round-trip of x^(n-1)",
                intt(ntt(xnm1)) == xnm1)

        # negacyclic: x^n ≡ -1.  So x · x^(n-1) = x^n ≡ q-1 (as constant poly)
        prod = poly_mul_ntt(x, xnm1)
        expected_xn = [Q - 1] + [0] * (n - 1)   # -1 = q-1 as constant polynomial
        _record(f"n={n}: x · x^(n-1) ≡ -1  (negacyclic property)",
                prod == expected_xn)

        # 1 · f == f
        rng = random.Random(13)
        f   = [rng.randrange(Q) for _ in range(n)]
        _record(f"n={n}: 1 · f == f",
                poly_mul_ntt(one, f) == f)


# ─────────────────────────────────────────────────────────────────────────────
# 9. Linearity of the NTT
# ─────────────────────────────────────────────────────────────────────────────

def test_linearity() -> None:
    print("\n── 9. Linearity ─────────────────────────────────────────────────")

    rng = random.Random(55)
    for n in (512, 1024):
        f = [rng.randrange(Q) for _ in range(n)]
        g = [rng.randrange(Q) for _ in range(n)]
        c = rng.randrange(1, Q)

        # ntt(f + g) == ntt(f) + ntt(g)  (pointwise in Z_q)
        lhs = ntt(_poly_add(f, g, Q))
        rhs = ntt_add(ntt(f), ntt(g))
        _record(f"n={n}: NTT(f+g) == NTT(f)+NTT(g)", lhs == rhs)

        # ntt(c·f) == c · ntt(f)
        lhs2 = ntt(_poly_scale(c, f, Q))
        rhs2 = [c * x % Q for x in ntt(f)]
        _record(f"n={n}: NTT(c·f) == c·NTT(f)", lhs2 == rhs2)


# ─────────────────────────────────────────────────────────────────────────────
# 10. Negacyclic convolution structure
# ─────────────────────────────────────────────────────────────────────────────

def test_negacyclic() -> None:
    print("\n── 10. Negacyclic structure ─────────────────────────────────────")

    for n in (512, 1024):
        # x^n mod (x^n + 1) = -1.  Verify directly.
        xn = [0] * (n + 1)
        xn[n] = 1            # x^n as length-(n+1) array
        # schoolbook reduce
        h = [0] * n
        for i in range(n + 1):
            if xn[i] == 0:
                continue
            idx  = i % n
            sign = -1 if i >= n else 1
            h[idx] = (h[idx] + sign * xn[i]) % Q
        # h should be [q-1, 0, 0, ...] = -1
        expected = [Q - 1] + [0] * (n - 1)
        _record(f"n={n}: x^n ≡ q-1 (negacyclic reduction)", h == expected)

        # Cauchy-product wrap-around: f*x^k shifts, negating wrapped terms
        rng = random.Random(77)
        f = [rng.randrange(Q) for _ in range(n)]
        # multiply by x (shift by 1, negate coefficient that wraps)
        f_x = [(-f[n - 1]) % Q] + f[:n - 1]
        got = poly_mul_ntt(f, [0, 1] + [0] * (n - 2))  # multiply by x
        _record(f"n={n}: f·x is correct cyclic shift with sign flip", got == f_x)


# ─────────────────────────────────────────────────────────────────────────────
# 11. Benchmark
# ─────────────────────────────────────────────────────────────────────────────

def run_benchmark(iters: int = 200) -> None:
    print("\n── 11. Benchmark ────────────────────────────────────────────────")

    for n in (512, 1024):
        rng = random.Random(0)
        f   = [rng.randrange(Q) for _ in range(n)]

        rows = []
        for label, fn in [
            ("ntt         (plain)",       ntt),
            ("intt        (plain)",       intt),
            ("ntt         (montgomery)",  ntt_montgomery),
            ("intt        (montgomery)",  intt_montgomery),
            ("poly_mul_ntt",              lambda p: poly_mul_ntt(p, f)),
        ]:
            t0  = time.perf_counter()
            for _ in range(iters):
                fn(f)
            ms = (time.perf_counter() - t0) / iters * 1000
            rows.append((label, ms))

        print(f"\n  n={n},  {iters} iterations:")
        print(f"  {'Function':<32}  {'ms/call':>9}")
        print(f"  {'-'*32}  {'-'*9}")
        for label, ms in rows:
            print(f"  {label:<32}  {ms:9.4f}")


# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

def _summary() -> None:
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
    print("NTT Test Suite")
    print("=" * 66)

    test_constants()
    test_montgomery()
    test_twiddle_tables()
    test_roundtrip()
    test_cross_backend()
    test_convolution()
    test_edge_cases()
    test_linearity()
    test_negacyclic()
    run_benchmark()

    _summary()
    sys.exit(0 if all(ok for _, ok, _ in _results) else 1)
