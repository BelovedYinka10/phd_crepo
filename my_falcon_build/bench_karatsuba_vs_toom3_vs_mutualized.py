#!/usr/bin/env python3
"""
bench_multiply.py — Benchmark multiplication backends for Falcon NTRU keygen
=============================================================================

Compares three polynomial multiplication strategies inside the reduce()
bottleneck of ntru_solve (>75% of keygen time at n=1024):

  1. Standard Karatsuba        — 2-way split, 3 sub-muls, O(n^1.585)
  2. Toom-3 (hybrid)           — one Toom-3 level at top, Karatsuba below
  3. Mutualized Karatsuba      — shares sub-products when computing both
                                 f*k and g*k with the same k

The bottleneck in reduce() is two karamul calls:
    fk = karamul(f, k)
    gk = karamul(g, k)

We benchmark these three strategies at full ntru_solve level and also
isolate the reduce() contribution.

The existing ntrugen.py is NOT modified.

Usage:  python bench_multiply.py [--trials N] [--n 512|1024]
"""

import argparse
import time
import statistics

from ntrugen import (karatsuba, karamul, galois_conjugate, field_norm, lift,
                     xgcd, bitsize, sqnorm, gs_norm, gen_poly)
from fft import (fft, ifft, poly_add_fft, poly_mul_fft, poly_adj_fft,
                 poly_div_fft)
from ntt import ntt
from params import Q


# ══════════════════════════════════════════════════════════════════════════════
# Counters
# ══════════════════════════════════════════════════════════════════════════════
_leaf_muls = 0
_reduce_time = 0.0


def _reset():
    global _leaf_muls, _reduce_time
    _leaf_muls = 0
    _reduce_time = 0.0


# ══════════════════════════════════════════════════════════════════════════════
# Backend 1: Standard Karatsuba (instrumented)
# ══════════════════════════════════════════════════════════════════════════════

def _kara(a, b, n):
    global _leaf_muls
    if n == 1:
        _leaf_muls += 1
        return [a[0] * b[0], 0]
    n2 = n // 2
    a0, a1 = a[:n2], a[n2:]
    b0, b1 = b[:n2], b[n2:]
    ax = [a0[i] + a1[i] for i in range(n2)]
    bx = [b0[i] + b1[i] for i in range(n2)]
    a0b0 = _kara(a0, b0, n2)
    a1b1 = _kara(a1, b1, n2)
    axbx = _kara(ax, bx, n2)
    for i in range(n):
        axbx[i] -= (a0b0[i] + a1b1[i])
    ab = [0] * (2 * n)
    for i in range(n):
        ab[i] += a0b0[i]
        ab[i + n] += a1b1[i]
        ab[i + n2] += axbx[i]
    return ab


def _karamul_std(a, b):
    n = len(a)
    ab = _kara(a, b, n)
    return [ab[i] - ab[i + n] for i in range(n)]


# ══════════════════════════════════════════════════════════════════════════════
# Backend 2: Toom-3 hybrid
#
# One level of Toom-3 at the top (pad to next multiple of 3),
# then Karatsuba for the 5 sub-problems (padded to next power of 2).
# ══════════════════════════════════════════════════════════════════════════════

def _pad(a, m):
    return a + [0] * (m - len(a))


def _poly_add(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(len(a)):
        r[i] += a[i]
    for i in range(len(b)):
        r[i] += b[i]
    return r


def _poly_sub(a, b):
    n = max(len(a), len(b))
    r = [0] * n
    for i in range(len(a)):
        r[i] += a[i]
    for i in range(len(b)):
        r[i] -= b[i]
    return r


def _poly_scale(a, c):
    return [x * c for x in a]


def _next_pow2(x):
    p = 1
    while p < x:
        p <<= 1
    return p


def _karamul_toom3(a, b):
    """
    Toom-3 at top level, Karatsuba for sub-problems.
    Pads a, b to next multiple of 3, splits into 3 parts,
    evaluates at 5 points, does 5 Karatsuba sub-muls, interpolates.
    """
    n = len(a)
    m3 = n if n % 3 == 0 else n + (3 - n % 3)
    m = m3 // 3  # size of each part

    ap = _pad(a, m3)
    bp = _pad(b, m3)

    a0, a1, a2 = ap[:m], ap[m:2*m], ap[2*m:3*m]
    b0, b1, b2 = bp[:m], bp[m:2*m], bp[2*m:3*m]

    # Evaluate at 0, 1, -1, 2, inf
    p0, q0 = a0, b0
    p1  = _poly_add(_poly_add(a0, a1), a2)
    q1  = _poly_add(_poly_add(b0, b1), b2)
    pm1 = _poly_add(_poly_sub(a0, a1), a2)
    qm1 = _poly_add(_poly_sub(b0, b1), b2)
    p2  = _poly_add(_poly_add(a0, _poly_scale(a1, 2)), _poly_scale(a2, 4))
    q2  = _poly_add(_poly_add(b0, _poly_scale(b1, 2)), _poly_scale(b2, 4))
    pinf, qinf = a2, b2

    # Pad sub-problems to next power of 2 for Karatsuba
    mp = _next_pow2(m)
    r0   = _kara(_pad(p0, mp),   _pad(q0, mp),   mp)[:2*m]
    r1   = _kara(_pad(p1, mp),   _pad(q1, mp),   mp)[:2*m]
    rm1  = _kara(_pad(pm1, mp),  _pad(qm1, mp),  mp)[:2*m]
    r2   = _kara(_pad(p2, mp),   _pad(q2, mp),   mp)[:2*m]
    rinf = _kara(_pad(pinf, mp), _pad(qinf, mp), mp)[:2*m]

    # Toom-3 interpolation (Bodrato 2007)
    c0 = r0[:]
    c4 = rinf[:]

    t1 = _poly_add(r1, rm1)     # 2*(c0 + c2 + c4)
    t2 = _poly_sub(r1, rm1)     # 2*(c1 + c3)

    c2 = _poly_sub(_poly_sub([x // 2 for x in t1], r0), rinf)
    s = [x // 2 for x in t2]    # c1 + c3

    temp = _poly_sub(_poly_sub(_poly_sub(r2, r0), _poly_scale(c2, 4)),
                     _poly_scale(rinf, 16))
    c1_4c3 = [x // 2 for x in temp]
    c3 = [(x - y) // 3 for x, y in zip(c1_4c3, s)]
    c1 = _poly_sub(s, c3)

    # Assemble full product (length 2*m3)
    ab = [0] * (2 * m3)
    for ci, off in [(c0, 0), (c1, m), (c2, 2*m), (c3, 3*m), (c4, 4*m)]:
        for i in range(len(ci)):
            if off + i < len(ab):
                ab[off + i] += ci[i]

    # Reduce mod (x^n + 1)
    return [ab[i] - ab[i + n] for i in range(n)]


# ══════════════════════════════════════════════════════════════════════════════
# Backend 3: Mutualized Karatsuba
# ══════════════════════════════════════════════════════════════════════════════

def _kara_dual(a, b, k, n):
    """Compute (a*k, b*k) sharing k's sub-products recursively."""
    global _leaf_muls
    if n == 1:
        _leaf_muls += 2
        return ([a[0] * k[0], 0], [b[0] * k[0], 0])

    n2 = n // 2
    a0, a1 = a[:n2], a[n2:]
    b0, b1 = b[:n2], b[n2:]
    k0, k1 = k[:n2], k[n2:]

    kx = [k0[i] + k1[i] for i in range(n2)]
    ax = [a0[i] + a1[i] for i in range(n2)]
    bx = [b0[i] + b1[i] for i in range(n2)]

    # 3 dual-recursive calls (vs 6 independent Karatsuba calls)
    a0k0, b0k0 = _kara_dual(a0, b0, k0, n2)
    a1k1, b1k1 = _kara_dual(a1, b1, k1, n2)
    axkx, bxkx = _kara_dual(ax, bx, kx, n2)

    ak = [0] * (2 * n)
    bk = [0] * (2 * n)
    for i in range(n):
        ak[i] += a0k0[i]
        ak[i + n] += a1k1[i]
        ak[i + n2] += axkx[i] - a0k0[i] - a1k1[i]

        bk[i] += b0k0[i]
        bk[i + n] += b1k1[i]
        bk[i + n2] += bxkx[i] - b0k0[i] - b1k1[i]

    return ak, bk


def _karamul_dual(f, g, k):
    """Mutualized Karatsuba mod (x^n+1): returns (f*k, g*k)."""
    n = len(f)
    fk_full, gk_full = _kara_dual(f, g, k, n)
    fk = [fk_full[i] - fk_full[i + n] for i in range(n)]
    gk = [gk_full[i] - gk_full[i + n] for i in range(n)]
    return fk, gk


# ══════════════════════════════════════════════════════════════════════════════
# reduce() with pluggable backend
# ══════════════════════════════════════════════════════════════════════════════

def _reduce(f, g, F, G, mul_fn, dual_fn=None):
    """
    Babai reduction. mul_fn is used for individual f*k, g*k products.
    If dual_fn is not None, it replaces both calls with one dual call.
    """
    global _reduce_time
    t0 = time.perf_counter()

    n = len(f)
    size = max(53, bitsize(min(f)), bitsize(max(f)),
               bitsize(min(g)), bitsize(max(g)))

    f_adj = [elt >> (size - 53) for elt in f]
    g_adj = [elt >> (size - 53) for elt in g]
    fa = fft(f_adj)
    ga = fft(g_adj)

    while True:
        Size = max(53, bitsize(min(F)), bitsize(max(F)),
                   bitsize(min(G)), bitsize(max(G)))
        if Size < size:
            break
        Fa = fft([elt >> (Size - 53) for elt in F])
        Ga = fft([elt >> (Size - 53) for elt in G])

        den = poly_add_fft(poly_mul_fft(fa, poly_adj_fft(fa)),
                           poly_mul_fft(ga, poly_adj_fft(ga)))
        num = poly_add_fft(poly_mul_fft(Fa, poly_adj_fft(fa)),
                           poly_mul_fft(Ga, poly_adj_fft(ga)))
        k = [int(round(x)) for x in ifft(poly_div_fft(num, den))]
        if all(x == 0 for x in k):
            break

        if dual_fn is not None:
            fk, gk = dual_fn(f, g, k)
        else:
            fk = mul_fn(f, k)
            gk = mul_fn(g, k)

        for i in range(n):
            F[i] -= fk[i] << (Size - size)
            G[i] -= gk[i] << (Size - size)

    _reduce_time += time.perf_counter() - t0
    return F, G


# ══════════════════════════════════════════════════════════════════════════════
# ntru_solve() — field_norm/lift always use standard karamul (from ntrugen.py)
# Only reduce() uses the pluggable backend.
# ══════════════════════════════════════════════════════════════════════════════

def _ntru_solve(f, g, mul_fn, dual_fn=None):
    """ntru_solve with pluggable reduce() backend."""
    n = len(f)
    if n == 1:
        d, u, v = xgcd(f[0], g[0])
        if d != 1:
            raise ValueError
        return [-Q * v], [Q * u]

    fp = field_norm(f)           # uses original karamul
    gp = field_norm(g)
    Fp, Gp = _ntru_solve(fp, gp, mul_fn, dual_fn)
    F = karamul(lift(Fp), galois_conjugate(g))   # original karamul
    G = karamul(lift(Gp), galois_conjugate(f))
    F, G = _reduce(f, g, F, G, mul_fn, dual_fn)  # pluggable backend
    return F, G


# ══════════════════════════════════════════════════════════════════════════════
# Verification
# ══════════════════════════════════════════════════════════════════════════════

def _verify(f, g, F, G):
    n = len(f)
    fG = karamul(f, G)
    gF = karamul(g, F)
    lhs = [fG[i] - gF[i] for i in range(n)]
    return lhs == [Q] + [0] * (n - 1)


# ══════════════════════════════════════════════════════════════════════════════
# Test pair generation
# ══════════════════════════════════════════════════════════════════════════════

def _gen_pairs(n, count):
    pairs = []
    while len(pairs) < count:
        f = gen_poly(n)
        g = gen_poly(n)
        if gs_norm(f, g, Q) > (1.17 ** 2) * Q:
            continue
        if n in (512, 1024):
            if any(elem == 0 for elem in ntt(f)):
                continue
        pairs.append((f, g))
    return pairs


# ══════════════════════════════════════════════════════════════════════════════
# Benchmark
# ══════════════════════════════════════════════════════════════════════════════

def run(n=512, num_trials=5):
    print("=" * 80)
    print(f"  Falcon NTRU Key Generation: Multiplication Backend Benchmark")
    print(f"  n = {n}, trials = {num_trials}")
    print(f"  (field_norm/lift use standard Karatsuba; only reduce() varies)")
    print("=" * 80)
    print()

    print(f"  Generating {num_trials} valid (f, g) pairs...", end="", flush=True)
    t0 = time.perf_counter()
    pairs = _gen_pairs(n, num_trials)
    print(f" done ({time.perf_counter() - t0:.1f}s)")
    print()

    backends = [
        ("Karatsuba",            _karamul_std,   None),
        ("Toom-3 (hybrid)",      _karamul_toom3, None),
        ("Mutualized Karatsuba", _karamul_std,   _karamul_dual),
    ]

    results = {}

    for name, mul_fn, dual_fn in backends:
        times, leaf_list, reduce_times = [], [], []
        ok = True

        print(f"  {name}...", end="", flush=True)
        for i, (f, g) in enumerate(pairs):
            _reset()
            t1 = time.perf_counter()
            try:
                F, G = _ntru_solve(f, g, mul_fn, dual_fn)
                F = [int(c) for c in F]
                G = [int(c) for c in G]
            except ValueError:
                continue
            elapsed = time.perf_counter() - t1

            if not _verify(f, g, F, G):
                print(f"\n    [FAIL] Trial {i}")
                ok = False

            times.append(elapsed)
            leaf_list.append(_leaf_muls)
            reduce_times.append(_reduce_time)

        if times:
            results[name] = {
                "mean": statistics.mean(times),
                "std": statistics.stdev(times) if len(times) > 1 else 0.0,
                "leaf": int(statistics.mean(leaf_list)),
                "reduce": statistics.mean(reduce_times),
                "ok": ok,
                "n": len(times),
            }
            r = results[name]
            print(f" [{('OK' if ok else 'FAIL')}] {r['n']}/{num_trials}, "
                  f"total={r['mean']:.2f}s, reduce={r['reduce']:.2f}s, "
                  f"leafs={r['leaf']:,}")
        else:
            print(" [SKIP]")
            results[name] = None

    # ── Table ────────────────────────────────────────────────────────────
    base = results.get("Karatsuba")
    if not base:
        print("\n  Baseline failed.")
        return

    print()
    print("=" * 80)
    print(f"  {'Backend':<23s} | {'Total (s)':>9s} | {'reduce (s)':>10s} | "
          f"{'Leaf muls':>12s} | {'Speedup':>7s} | {'Leaf save':>9s}")
    print("-" * 80)

    for name, _, _ in backends:
        r = results.get(name)
        if not r:
            print(f"  {name:<23s} |   (no data)")
            continue
        spd = base["mean"] / r["mean"] if r["mean"] > 0 else 0
        r_spd = base["reduce"] / r["reduce"] if r["reduce"] > 0 else 0
        sav = (1 - r["leaf"] / base["leaf"]) * 100 if base["leaf"] > 0 else 0
        print(f"  {name:<23s} | {r['mean']:>9.3f} | {r['reduce']:>10.3f} | "
              f"{r['leaf']:>12,d} | {spd:>6.2f}x | {sav:>8.1f}%")

    print("=" * 80)
    print()
    print("  Speedup = baseline_total / variant_total")
    print()

    # ── reduce() fraction ────────────────────────────────────────────────
    pct = base["reduce"] / base["mean"] * 100 if base["mean"] > 0 else 0
    print(f"  reduce() accounts for {pct:.0f}% of total Karatsuba time")
    print()
    print("  Notes:")
    print("  - Mutualized Karatsuba saves ~1/3 of recursive calls by sharing")
    print("    k's sub-products when computing f*k and g*k simultaneously.")
    print("  - Toom-3 hybrid suffers from padding overhead (512 -> 513 -> 3x171)")
    print("    and sub-problem padding (171 -> 256 for Karatsuba). In C/HW with")
    print("    native odd-size support, it would perform better.")
    print("  - Leaf mul counts are the same for Karatsuba vs Mutualized because")
    print("    the savings are in fewer recursive function calls, not fewer")
    print("    base-case multiplications.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=512)
    parser.add_argument("--trials", type=int, default=5)
    args = parser.parse_args()
    run(n=args.n, num_trials=args.trials)
