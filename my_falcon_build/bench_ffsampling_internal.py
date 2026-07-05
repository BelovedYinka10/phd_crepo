"""
bench_ffsampling_internal.py — RDTSC drill-down inside ffsampling_fft
=======================================================================

The outer bench_rdtsc.py treats ffsampling_fft as a black box.
This script breaks it open:

  Per-operation accumulators:
    samplerz        — leaf Gaussian draws  (called 2n times per sign)
    poly_split      — halve FFT array       (called 2×(n-1) times, all levels)
    poly_merge      — double FFT array      (called 2×(n-1) times, all levels)
    l10_mul         — poly_mul_fft(residual, L10)   (n-1 per sign)
    l10_add_sub     — poly_add_fft + poly_sub_fft   (2×(n-1) per sign)

  Per-depth breakdown:
    Timing accumulated by tree depth (0 = root, log2(n)-1 = leaf level)

Usage:
  python bench_ffsampling_internal.py [--n 512|1024] [--trials 30]
"""

import argparse
import ctypes
import os
import subprocess
import sys
import tempfile
import statistics
from math import log2
from collections import defaultdict

# ── RDTSC ─────────────────────────────────────────────────────────────────────

_RDTSC_C = """\
#include <stdint.h>
uint64_t rdtsc_s(void) {
    uint64_t lo, hi;
    __asm__ __volatile__("mfence\\n\\trdtsc" : "=a"(lo), "=d"(hi));
    return (hi << 32) | lo;
}
uint64_t rdtsc_e(void) {
    uint64_t lo, hi;
    unsigned int aux;
    __asm__ __volatile__("rdtscp" : "=a"(lo), "=d"(hi), "=c"(aux));
    __asm__ __volatile__("mfence");
    return (hi << 32) | lo;
}
"""

def _build_rdtsc():
    td  = tempfile.mkdtemp(prefix="ffsamp_rdtsc_")
    src = os.path.join(td, "r.c")
    lib = os.path.join(td, "r.so")
    with open(src, "w") as f:
        f.write(_RDTSC_C)
    try:
        r = subprocess.run(["cc", "-O2", "-shared", "-fPIC", "-o", lib, src],
                           capture_output=True, timeout=10)
        if r.returncode != 0:
            raise RuntimeError(r.stderr.decode())
        dll = ctypes.CDLL(lib)
        dll.rdtsc_s.restype = ctypes.c_uint64
        dll.rdtsc_e.restype = ctypes.c_uint64
        return dll.rdtsc_s, dll.rdtsc_e, "hardware RDTSC"
    except Exception as e:
        import time
        try:
            hz = int(subprocess.run(["sysctl", "-n", "hw.cpufrequency_max"],
                                    capture_output=True, text=True, timeout=3).stdout.strip())
        except Exception:
            hz = 2_600_000_000
        ghz = hz / 1e9
        print(f"[warn] RDTSC compile failed ({e}); using perf_counter_ns × {ghz:.2f} GHz")
        def _fake():
            return int(time.perf_counter_ns() * ghz)
        return _fake, _fake, f"perf_counter_ns×{ghz:.2f}"


_rs, _re, RDTSC_SRC = _build_rdtsc()

sys.path.insert(0, os.path.dirname(__file__))

from fft import (fft, ifft, poly_split_fft, poly_merge_fft,
                 poly_add_fft, poly_sub_fft, poly_mul_fft)
from samplerz import samplerz as _samplerz_orig
from ffsampling import ldl_fft, ffldl_fft, normalize_tree
from params import FALCON_PARAMS

CPU_HZ = 2_600_000_000
try:
    CPU_HZ = int(subprocess.run(["sysctl", "-n", "hw.cpufrequency_max"],
                                capture_output=True, text=True, timeout=3).stdout.strip())
except Exception:
    pass
CPU_GHZ = CPU_HZ / 1e9


# ── Instrumented ffsampling ────────────────────────────────────────────────────
#
# We reimplement ffsampling_fft here with RDTSC wrappers around each primitive.
# Accumulators are global dicts filled per signature, then appended to lists
# for statistical analysis.

class _Acc:
    """Holds per-signature cycle tallies and per-depth timing."""
    def __init__(self, max_depth):
        self.ops   = defaultdict(int)   # op_name → cycles this signature
        self.depth = defaultdict(int)   # depth   → cycles this signature
        self.max_depth = max_depth


def _ffsampling_instrumented(t, T, sigmin, randombytes, acc, depth=0):
    """
    Drop-in replacement for ffsampling_fft that accumulates RDTSC cycle counts
    into `acc` (an _Acc instance).

    Operations timed:
      samplerz    — leaf Gaussian draws
      poly_split  — poly_split_fft calls (halving the FFT array)
      poly_merge  — poly_merge_fft calls (doubling back)
      l10_mulsub  — poly_mul_fft(residual, l10) + poly_sub_fft + poly_add_fft
    """
    n = len(t[0])
    z = [None, None]
    d0 = _rs()

    if len(T) == 3 and len(T[1]) == 2:
        # ── Leaf level ─────────────────────────────────────────────────────
        l10    = T[0]
        sigma1 = T[2][0]
        sigma0 = T[1][0]

        # samplerz for z[1]
        t0 = _rs()
        z1_re = _samplerz_orig(t[1][0], sigma1, sigmin, randombytes)
        t1 = _re()
        acc.ops["samplerz"] += t1 - t0

        z[1] = [float(z1_re), 0.0]

        # L10 adjustment (scalar ops — no FFT, just 2-element arithmetic)
        t0 = _rs()
        diff_re = t[1][0] - z[1][0]
        diff_im = t[1][1] - z[1][1]
        l_re, l_im = l10[0], l10[1]
        adj_re = diff_re * l_re - diff_im * l_im
        adj_im = diff_re * l_im + diff_im * l_re
        t0b_re = t[0][0] + adj_re
        t0b_im = t[0][1] + adj_im
        t1 = _re()
        acc.ops["leaf_l10_adj"] += t1 - t0

        # samplerz for z[0]
        t0 = _rs()
        z0_re = _samplerz_orig(t0b_re, sigma0, sigmin, randombytes)
        t1 = _re()
        acc.ops["samplerz"] += t1 - t0

        z[0] = [float(z0_re), 0.0]

    else:
        # ── Internal node ───────────────────────────────────────────────────
        l10, T0, T1 = T

        # split t[1]
        t0 = _rs()
        t1_split = poly_split_fft(t[1])
        t1 = _re()
        acc.ops["poly_split"] += t1 - t0

        # recurse T1
        z1_halves = _ffsampling_instrumented(t1_split, T1, sigmin, randombytes, acc, depth + 1)

        # merge z[1]
        t0 = _rs()
        z[1] = poly_merge_fft(*z1_halves)
        t1 = _re()
        acc.ops["poly_merge"] += t1 - t0

        # compute t0_adjusted = t[0] + (t[1] - z[1]) * l10
        t0 = _rs()
        residual = poly_sub_fft(t[1], z[1])
        correction = poly_mul_fft(residual, l10)
        t0b = poly_add_fft(t[0], correction)
        t1 = _re()
        acc.ops["l10_mulsub"] += t1 - t0

        # split t0b
        t0 = _rs()
        t0b_split = poly_split_fft(t0b)
        t1 = _re()
        acc.ops["poly_split"] += t1 - t0

        # recurse T0
        z0_halves = _ffsampling_instrumented(t0b_split, T0, sigmin, randombytes, acc, depth + 1)

        # merge z[0]
        t0 = _rs()
        z[0] = poly_merge_fft(*z0_halves)
        t1 = _re()
        acc.ops["poly_merge"] += t1 - t0

    d1 = _re()
    acc.depth[depth] += d1 - d0
    return z


# ── Formatting ────────────────────────────────────────────────────────────────

def fmt_cyc(c):
    if c >= 1e9: return f"{c/1e9:8.3f} Gc"
    if c >= 1e6: return f"{c/1e6:8.3f} Mc"
    if c >= 1e3: return f"{c/1e3:8.1f} Kc"
    return f"{c:8.0f}  c"

def fmt_us(c):
    us = c / CPU_HZ * 1e6
    if us >= 1e3: return f"{us/1e3:8.2f} ms"
    return f"{us:8.1f} µs"


# ── Main benchmark ─────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n",      type=int, default=512, choices=[512, 1024])
    ap.add_argument("--trials", type=int, default=20)
    args = ap.parse_args()

    n      = args.n
    trials = args.trials
    p      = FALCON_PARAMS[n]
    sigmin = p["sigma_min"]
    sigma  = p["sigma"]
    max_depth = int(log2(n))    # depth 0 = root, max_depth-1 = leaf level

    print(f"\nffsampling_fft internal RDTSC profiler  — Falcon-{n}")
    print(f"  Timer  : {RDTSC_SRC}")
    print(f"  CPU    : {CPU_GHZ:.3f} GHz")
    print(f"  Trials : {trials}")
    print(f"  Tree depth: {max_depth} levels  (2n = {2*n} leaf samplerz calls/sign)")

    # ── Build one LDL tree (keygen-lite: skip ntru_gen for speed) ────────────
    print("\nBuilding LDL tree (using ntrugen, please wait) …", end="", flush=True)
    from ntrugen import ntru_gen
    from ffsampling import gram
    f, g, F, G = ntru_gen(n)
    neg_f, neg_F = [-c for c in f], [-c for c in F]
    B = [[g, neg_f], [G, neg_F]]
    G_mat = gram(B)
    G_fft = [[fft(G_mat[i][j]) for j in range(2)] for i in range(2)]
    T = ffldl_fft(G_fft)
    normalize_tree(T, sigma)
    print(" done.")

    # ── Sample random targets (flat random floats, good enough for timing) ───
    import random
    rng = random.Random(42)

    # per-op and per-depth accumulators across all trials
    all_ops   = defaultdict(list)   # op_name → [total_cycles_per_trial]
    all_depth = defaultdict(list)   # depth   → [total_cycles_per_trial]
    total_cyc = []

    for _ in range(trials):
        t0v = [rng.gauss(0, 1) for _ in range(n)]
        t1v = [rng.gauss(0, 1) for _ in range(n)]

        acc = _Acc(max_depth)
        ta = _rs()
        _ffsampling_instrumented([t0v, t1v], T, sigmin, os.urandom, acc, depth=0)
        tb = _re()
        total_cyc.append(tb - ta)

        for op, cyc in acc.ops.items():
            all_ops[op].append(cyc)
        for d, cyc in acc.depth.items():
            all_depth[d].append(cyc)

    # ── Report: per-operation breakdown ─────────────────────────────────────
    OP_ORDER = ["samplerz", "l10_mulsub", "poly_split", "poly_merge", "leaf_l10_adj"]
    med_total = statistics.median(total_cyc)

    print(f"\n  ── Per-operation breakdown (median over {trials} trials) ──")
    W = max(len(k) for k in OP_ORDER) + 2
    print(f"  {'Operation':<{W}}  {'Median':>12}  {'µs':>10}  {'Min':>12}  {'Max':>12}  {'%':>6}  {'calls/sign':>10}")
    print(f"  {'-'*W}  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*6}  {'-'*10}")

    # Expected call counts per sign for n-degree:
    #   samplerz    : 2n     (one per coefficient, two per leaf node)
    #   poly_split  : 2*(n-1) (split t[1] and t[0b] at every internal node)
    #   poly_merge  : 2*(n-1)
    #   l10_mulsub  : n-1    (one per internal node, covers sub+mul+add)
    #   leaf_l10_adj: n      (one per leaf node pair? actually n/2*2 = n leaves ... )
    # Actually:
    #   internal nodes at depth d: 2^d nodes, each has 2 children
    #   leaf nodes (bottom internal, len(T[1])==2): 2^(log2(n)-1) = n/2 nodes, 2 samplerz each → n
    #   BUT ffsampling_fft is called recursively on halved arrays, and T root covers 2 dimensions
    #   Root: processes [t0, t1] each of length n → calls itself twice on length n/2
    #   Each level doubles calls; leaves call samplerz 2x → total 2n samplerz calls
    call_counts = {
        "samplerz":     2 * n,
        "l10_mulsub":   n - 1,
        "poly_split":   2 * (n - 1),
        "poly_merge":   2 * (n - 1),
        "leaf_l10_adj": n,
    }

    accounted = 0
    for op in OP_ORDER:
        if op not in all_ops:
            continue
        s   = all_ops[op]
        med = statistics.median(s)
        mn  = min(s)
        mx  = max(s)
        pct = 100.0 * med / med_total
        cnt = call_counts.get(op, "?")
        per = f"{med/cnt:.0f} c" if isinstance(cnt, int) else "?"
        print(f"  {op:<{W}}  {fmt_cyc(med)}  {fmt_us(med)}  {fmt_cyc(mn)}  {fmt_cyc(mx)}  {pct:5.1f}%  {cnt:>5}× @ {per}")
        accounted += med

    unaccounted = med_total - accounted
    unaccounted_pct = 100.0 * unaccounted / med_total
    print(f"  {'overhead/recursion':<{W}}  {fmt_cyc(unaccounted)}  {fmt_us(unaccounted)}  {'':>12}  {'':>12}  {unaccounted_pct:5.1f}%")
    print(f"  {'TOTAL (wall)':<{W}}  {fmt_cyc(med_total)}  {fmt_us(med_total)}")

    # ── Bottleneck callout ────────────────────────────────────────────────────
    meds_ops = {op: statistics.median(all_ops[op]) for op in OP_ORDER if op in all_ops}
    worst_op = max(meds_ops, key=meds_ops.get)
    print(f"\n  *** Bottleneck inside ffsampling: [{worst_op}]  {100*meds_ops[worst_op]/med_total:.1f}% ***")

    # ── Per-depth breakdown ───────────────────────────────────────────────────
    print(f"\n  ── Per-depth wall time (depth 0 = root, {max_depth-1} = leaf) ──")
    print(f"  {'Depth':>6}  {'n at level':>10}  {'Median':>12}  {'µs':>10}  {'%':>6}")
    print(f"  {'-'*6}  {'-'*10}  {'-'*12}  {'-'*10}  {'-'*6}")
    depth_sum = 0
    for d in range(max_depth):
        if d not in all_depth:
            continue
        s   = all_depth[d]
        med = statistics.median(s)
        n_at_d = n >> d           # array length at this depth
        pct = 100.0 * med / med_total
        depth_sum += med
        print(f"  {d:>6}  {n_at_d:>10}  {fmt_cyc(med)}  {fmt_us(med)}  {pct:5.1f}%")

    print(f"\n  samplerz throughput: "
          f"{2*n / (med_total / CPU_HZ) / 1e6:.2f} M calls/sec  "
          f"({(med_total / CPU_HZ * 1e9 / (2*n)):.0f} ns/call)")
    print(f"  Full ffsampling throughput: "
          f"{1.0 / (med_total / CPU_HZ):.1f} calls/sec  "
          f"({med_total / CPU_HZ * 1000:.2f} ms/sign)\n")


if __name__ == "__main__":
    main()
