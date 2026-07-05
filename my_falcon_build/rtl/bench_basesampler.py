#!/usr/bin/env python3
"""
bench_basesampler.py — software BaseSampler throughput vs. a hardware projection.

What this measures (real):
  * software throughput of basesampler()           [samples/sec on THIS machine]
  * software throughput of full samplerz()         [samples/sec]
  * the fraction of samplerz() time spent in BaseSampler  (the "~30%" claim)

What this projects (estimate, NOT a measurement):
  * hardware throughput of base_sampler.v.  The RTL core is combinational and
    emits one z0 per clock, so HW throughput = clock frequency.  We print the
    projection at several Fmax values.  Real Fmax requires synthesis (Yosys /
    Vivado / Quartus) — until then these are ceilings, not results.

HONESTY CAVEAT printed at the end:
  Python is ~50-100x slower than the optimized C the literature benchmarks
  against, so a Python-vs-hardware ratio OVERSTATES the true speedup.  Use this
  to understand the *shape* of the bottleneck, not to quote a headline number.

Usage:
  python3 bench_basesampler.py [--n N] [--fmax MHZ [MHZ ...]]
"""
import argparse
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import samplerz                                   # noqa: E402
from samplerz import basesampler, samplerz as samplerz_fn  # noqa: E402

# Falcon-512 sampler parameters (params.py): sigma in (sigmin, MAX_SIGMA).
SIGMIN = 1.2778336969128337     # Falcon-512 sigma_min
SIGMA  = 1.55                   # a representative per-leaf sigma
MU     = 0.5                    # a representative center


def time_calls(fn, n):
    """Return (seconds_total, seconds_per_call) for n invocations of fn()."""
    fn()                                 # warm up
    t0 = time.perf_counter()
    for _ in range(n):
        fn()
    dt = time.perf_counter() - t0
    return dt, dt / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=200_000,
                    help="iterations per benchmark (default 200000)")
    ap.add_argument("--fmax", type=float, nargs="+", default=[100, 250, 500],
                    help="projected hardware clock(s) in MHz (default 100 250 500)")
    args = ap.parse_args()
    n = args.n

    print(f"BaseSampler benchmark  (n = {n:,} iterations)")
    print("=" * 64)

    # ---- 1. software BaseSampler ------------------------------------------
    _, t_base = time_calls(basesampler, n)
    sw_base_rate = 1.0 / t_base
    print(f"\n[software] basesampler()")
    print(f"  per call : {t_base * 1e9:10.1f} ns")
    print(f"  rate     : {sw_base_rate / 1e6:10.3f} M samples/sec")

    # ---- 2. software full samplerz(), and BaseSampler's share -------------
    # Count how many basesampler() calls each samplerz() makes (rejection loop).
    call_count = {"n": 0}
    orig = samplerz.basesampler

    def counting_basesampler(*a, **k):
        call_count["n"] += 1
        return orig(*a, **k)

    samplerz.basesampler = counting_basesampler
    try:
        _, t_samp = time_calls(
            lambda: samplerz_fn(MU, SIGMA, SIGMIN), n)
    finally:
        samplerz.basesampler = orig

    avg_iters = call_count["n"] / (n + 1)        # +1 for the warm-up call
    base_share = (avg_iters * t_base) / t_samp
    print(f"\n[software] samplerz()  (mu={MU}, sigma={SIGMA})")
    print(f"  per call         : {t_samp * 1e9:10.1f} ns")
    print(f"  rate             : {1e-6 / t_samp:10.3f} M samples/sec")
    print(f"  basesampler calls/sample : {avg_iters:6.3f}  (rejection loop)")
    print(f"  BaseSampler share of samplerz time : {base_share * 100:5.1f} %")

    # ---- 3. hardware projection -------------------------------------------
    print(f"\n[hardware projection] base_sampler.v  (1 sample / clock)")
    print(f"  {'Fmax':>8}   {'HW rate':>16}   {'vs this Python sw':>20}")
    for mhz in args.fmax:
        hw_rate = mhz * 1e6                       # 1 sample per cycle
        speedup = hw_rate / sw_base_rate
        print(f"  {mhz:6.0f}MHz   {hw_rate / 1e6:10.1f} Msps   "
              f"{speedup:14,.0f}x")

    # ---- caveats ----------------------------------------------------------
    print("\n" + "=" * 64)
    print("CAVEATS")
    print("  * HW numbers are PROJECTIONS (1 sample/clk). Real Fmax needs")
    print("    synthesis; the 72-bit compare -> popcount path sets the ceiling.")
    print("  * Python sw is ~50-100x slower than optimized C. The 'vs Python'")
    print("    column therefore OVERSTATES the true HW speedup -- treat it as")
    print("    an upper bound, not a headline figure.")
    print("  * The honest takeaway is the SHAPE: BaseSampler is a large, fully")
    print(f"    parallelisable slice (~{base_share*100:.0f}% of samplerz here), and in HW it")
    print("    collapses to one combinational cone -> one sample per cycle.")


if __name__ == "__main__":
    main()
