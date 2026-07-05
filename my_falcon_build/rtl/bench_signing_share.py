#!/usr/bin/env python3
"""
bench_signing_share.py — BaseSampler's share of the SIGNING-dominant path.

The paper's "~30%" is 30% of *signing*, not 30% of samplerz.  Full sign() isn't
built in this repo yet (keygen done; sign/verify next), so we use the honest
proxy: ffsampling_fft over a real Falcon-512 LDL tree.  ffsampling IS the
dominant cost of signing — the surrounding hash / s=c-z*B / compression are
comparatively cheap — so

    BaseSampler% of ffsampling   ~~   BaseSampler% of signing.

Method:
  1. keygen(512, "tree")  -> normalized LDL tree T  (the thing sign uses)
  2. random FFT-domain target t = [t0, t1]
  3. wrap samplerz.basesampler with a timer to accumulate its total time
  4. run ffsampling_fft N times; report BaseSampler time / ffsampling time

Each ffsampling pass over n=512 issues n=512 samplerz calls (one per tree leaf
coordinate), each averaging ~1.7 basesampler calls — a signing-scale workload.

Usage:  python3 bench_signing_share.py [--passes P] [--n 512|1024]
"""
import argparse
import os
import sys
import time
import random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import samplerz                                  # noqa: E402
from params import FALCON_PARAMS                 # noqa: E402
from keygen import keygen                        # noqa: E402
from ffsampling import ffsampling_fft            # noqa: E402


def random_fft_target(n, rng):
    """A plausible FFT-domain target [t0, t1] (split-halves layout, length n)."""
    return [[rng.gauss(0.0, 1.0) for _ in range(n)] for _ in range(2)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--passes", type=int, default=200,
                    help="number of ffsampling_fft passes to time (default 200)")
    ap.add_argument("--n", type=int, default=512, choices=(512, 1024))
    args = ap.parse_args()
    n = args.n

    print(f"Signing-share benchmark  (Falcon-{n}, {args.passes} ffsampling passes)")
    print("=" * 66)
    print("keygen (building LDL tree) ...", flush=True)
    sk, _ = keygen(n, mode="tree")
    sigmin = FALCON_PARAMS[n]["sigma_min"]
    rng = random.Random(0xFA1C0)

    # --- instrument basesampler with a timer ------------------------------
    acc = {"t": 0.0, "calls": 0}
    orig = samplerz.basesampler

    def timed_basesampler(*a, **k):
        t0 = time.perf_counter()
        r = orig(*a, **k)
        acc["t"] += time.perf_counter() - t0
        acc["calls"] += 1
        return r

    # warm up (also triggers any import-time costs)
    samplerz.basesampler = timed_basesampler
    ffsampling_fft(random_fft_target(n, rng), sk.T, sigmin, os.urandom)
    acc["t"] = 0.0
    acc["calls"] = 0

    # --- timed loop -------------------------------------------------------
    t_start = time.perf_counter()
    for _ in range(args.passes):
        t = random_fft_target(n, rng)
        ffsampling_fft(t, sk.T, sigmin, os.urandom)
    t_ffs = time.perf_counter() - t_start
    samplerz.basesampler = orig

    base_t = acc["t"]
    share = base_t / t_ffs
    per_pass = t_ffs / args.passes
    calls_per_pass = acc["calls"] / args.passes

    print("\nResults")
    print("-" * 66)
    print(f"  ffsampling_fft per pass     : {per_pass * 1e3:9.3f} ms")
    print(f"  basesampler calls per pass  : {calls_per_pass:9.1f}")
    print(f"  total ffsampling time       : {t_ffs:9.3f} s")
    print(f"  total BaseSampler time      : {base_t:9.3f} s")
    print(f"  ---------------------------------------------")
    print(f"  BaseSampler share of ffsampling : {share * 100:5.1f} %")
    print("\n  (this is the honest proxy for the paper's 'BaseSampler ~ 30% of")
    print("   signing'; the surrounding hash/recombine/compress steps that the")
    print("   full sign() will add are comparatively cheap, so the true")
    print("   signing-share will be slightly LOWER than the number above.)")
    print("\n  NOTE: the per-call timer adds a little overhead to BaseSampler,")
    print("        so this share is a mild OVER-estimate.")


if __name__ == "__main__":
    main()
