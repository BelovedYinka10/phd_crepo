"""
bench_sign.py  —  Falcon.py reference: true CPU-cycle benchmark.

Uses RDTSC (via rdtsc.so) for hardware cycle counting — reads the CPU
timestamp counter directly, more accurate than wall_time × freq.

Note: Intel TSC is invariant (fixed at base freq = 2.60 GHz regardless
of Turbo). These are *reference cycles*, not *core cycles*. To get core
cycles you need `perf stat` (Linux PMC).

Usage:
  python bench_sign.py          # both n=512 and n=1024
  python bench_sign.py 512      # just n=512
"""
import sys
import ctypes
import statistics
import os

from ntrugen import ntru_gen
from falcon import Falcon

# ── Load RDTSC shim ──────────────────────────────────────────────────────────
_here = os.path.dirname(os.path.abspath(__file__))
_lib  = ctypes.CDLL(os.path.join(_here, "rdtsc.so"))
_lib.rdtsc_start.restype = ctypes.c_uint64
_lib.rdtsc_stop.restype  = ctypes.c_uint64

def tsc_start():
    return _lib.rdtsc_start()

def tsc_stop():
    return _lib.rdtsc_stop()

# ── Benchmark parameters ─────────────────────────────────────────────────────
N_KEYGEN = 5
N_SIGN   = 50
N_VER    = 100


def bench(n):
    falcon = Falcon(n)
    msg    = b"benchmark message"

    # ── KeyGen (full: ntru_gen + LDL tree) ───────────────────────────────────
    kg_cycles = []
    last_fgFG = None
    for _ in range(N_KEYGEN):
        t0 = tsc_start()
        fgFG = ntru_gen(n)
        falcon.keygen(list(fgFG))
        t1 = tsc_stop()
        kg_cycles.append(t1 - t0)
        last_fgFG = fgFG

    sk, vk = falcon.keygen(list(last_fgFG))

    # ── Sign ──────────────────────────────────────────────────────────────────
    for _ in range(3):                      # warm-up
        falcon.sign(sk, msg)

    sign_cycles = []
    sigs = []
    for _ in range(N_SIGN):
        t0 = tsc_start()
        sig = falcon.sign(sk, msg)
        t1 = tsc_stop()
        sign_cycles.append(t1 - t0)
        sigs.append(sig)

    # ── Verify ────────────────────────────────────────────────────────────────
    for _ in range(5):
        falcon.verify(vk, msg, sigs[0])

    ver_cycles = []
    pool = (sigs * (N_VER // N_SIGN + 1))[:N_VER]
    for sig in pool:
        t0 = tsc_start()
        assert falcon.verify(vk, msg, sig)
        t1 = tsc_stop()
        ver_cycles.append(t1 - t0)

    return {
        "keygen_med": int(statistics.median(kg_cycles)),
        "keygen_min": min(kg_cycles),
        "sign_med":   int(statistics.median(sign_cycles)),
        "sign_min":   min(sign_cycles),
        "verify_med": int(statistics.median(ver_cycles)),
        "verify_min": min(ver_cycles),
        "siglen":     len(sigs[0]),
    }


def print_table(all_results):
    ops = [
        ("KeyGen", "keygen_med", "keygen_min"),
        ("Sign",   "sign_med",   "sign_min"),
        ("Verify", "verify_med", "verify_min"),
    ]
    ns = sorted(all_results.keys())
    cw = 22

    # Header
    header = f"{'Operation':<10}" + "".join(
        f"{'Falcon-'+str(n)+' (med)':>{cw}}{'Falcon-'+str(n)+' (min)':>{cw}}"
        for n in ns
    )
    sep = "-" * len(header)
    print(header)
    print(sep)

    for label, med_key, min_key in ops:
        row = f"{label:<10}"
        for n in ns:
            row += f"{all_results[n][med_key]:>{cw},}{all_results[n][min_key]:>{cw},}"
        print(row)

    print(sep)
    for n in ns:
        r = all_results[n]
        print(f"  Falcon-{n}  sig={r['siglen']}B  "
              f"keygen×{N_KEYGEN}  sign×{N_SIGN}  verify×{N_VER}")


if __name__ == "__main__":
    print("Platform : Intel i7-9750H @ 2.60 GHz (invariant TSC)")
    print("Metric   : RDTSC cycles  (median and min across trials)\n")

    ns = [int(sys.argv[1])] if len(sys.argv) > 1 else [512, 1024]
    all_results = {}
    for n in ns:
        print(f"Benchmarking Falcon-{n} ...", flush=True)
        all_results[n] = bench(n)

    print()
    print_table(all_results)
    print()
