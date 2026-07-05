"""
bench_rdtsc.py — RDTSC-based Falcon bottleneck profiler
=========================================================

Measures cycle counts for every major phase of Falcon-512 and Falcon-1024:

  keygen  : ntru_gen │ poly_div_ntt │ gram │ fft_gram │ ffldl │ normalize │ B_fft
  sign    : hash_to_point │ fft_c │ targets │ ffsampling │ v_fft │ ifft+round │ compress
  verify  : decompress │ hash_to_point │ poly_mul │ norm_check

Results: median cycles per trial, min/max, µs, and % of phase total.

Usage:
  python bench_rdtsc.py [--n 512|1024] [--trials 20] [--keygen-trials 5]
"""

import argparse
import ctypes
import os
import subprocess
import sys
import tempfile
import time
import statistics
from math import sqrt

# ── RDTSC ─────────────────────────────────────────────────────────────────────

_RDTSC_C = """\
#include <stdint.h>
uint64_t rdtsc(void) {
    uint64_t lo, hi;
    __asm__ __volatile__(
        "mfence\\n\\t"
        "rdtsc"
        : "=a"(lo), "=d"(hi)
    );
    return (hi << 32) | lo;
}
uint64_t rdtscp(void) {
    uint64_t lo, hi;
    unsigned int aux;
    __asm__ __volatile__(
        "rdtscp"
        : "=a"(lo), "=d"(hi), "=c"(aux)
    );
    __asm__ __volatile__("mfence");
    return (hi << 32) | lo;
}
"""

def _build_rdtsc():
    """Compile a tiny RDTSC shared lib; return (start_fn, stop_fn, source_label)."""
    td = tempfile.mkdtemp(prefix="falcon_rdtsc_")
    src = os.path.join(td, "rdtsc.c")
    lib = os.path.join(td, "rdtsc.so")
    with open(src, "w") as f:
        f.write(_RDTSC_C)
    try:
        r = subprocess.run(
            ["cc", "-O2", "-shared", "-fPIC", "-o", lib, src],
            capture_output=True, timeout=10
        )
        if r.returncode != 0:
            raise RuntimeError(r.stderr.decode())
        dll = ctypes.CDLL(lib)
        dll.rdtsc.restype  = ctypes.c_uint64
        dll.rdtscp.restype = ctypes.c_uint64
        # sanity-check: two consecutive calls should differ
        a, b = dll.rdtsc(), dll.rdtsc()
        if b <= a:
            raise RuntimeError("RDTSC not advancing")
        return dll.rdtsc, dll.rdtscp, "hardware RDTSC (mfence·rdtsc / rdtscp·mfence)"
    except Exception as e:
        # ── Fallback: perf_counter_ns → cycles ──────────────────────────────
        try:
            res = subprocess.run(
                ["sysctl", "-n", "hw.cpufrequency_max"],
                capture_output=True, text=True, timeout=5
            )
            hz = int(res.stdout.strip())
        except Exception:
            hz = 2_600_000_000  # i7-9750H default
        ghz = hz / 1e9
        print(f"[RDTSC] compile failed ({e}); falling back to perf_counter_ns × {ghz:.2f} GHz")
        def _fake_rdtsc():
            return int(time.perf_counter_ns() * ghz)
        return _fake_rdtsc, _fake_rdtsc, f"perf_counter_ns × {ghz:.2f} GHz"


_rdtsc_start, _rdtsc_stop, RDTSC_SOURCE = _build_rdtsc()

def rdtsc_start() -> int:
    """Use before the timed region (serialising RDTSC)."""
    return _rdtsc_start()

def rdtsc_stop() -> int:
    """Use after the timed region (serialising RDTSCP)."""
    return _rdtsc_stop()

# ── Project imports ────────────────────────────────────────────────────────────

sys.path.insert(0, os.path.dirname(__file__))

from params import Q, FALCON_PARAMS, SALT_LEN
from fft import fft, ifft, poly_mul_fft, poly_add_fft
from ntrugen import ntru_gen
from ffsampling import gram, ffldl_fft, normalize_tree, ffsampling_fft
from keygen import poly_div_ntt, _build_tree, SecretKey, PublicKey
from sign import hash_to_point, compress, decompress
from poly_arith import poly_mul, poly_center

# ── Helpers ────────────────────────────────────────────────────────────────────

CPU_HZ = 2_600_000_000  # Intel i7-9750H — update if running on different HW
try:
    _r = subprocess.run(["sysctl", "-n", "hw.cpufrequency_max"],
                        capture_output=True, text=True, timeout=3)
    CPU_HZ = int(_r.stdout.strip())
except Exception:
    pass

CPU_GHZ = CPU_HZ / 1e9


def cycles_to_us(c: float) -> float:
    return c / CPU_HZ * 1e6


def fmt_cyc(c: float) -> str:
    if c >= 1e9:  return f"{c/1e9:8.3f} Gc"
    if c >= 1e6:  return f"{c/1e6:8.3f} Mc"
    if c >= 1e3:  return f"{c/1e3:8.1f} Kc"
    return f"{c:8.0f}  c"


def fmt_us(c: float) -> str:
    us = cycles_to_us(c)
    if us >= 1e6:  return f"{us/1e6:8.2f} s "
    if us >= 1e3:  return f"{us/1e3:8.2f} ms"
    return f"{us:8.1f} µs"


def median_cycles(samples):
    return statistics.median(samples)


def print_phase_table(title: str, phases: list, phase_samples: dict):
    """
    phases      : list of phase-name strings
    phase_samples: {phase: [cycle_count, ...]}
    """
    medians = {p: median_cycles(phase_samples[p]) for p in phases}
    total   = sum(medians.values())

    W = max(len(p) for p in phases) + 2
    print(f"\n  {'Phase':<{W}}  {'Median':>12}  {'µs':>10}  {'Min':>12}  {'Max':>12}  {'%':>6}")
    print(f"  {'-'*W}  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*6}")
    for p in phases:
        s   = phase_samples[p]
        med = medians[p]
        mn  = min(s)
        mx  = max(s)
        pct = 100.0 * med / total if total else 0
        print(f"  {p:<{W}}  {fmt_cyc(med)}  {fmt_us(med)}  {fmt_cyc(mn)}  {fmt_cyc(mx)}  {pct:5.1f}%")
    print(f"  {'TOTAL':<{W}}  {fmt_cyc(total)}  {fmt_us(total)}")
    print()

    # Bottleneck callout
    worst = max(phases, key=lambda p: medians[p])
    pct   = 100.0 * medians[worst] / total
    print(f"  *** Bottleneck: [{worst}]  {pct:.1f}% of {title} ***")


# ── Keygen phases ──────────────────────────────────────────────────────────────

def bench_keygen(n: int, trials: int):
    print(f"\n{'='*70}")
    print(f"  KEYGEN  Falcon-{n}   ({trials} trials)")
    print(f"{'='*70}")

    PHASES = [
        "ntru_gen",
        "poly_div_ntt",
        "gram",
        "fft_gram",
        "ffldl",
        "normalize",
        "B_fft",
    ]
    samples = {p: [] for p in PHASES}
    sigma = FALCON_PARAMS[n]["sigma"]

    for _ in range(trials):
        # 1. ntru_gen
        t0 = rdtsc_start(); f, g, F, G = ntru_gen(n); t1 = rdtsc_stop()
        samples["ntru_gen"].append(t1 - t0)

        # 2. poly_div_ntt  (public key h = g/f mod q)
        t0 = rdtsc_start(); h = poly_div_ntt(g, f); t1 = rdtsc_stop()
        samples["poly_div_ntt"].append(t1 - t0)

        # 3. gram
        neg_f = [-c for c in f]; neg_F = [-c for c in F]
        B = [[g, neg_f], [G, neg_F]]
        t0 = rdtsc_start(); G_mat = gram(B); t1 = rdtsc_stop()
        samples["gram"].append(t1 - t0)

        # 4. fft of gram matrix
        t0 = rdtsc_start()
        G_fft = [[fft(G_mat[i][j]) for j in range(2)] for i in range(2)]
        t1 = rdtsc_stop()
        samples["fft_gram"].append(t1 - t0)

        # 5. ffldl tree
        t0 = rdtsc_start(); T = ffldl_fft(G_fft); t1 = rdtsc_stop()
        samples["ffldl"].append(t1 - t0)

        # 6. normalize
        import copy
        T_copy = copy.deepcopy(T)
        t0 = rdtsc_start(); normalize_tree(T_copy, sigma); t1 = rdtsc_stop()
        samples["normalize"].append(t1 - t0)

        # 7. FFT basis B_fft
        t0 = rdtsc_start()
        B_fft = [[fft([float(c) for c in poly]) for poly in row] for row in B]
        t1 = rdtsc_stop()
        samples["B_fft"].append(t1 - t0)

    print_phase_table(f"Falcon-{n} keygen", PHASES, samples)
    return samples


# ── Sign phases ────────────────────────────────────────────────────────────────

def bench_sign(n: int, trials: int, sk, msg: bytes):
    print(f"\n{'='*70}")
    print(f"  SIGN  Falcon-{n}   ({trials} trials)")
    print(f"{'='*70}")

    PHASES = [
        "hash_to_point",
        "fft_c",
        "compute_targets",
        "ffsampling",
        "v_fft_muls",
        "ifft_round",
        "norm_check",
        "compress",
    ]
    samples = {p: [] for p in PHASES}

    p      = FALCON_PARAMS[n]
    sigmin = p["sigma_min"]
    beta_sq     = p["beta_sq"]
    payload_len = p["sig_bytelen"] - 1 - SALT_LEN

    [[a, b], [c_mat, d]] = sk.B_fft

    attempts = 0
    completed = 0
    while completed < trials:
        attempts += 1
        salt = os.urandom(SALT_LEN)

        # 1. hash_to_point
        t0 = rdtsc_start()
        c_poly = hash_to_point(msg, salt, n)
        t1 = rdtsc_stop()
        _hash_cyc = t1 - t0

        # 2. fft(c)
        t0 = rdtsc_start()
        c_fft = fft([float(x) for x in c_poly])
        t1 = rdtsc_stop()
        _fft_c_cyc = t1 - t0

        # 3. compute targets t0, t1
        t0 = rdtsc_start()
        tgt0 = [x / Q for x in poly_mul_fft(c_fft, d)]
        tgt1 = [-x / Q for x in poly_mul_fft(c_fft, b)]
        t1 = rdtsc_stop()
        _targets_cyc = t1 - t0

        # 4. ffsampling
        t0 = rdtsc_start()
        z = ffsampling_fft([tgt0, tgt1], sk.T, sigmin, os.urandom)
        t1 = rdtsc_stop()
        _ffsamp_cyc = t1 - t0

        # 5. v = z · B  (4 poly_mul_fft + 2 poly_add_fft)
        t0 = rdtsc_start()
        v0_fft = poly_add_fft(poly_mul_fft(z[0], a), poly_mul_fft(z[1], c_mat))
        v1_fft = poly_add_fft(poly_mul_fft(z[0], b), poly_mul_fft(z[1], d))
        t1 = rdtsc_stop()
        _vmul_cyc = t1 - t0

        # 6. ifft + round  (back to integer coefficients)
        t0 = rdtsc_start()
        v0 = [int(round(x)) for x in ifft(v0_fft)]
        v1 = [int(round(x)) for x in ifft(v1_fft)]
        t1 = rdtsc_stop()
        _ifft_cyc = t1 - t0

        # (s0, s1) = (c − v0, −v1)
        s0 = [c_poly[i] - v0[i] for i in range(n)]
        s1 = [-v1[i] for i in range(n)]

        # 7. norm check
        t0 = rdtsc_start()
        norm_sq = sum(x * x for x in s0) + sum(x * x for x in s1)
        norm_ok = (norm_sq <= beta_sq)
        t1 = rdtsc_stop()
        _norm_cyc = t1 - t0

        if not norm_ok:
            continue

        # 8. compress
        t0 = rdtsc_start()
        enc = compress(s1, payload_len)
        t1 = rdtsc_stop()
        _comp_cyc = t1 - t0

        if enc is False:
            continue

        # All checks passed — record this trial
        samples["hash_to_point"].append(_hash_cyc)
        samples["fft_c"].append(_fft_c_cyc)
        samples["compute_targets"].append(_targets_cyc)
        samples["ffsampling"].append(_ffsamp_cyc)
        samples["v_fft_muls"].append(_vmul_cyc)
        samples["ifft_round"].append(_ifft_cyc)
        samples["norm_check"].append(_norm_cyc)
        samples["compress"].append(_comp_cyc)
        completed += 1

    if attempts > trials:
        print(f"  (sign rejection rate: {(attempts-trials)/attempts*100:.1f}%  — {attempts} attempts for {trials} valid sigs)")

    print_phase_table(f"Falcon-{n} sign", PHASES, samples)
    return samples


# ── Verify phases ──────────────────────────────────────────────────────────────

def bench_verify(n: int, trials: int, pk, sigs: list, msg: bytes):
    print(f"\n{'='*70}")
    print(f"  VERIFY  Falcon-{n}   ({trials} trials)")
    print(f"{'='*70}")

    PHASES = [
        "decompress",
        "hash_to_point",
        "poly_mul_ntt",
        "norm_check",
    ]
    samples = {p: [] for p in PHASES}

    p           = FALCON_PARAMS[n]
    beta_sq     = p["beta_sq"]
    payload_len = p["sig_bytelen"] - 1 - SALT_LEN
    h           = pk.h

    for i in range(trials):
        sig = sigs[i % len(sigs)]
        salt   = sig[1 : 1 + SALT_LEN]
        enc_s1 = sig[1 + SALT_LEN :]

        # 1. decompress
        t0 = rdtsc_start()
        s1 = decompress(enc_s1, payload_len, n)
        t1 = rdtsc_stop()
        samples["decompress"].append(t1 - t0)

        # 2. hash_to_point
        t0 = rdtsc_start()
        c = hash_to_point(msg, salt, n)
        t1 = rdtsc_stop()
        samples["hash_to_point"].append(t1 - t0)

        # 3. NTT poly_mul: s0 = c − s1·h mod q
        t0 = rdtsc_start()
        s1_mod = [x % Q for x in s1]
        s1h    = poly_mul(s1_mod, h)
        s0     = [(c[i] - s1h[i]) % Q for i in range(n)]
        s0     = poly_center(s0)
        t1 = rdtsc_stop()
        samples["poly_mul_ntt"].append(t1 - t0)

        # 4. norm check
        t0 = rdtsc_start()
        norm_sq = sum(x * x for x in s0) + sum(x * x for x in s1)
        _ = norm_sq <= beta_sq
        t1 = rdtsc_stop()
        samples["norm_check"].append(t1 - t0)

    print_phase_table(f"Falcon-{n} verify", PHASES, samples)
    return samples


# ── Throughput summary ─────────────────────────────────────────────────────────

def throughput_summary(n: int, kg_s: dict, sg_s: dict, vf_s: dict):
    KG_PHASES = list(kg_s.keys())
    SG_PHASES = list(sg_s.keys())
    VF_PHASES = list(vf_s.keys())

    kg_med = sum(median_cycles(kg_s[p]) for p in KG_PHASES)
    sg_med = sum(median_cycles(sg_s[p]) for p in SG_PHASES)
    vf_med = sum(median_cycles(vf_s[p]) for p in VF_PHASES)

    def ops_per_sec(c):
        return CPU_HZ / c

    print(f"\n{'='*70}")
    print(f"  THROUGHPUT SUMMARY   Falcon-{n}   @ {CPU_GHZ:.2f} GHz")
    print(f"{'='*70}")
    print(f"  {'Operation':<12}  {'Cycles (median)':>16}  {'Time':>10}  {'ops/s':>12}")
    print(f"  {'-'*12}  {'-'*16}  {'-'*10}  {'-'*12}")
    for label, cyc in [("keygen", kg_med), ("sign", sg_med), ("verify", vf_med)]:
        print(f"  {label:<12}  {fmt_cyc(cyc)}  {fmt_us(cyc)}  {ops_per_sec(cyc):>10.1f}/s")
    print()


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="RDTSC-based Falcon bottleneck profiler")
    ap.add_argument("--n",             type=int, default=512,  choices=[512, 1024])
    ap.add_argument("--trials",        type=int, default=20,   help="sign/verify trial count")
    ap.add_argument("--keygen-trials", type=int, default=5,    help="keygen trial count (slow)")
    args = ap.parse_args()

    n  = args.n
    ST = args.trials
    KT = args.keygen_trials

    print(f"\nFalcon-{n} RDTSC bottleneck profiler")
    print(f"  Timer source : {RDTSC_SOURCE}")
    print(f"  CPU          : {CPU_GHZ:.3f} GHz  ({CPU_HZ:,} Hz)")
    print(f"  keygen trials: {KT}")
    print(f"  sign/verify  : {ST} valid signatures each")

    # ── Keygen ────────────────────────────────────────────────────────────────
    kg_samples = bench_keygen(n, KT)

    # ── Build one tree-mode key for sign/verify ────────────────────────────────
    print(f"\nBuilding tree-mode key for sign/verify benchmarks …", end="", flush=True)
    from keygen import keygen
    sk, pk = keygen(n, mode="tree")
    print(" done.")

    msg = b"benchmark message for Falcon sign/verify"

    # ── Sign ──────────────────────────────────────────────────────────────────
    sg_samples = bench_sign(n, ST, sk, msg)

    # ── Collect valid signatures for verify ───────────────────────────────────
    from sign import sign as falcon_sign
    print(f"\nGenerating {ST} signatures for verify bench …", end="", flush=True)
    sigs = [falcon_sign(sk, msg) for _ in range(ST)]
    print(" done.")

    # ── Verify ────────────────────────────────────────────────────────────────
    vf_samples = bench_verify(n, ST, pk, sigs, msg)

    # ── Throughput summary ────────────────────────────────────────────────────
    throughput_summary(n, kg_samples, sg_samples, vf_samples)


if __name__ == "__main__":
    main()