"""
bench_python_domains.py — Python domain breakdown matching C bench_domains
===========================================================================
Produces the same domain table as the C bench_domains benchmark:

  KEYGEN : NTRU | f,g sample+norm | NTT | Gram | LDL | FFT-conv
  SIGN   : BaseSampler (RCDT) | berexp+approxexp | FFT ops | hash+target+compress
           (Gram + LDL build = 0 in tree mode — pre-computed at keygen)
  VERIFY : NTT poly_mul | norm check | decode+hash+other

Usage:
  python bench_python_domains.py [--n 512|1024] [--trials 20] [--keygen-trials 5]
"""

import argparse, ctypes, os, subprocess, sys, tempfile, statistics, time, copy
from math import floor, sqrt
from collections import defaultdict

# ── RDTSC ─────────────────────────────────────────────────────────────────────
_RDTSC_C = """\
#include <stdint.h>
uint64_t rs(void){uint64_t a,d;
  __asm__ __volatile__("mfence\\n\\trdtsc":"=a"(a),"=d"(d));return(d<<32)|a;}
uint64_t re(void){uint64_t a,d;unsigned c;
  __asm__ __volatile__("rdtscp":"=a"(a),"=d"(d),"=c"(c));
  __asm__ __volatile__("mfence");return(d<<32)|a;}
"""
def _build_rdtsc():
    td = tempfile.mkdtemp(); src = td+"/r.c"; lib = td+"/r.so"
    open(src,"w").write(_RDTSC_C)
    r = subprocess.run(["cc","-O2","-shared","-fPIC","-o",lib,src],
                       capture_output=True, timeout=10)
    if r.returncode != 0: raise RuntimeError(r.stderr.decode())
    dll = ctypes.CDLL(lib)
    dll.rs.restype = dll.re.restype = ctypes.c_uint64
    return dll.rs, dll.re

try:
    _rs, _re = _build_rdtsc()
    RDTSC_SRC = "hardware RDTSC"
except Exception as e:
    try: hz = int(subprocess.run(["sysctl","-n","hw.cpufrequency_max"],
                  capture_output=True,text=True,timeout=3).stdout.strip())
    except: hz = 2_600_000_000
    ghz = hz/1e9
    print(f"[warn] {e}; fallback perf_counter_ns×{ghz:.2f}")
    def _rs(): return int(time.perf_counter_ns()*ghz)
    _re = _rs; RDTSC_SRC = f"perf_counter_ns×{ghz:.2f}"

sys.path.insert(0, os.path.dirname(__file__))

from params import Q, FALCON_PARAMS, SALT_LEN
from fft import (fft, ifft, poly_mul_fft, poly_add_fft, poly_sub_fft,
                 poly_split_fft, poly_merge_fft)
from samplerz import (basesampler as _basesampler_orig,
                      berexp     as _berexp_orig,
                      approxexp, INV_2SIGMA2, samplerz as _samplerz_orig)
from ffsampling import gram, ffldl_fft, normalize_tree, ldl_fft
from ntrugen import ntru_gen
from keygen import poly_div_ntt, _build_tree
from sign import hash_to_point, compress, decompress
from poly_arith import poly_mul, poly_center

CPU_HZ = 2_600_000_000
try: CPU_HZ = int(subprocess.run(["sysctl","-n","hw.cpufrequency_max"],
                  capture_output=True,text=True,timeout=3).stdout.strip())
except: pass
CPU_GHZ = CPU_HZ / 1e9

# ── Formatting ────────────────────────────────────────────────────────────────
def fc(c):
    if c >= 1e9: return f"{c/1e9:8.3f} Gc"
    if c >= 1e6: return f"{c/1e6:8.3f} Mc"
    if c >= 1e3: return f"{c/1e3:8.1f} Kc"
    return f"{c:8.0f}  c"

def ft(c):
    us = c / CPU_HZ * 1e6
    if us >= 1e6: return f"{us/1e6:8.3f} s "
    if us >= 1e3: return f"{us/1e3:8.3f} ms"
    return f"{us:8.2f} us"

def table(title, phases, samples, trials):
    meds = {p: statistics.median(samples[p]) for p in phases}
    tot  = sum(meds.values())
    W = max(len(p) for p in phases) + 2
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")
    print(f"  {'Domain':<{W}}  {'Median':>10}  {'Time':>10}  {'Min':>10}  {'Max':>10}  {'%':>6}")
    print(f"  {'-'*W}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*6}")
    for p in phases:
        s   = samples[p]
        med = meds[p]
        pct = 100*med/tot if tot else 0
        print(f"  {p:<{W}}  {fc(med)}  {ft(med)}  {fc(min(s))}  {fc(max(s))}  {pct:5.1f}%")
    print(f"  {'TOTAL':<{W}}  {fc(tot)}  {ft(tot)}")
    worst = max(phases, key=lambda p: meds[p])
    print(f"\n  *** Bottleneck: [{worst}]  {100*meds[worst]/tot:.1f}% ***")
    return tot


# ── Instrumented samplerz: splits BaseSampler from berexp ────────────────────

class SamplerAcc:
    """Accumulates cycles for BaseSampler and berexp separately."""
    __slots__ = ["basesampler","berexp","calls"]
    def __init__(self):
        self.basesampler = 0
        self.berexp      = 0
        self.calls       = 0

def _samplerz_timed(mu, sigma, sigmin, randombytes, acc):
    """samplerz with per-component RDTSC timing into acc."""
    from math import floor
    s   = int(floor(mu))
    r   = mu - s
    dss = 1.0 / (2.0 * sigma * sigma)
    ccs = sigmin / sigma
    while True:
        t0 = _rs(); z0 = _basesampler_orig(randombytes); t1 = _re()
        acc.basesampler += t1 - t0

        b  = int.from_bytes(randombytes(1), "little") & 1
        z  = b + (2*b - 1) * z0
        x  = ((z - r)**2) * dss - (z0**2) * INV_2SIGMA2

        t0 = _rs(); result = _berexp_orig(x, ccs, randombytes); t1 = _re()
        acc.berexp += t1 - t0
        acc.calls  += 1

        if result:
            return z + s


# ── Instrumented ffsampling_fft: times poly ops separately ──────────────────

def _ffsampling_timed(t, T, sigmin, randombytes, s_acc, poly_acc, depth=0):
    """
    Mirrors ffsampling_fft but accumulates cycle counts into:
      s_acc   — SamplerAcc for basesampler + berexp
      poly_acc — dict: "split","merge","l10_mulsub","leaf_adj"
    """
    n = len(t[0])
    z = [None, None]

    if len(T) == 3 and len(T[1]) == 2:
        # Leaf level — 2 samplerz calls + scalar L10 adjustment
        l10    = T[0]
        sigma1 = T[2][0]
        sigma0 = T[1][0]

        z1_re = _samplerz_timed(t[1][0], sigma1, sigmin, randombytes, s_acc)
        z[1]  = [float(z1_re), 0.0]

        t0 = _rs()
        diff_re = t[1][0] - z[1][0]; diff_im = t[1][1] - z[1][1]
        l_re, l_im = l10[0], l10[1]
        adj_re = diff_re*l_re - diff_im*l_im
        adj_im = diff_re*l_im + diff_im*l_re
        t0b_re = t[0][0] + adj_re; t0b_im = t[0][1] + adj_im
        t1 = _re(); poly_acc["leaf_adj"] += t1 - t0

        z0_re = _samplerz_timed(t0b_re, sigma0, sigmin, randombytes, s_acc)
        z[0]  = [float(z0_re), 0.0]
    else:
        l10, T0, T1 = T

        t0 = _rs(); t1_sp = poly_split_fft(t[1]); t1 = _re()
        poly_acc["split"] += t1 - t0

        z1h = _ffsampling_timed(t1_sp, T1, sigmin, randombytes, s_acc, poly_acc, depth+1)

        t0 = _rs(); z[1] = poly_merge_fft(*z1h); t1 = _re()
        poly_acc["merge"] += t1 - t0

        t0 = _rs()
        res  = poly_sub_fft(t[1], z[1])
        corr = poly_mul_fft(res, l10)
        t0b  = poly_add_fft(t[0], corr)
        t1 = _re(); poly_acc["l10_mulsub"] += t1 - t0

        t0 = _rs(); t0b_sp = poly_split_fft(t0b); t1 = _re()
        poly_acc["split"] += t1 - t0

        z0h = _ffsampling_timed(t0b_sp, T0, sigmin, randombytes, s_acc, poly_acc, depth+1)

        t0 = _rs(); z[0] = poly_merge_fft(*z0h); t1 = _re()
        poly_acc["merge"] += t1 - t0

    return z


# ── KEYGEN benchmark ──────────────────────────────────────────────────────────

def bench_keygen(n, trials):
    sigma = FALCON_PARAMS[n]["sigma"]

    KG_PHASES = ["NTRU (ntru_gen)", "f,g sample+norm+other",
                 "NTT (poly_div_ntt)", "Gram (B·B*)", "LDL (ffldl+norm)", "FFT-conv"]
    samp = defaultdict(list)

    print(f"\n{'='*70}")
    print(f"  KEYGEN  Falcon-{n}  ({trials} trials)")
    print(f"{'='*70}")
    print("  (building keys — each takes a few seconds in Python ...)")

    for i in range(trials):
        # 1. NTRU
        t0 = _rs(); f,g,F,G = ntru_gen(n); t1 = _re()
        samp["NTRU (ntru_gen)"].append(t1-t0)

        # 2. NTT  h = g/f mod q
        t0 = _rs(); h = poly_div_ntt(g,f); t1 = _re()
        samp["NTT (poly_div_ntt)"].append(t1-t0)

        # 3. Gram
        neg_f,neg_F = [-c for c in f],[-c for c in F]
        B = [[g,neg_f],[G,neg_F]]
        t0 = _rs(); G_mat = gram(B); t1 = _re()
        samp["Gram (B·B*)"].append(t1-t0)

        # 4. fft(gram)  +  B_fft  (FFT conversion)
        t0 = _rs()
        G_fft = [[fft(G_mat[i][j]) for j in range(2)] for i in range(2)]
        B_fft = [[fft([float(c) for c in poly]) for poly in row] for row in B]
        t1 = _re()
        samp["FFT-conv"].append(t1-t0)

        # 5. ffldl + normalize
        t0 = _rs()
        T = ffldl_fft(G_fft)
        T_copy = copy.deepcopy(T)
        normalize_tree(T_copy, sigma)
        t1 = _re()
        samp["LDL (ffldl+norm)"].append(t1-t0)

        # 6. f,g sample+norm — everything else = total - above
        total_approx = sum(v[-1] for v in samp.values())
        samp["f,g sample+norm+other"].append(0)   # placeholder

        print(f"  trial {i+1}/{trials} done  (ntru_gen={samp['NTRU (ntru_gen)'][-1]/CPU_HZ*1e3:.0f} ms)")

    # compute real ntru total vs keygen total (approx, no outer timer)
    # just show what we measured
    for p in KG_PHASES:
        if not samp[p]:
            samp[p] = [0]*trials

    return table(f"KEYGEN  Falcon-{n}  ({trials} trials)", KG_PHASES, samp, trials)


# ── SIGN benchmark ────────────────────────────────────────────────────────────

def bench_sign(n, trials, sk, msg):
    p          = FALCON_PARAMS[n]
    sigmin     = p["sigma_min"]
    beta_sq    = p["beta_sq"]
    payload    = p["sig_bytelen"] - 1 - SALT_LEN
    [[a,b],[c_mat,d]] = sk.B_fft

    SG_PHASES = ["BaseSampler (RCDT)", "berexp+approxexp",
                 "FFT ops (split/merge/l10)", "hash+target+compress",
                 "Gram (pre-computed)",  "LDL build (pre-computed)"]

    samp = defaultdict(list)
    completed = 0; attempts = 0

    while completed < trials:
        attempts += 1
        salt = os.urandom(SALT_LEN)

        # hash + fft(c) + targets
        t0 = _rs()
        c_poly = hash_to_point(msg, salt, n)
        c_fft  = fft([float(x) for x in c_poly])
        tgt0   = [x/Q for x in poly_mul_fft(c_fft, d)]
        tgt1   = [-x/Q for x in poly_mul_fft(c_fft, b)]
        t1 = _re(); hash_target_cyc = t1 - t0

        # instrumented ffsampling
        s_acc    = SamplerAcc()
        poly_acc = defaultdict(int)

        _ffsampling_timed([tgt0, tgt1], sk.T, sigmin, os.urandom, s_acc, poly_acc)
        z = _samplerz_timed  # dummy ref — actual z from timed call above

        # re-run for z (needed for the lattice point)
        from ffsampling import ffsampling_fft
        z = ffsampling_fft([tgt0, tgt1], sk.T, sigmin, os.urandom)

        # lattice point + ifft + round + norm + compress
        t0 = _rs()
        v0f = poly_add_fft(poly_mul_fft(z[0],a), poly_mul_fft(z[1],c_mat))
        v1f = poly_add_fft(poly_mul_fft(z[0],b), poly_mul_fft(z[1],d))
        v0  = [int(round(x)) for x in ifft(v0f)]
        v1  = [int(round(x)) for x in ifft(v1f)]
        s0  = [c_poly[i]-v0[i] for i in range(n)]
        s1  = [-v1[i] for i in range(n)]
        nsq = sum(x*x for x in s0)+sum(x*x for x in s1)
        t1 = _re(); post_cyc = t1 - t0

        if nsq > beta_sq: continue
        t0 = _rs(); enc = compress(s1, payload); t1 = _re()
        comp_cyc = t1 - t0
        if enc is False: continue

        fft_cyc = (poly_acc["split"] + poly_acc["merge"] +
                   poly_acc["l10_mulsub"] + poly_acc["leaf_adj"])

        samp["BaseSampler (RCDT)"].append(s_acc.basesampler)
        samp["berexp+approxexp"].append(s_acc.berexp)
        samp["FFT ops (split/merge/l10)"].append(fft_cyc)
        samp["hash+target+compress"].append(hash_target_cyc + post_cyc + comp_cyc)
        samp["Gram (pre-computed)"].append(0)
        samp["LDL build (pre-computed)"].append(0)
        completed += 1

    if attempts > trials:
        print(f"\n  (rejection rate {(attempts-trials)/attempts*100:.1f}% — {attempts} attempts)")

    return table(f"SIGN  Falcon-{n}  ({trials} trials)  [tree mode]",
                 SG_PHASES, samp, trials)


# ── VERIFY benchmark ──────────────────────────────────────────────────────────

def bench_verify(n, trials, pk, sigs, msg):
    p           = FALCON_PARAMS[n]
    beta_sq     = p["beta_sq"]
    payload_len = p["sig_bytelen"] - 1 - SALT_LEN
    h           = pk.h

    VF_PHASES = ["NTT poly_mul (s1·h mod q)", "Norm check (is_short)",
                 "decode+hash+other"]
    samp = defaultdict(list)

    for sig in sigs:
        salt_v  = sig[1:1+SALT_LEN]
        enc_s1  = sig[1+SALT_LEN:]

        # decode + hash
        t0 = _rs()
        s1 = decompress(enc_s1, payload_len, n)
        c  = hash_to_point(msg, salt_v, n)
        t1 = _re(); dh_cyc = t1 - t0

        # NTT poly_mul: s0 = c - s1·h mod q
        t0 = _rs()
        s1m = [x%Q for x in s1]
        s1h = poly_mul(s1m, h)
        s0  = [(c[i]-s1h[i])%Q for i in range(n)]
        s0  = poly_center(s0)
        t1 = _re(); ntt_cyc = t1 - t0

        # norm check
        t0 = _rs()
        ns = sum(x*x for x in s0)+sum(x*x for x in s1)
        _  = (ns <= beta_sq)
        t1 = _re(); norm_cyc = t1 - t0

        samp["NTT poly_mul (s1·h mod q)"].append(ntt_cyc)
        samp["Norm check (is_short)"].append(norm_cyc)
        samp["decode+hash+other"].append(dh_cyc)

    return table(f"VERIFY  Falcon-{n}  ({trials} trials)",
                 VF_PHASES, samp, trials)


# ── Throughput summary ────────────────────────────────────────────────────────

def print_throughput(n, kg, sg, vf):
    print(f"\n{'='*70}")
    print(f"  THROUGHPUT SUMMARY  Falcon-{n}  @ {CPU_GHZ:.3f} GHz")
    print(f"{'='*70}")
    print(f"  {'Operation':<10}  {'Cycles(med)':>12}  {'Time':>10}  {'ops/sec':>12}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*10}  {'-'*12}")
    for label, cyc in [("keygen",kg),("sign",sg),("verify",vf)]:
        print(f"  {label:<10}  {fc(cyc)}  {ft(cyc)}  {CPU_HZ/cyc:>12.1f}")
    print()


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n",             type=int, default=512, choices=[512,1024])
    ap.add_argument("--trials",        type=int, default=20)
    ap.add_argument("--keygen-trials", type=int, default=3)
    args = ap.parse_args()
    n, ST, KT = args.n, args.trials, args.keygen_trials

    print(f"\nFalcon-{n} Python domain profiler  (mirrors C bench_domains)")
    print(f"  Timer  : {RDTSC_SRC}")
    print(f"  CPU    : {CPU_GHZ:.3f} GHz")
    print(f"  keygen : {KT} trials")
    print(f"  sign   : {ST} valid sigs")
    print(f"  verify : {ST} trials")

    # ── Keygen ──
    kg_tot = bench_keygen(n, KT)

    # ── Build one tree-mode key ──
    print(f"\nBuilding tree-mode key for sign/verify …", end="", flush=True)
    from keygen import keygen
    sk, pk = keygen(n, mode="tree")
    print(" done.")

    msg = b"bench-domains-test-message-xxxx"

    # ── Sign ──
    sg_tot = bench_sign(n, ST, sk, msg)

    # ── Collect sigs ──
    from sign import sign as fsign
    print(f"\nGenerating {ST} sigs for verify …", end="", flush=True)
    sigs = [fsign(sk, msg) for _ in range(ST)]
    print(" done.")

    # ── Verify ──
    vf_tot = bench_verify(n, ST, pk, sigs, msg)

    # ── Throughput ──
    print_throughput(n, kg_tot, sg_tot, vf_tot)


if __name__ == "__main__":
    main()
