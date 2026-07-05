"""
bench_domains.py — RDTSC domain profiler for falcon.py-master 3
================================================================
Mirrors the C bench_domains output exactly:

  KEYGEN : NTRU | Gram | LDL | FFT-conv | NTT (div_zq)
  SIGN   : BaseSampler (RCDT) | berexp+approxexp | FFT ops | hash+target+compress
  VERIFY : NTT poly_mul (mul_zq) | norm check | decode+hash+other

Uses the rdtsc.so already present in this directory.

Usage:
  python bench_domains.py [--n 512|1024] [--trials 20] [--keygen-trials 3]
"""

import argparse, ctypes, os, sys, statistics, time, copy
from math import floor
from collections import defaultdict
from os import urandom

# ── RDTSC — use existing rdtsc.so in this directory ───────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))

def _load_rdtsc():
    lib_path = os.path.join(_HERE, "rdtsc.so")
    if not os.path.exists(lib_path):
        # compile it
        import subprocess, tempfile
        r = subprocess.run(
            ["cc", "-O2", "-shared", "-fPIC", "-o", lib_path,
             os.path.join(_HERE, "rdtsc.c")],
            capture_output=True, timeout=10)
        if r.returncode != 0:
            raise RuntimeError(r.stderr.decode())
    dll = ctypes.CDLL(lib_path)
    dll.rdtsc_start.restype = ctypes.c_uint64
    dll.rdtsc_stop.restype  = ctypes.c_uint64
    # sanity check
    a, b = dll.rdtsc_start(), dll.rdtsc_stop()
    if b <= a:
        raise RuntimeError("RDTSC not advancing")
    return dll.rdtsc_start, dll.rdtsc_stop

try:
    _rs, _re = _load_rdtsc()
    RDTSC_SRC = "hardware RDTSC (rdtsc.so)"
except Exception as e:
    import subprocess
    try:
        hz = int(subprocess.run(["sysctl","-n","hw.cpufrequency_max"],
                  capture_output=True, text=True, timeout=3).stdout.strip())
    except Exception:
        hz = 2_600_000_000
    ghz = hz / 1e9
    print(f"[warn] {e}; fallback perf_counter_ns×{ghz:.2f}")
    def _rs(): return int(time.perf_counter_ns() * ghz)
    _re = _rs
    RDTSC_SRC = f"perf_counter_ns×{ghz:.2f}"

# ── falcon.py-master imports ──────────────────────────────────────────────────
sys.path.insert(0, _HERE)

from common import q
from fft import fft, ifft, add_fft, sub_fft, mul_fft, split_fft, merge_fft, fft_ratio
from ntt import div_zq, mul_zq, sub_zq
from ffsampling import gram, ffldl_fft
from math import sqrt as _sqrt

def normalize_tree(tree, sigma):
    if len(tree) == 3:
        normalize_tree(tree[1], sigma)
        normalize_tree(tree[2], sigma)
    else:
        tree[0] = sigma / _sqrt(tree[0].real)
        tree[1] = 0
from ntrugen import ntru_gen
from samplerz import (basesampler as _bs_orig, berexp as _be_orig,
                      samplerz as _sz_orig, INV_2SIGMA2)
from encoding import compress, decompress
from Crypto.Hash import SHAKE256

CPU_HZ = 2_600_000_000
try:
    import subprocess
    CPU_HZ = int(subprocess.run(["sysctl","-n","hw.cpufrequency_max"],
                  capture_output=True, text=True, timeout=3).stdout.strip())
except Exception:
    pass
CPU_GHZ = CPU_HZ / 1e9

SALT_LEN = 40
HEAD_LEN = 1

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

def print_table(title, phases, samples):
    meds = {p: statistics.median(samples[p]) for p in phases}
    tot  = sum(meds.values())
    W    = max(len(p) for p in phases) + 2
    print(f"\n{'='*72}")
    print(f"  {title}")
    print(f"{'='*72}")
    print(f"  {'Domain':<{W}}  {'Median':>10}  {'Time':>10}  {'Min':>10}  {'Max':>10}  {'%':>6}")
    print(f"  {'-'*W}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*6}")
    for p in phases:
        s = samples[p]; med = meds[p]; pct = 100*med/tot if tot else 0
        print(f"  {p:<{W}}  {fc(med)}  {ft(med)}  {fc(min(s))}  {fc(max(s))}  {pct:5.1f}%")
    print(f"  {'TOTAL':<{W}}  {fc(tot)}  {ft(tot)}")
    worst = max((p for p in phases if meds[p] > 0), key=lambda p: meds[p], default=phases[0])
    print(f"\n  *** Bottleneck: [{worst}]  {100*meds[worst]/tot:.1f}% ***")
    return tot

# ── hash_to_point (replicated from Falcon class) ─────────────────────────────
def hash_to_point(message, salt, n):
    k = (1 << 16) // q
    shake = SHAKE256.new()
    shake.update(salt)
    shake.update(message)
    hashed = [0]*n; i = 0
    while i < n:
        tb = shake.read(2)
        elt = (tb[0] << 8) + tb[1]
        if elt < k * q:
            hashed[i] = elt % q; i += 1
    return hashed

# ── Instrumented samplerz ─────────────────────────────────────────────────────
class SamplerAcc:
    __slots__ = ["bs","be","calls"]
    def __init__(self): self.bs = 0; self.be = 0; self.calls = 0

def _sz_timed(mu, sigma, sigmin, randombytes, acc):
    s   = int(floor(mu)); r = mu - s
    dss = 1.0 / (2.0 * sigma * sigma); ccs = sigmin / sigma
    while True:
        t0 = _rs(); z0 = _bs_orig(randombytes); t1 = _re()
        acc.bs += t1 - t0
        b  = int.from_bytes(randombytes(1), "little") & 1
        z  = b + (2*b - 1)*z0
        x  = ((z - r)**2)*dss - z0**2 * INV_2SIGMA2
        t0 = _rs(); ok = _be_orig(x, ccs, randombytes); t1 = _re()
        acc.be += t1 - t0; acc.calls += 1
        if ok: return z + s

# ── Instrumented ffsampling_fft ───────────────────────────────────────────────
def _ffsamp_timed(t, T, sigmin, randombytes, s_acc, poly_acc):
    n = len(t[0]) * fft_ratio
    z = [0, 0]
    if n > 1:
        l10, T0, T1 = T

        t0 = _rs(); sp = split_fft(t[1]); t1 = _re()
        poly_acc["split"] += t1 - t0

        z1h = _ffsamp_timed(sp, T1, sigmin, randombytes, s_acc, poly_acc)

        t0 = _rs(); z[1] = merge_fft(z1h); t1 = _re()
        poly_acc["merge"] += t1 - t0

        t0 = _rs()
        t0b = add_fft(t[0], mul_fft(sub_fft(t[1], z[1]), l10))
        t1 = _re(); poly_acc["l10_mulsub"] += t1 - t0

        t0 = _rs(); sp0 = split_fft(t0b); t1 = _re()
        poly_acc["split"] += t1 - t0

        z0h = _ffsamp_timed(sp0, T0, sigmin, randombytes, s_acc, poly_acc)

        t0 = _rs(); z[0] = merge_fft(z0h); t1 = _re()
        poly_acc["merge"] += t1 - t0

    elif n == 1:
        z[0] = [_sz_timed(t[0][0].real, T[0], sigmin, randombytes, s_acc)]
        z[1] = [_sz_timed(t[1][0].real, T[0], sigmin, randombytes, s_acc)]
    return z

# ── KEYGEN benchmark ──────────────────────────────────────────────────────────
def bench_keygen(n, trials, sigma):
    KG = ["NTRU (ntru_gen)", "Gram (B·B*)", "FFT-conv (B0+G0)",
          "LDL (ffldl_fft+norm)", "NTT (div_zq h=g/f)",
          "f,g sample+norm+other"]
    samp = defaultdict(list)

    print(f"\n{'='*72}")
    print(f"  KEYGEN  Falcon-{n}  ({trials} trials)  [falcon.py-master 3]")
    print(f"{'='*72}")
    print("  (each trial runs ntru_gen — takes several seconds) ...")

    for i in range(trials):
        t0 = _rs(); f,g,F,G = ntru_gen(n); t1 = _re()
        samp["NTRU (ntru_gen)"].append(t1-t0)

        from fft import neg
        B0 = [[g, neg(f)], [G, neg(F)]]

        t0 = _rs(); G0 = gram(B0); t1 = _re()
        samp["Gram (B·B*)"].append(t1-t0)

        t0 = _rs()
        B0_fft = [[fft(elt) for elt in row] for row in B0]
        G0_fft = [[fft(elt) for elt in row] for row in G0]
        t1 = _re(); samp["FFT-conv (B0+G0)"].append(t1-t0)

        t0 = _rs()
        T = ffldl_fft(G0_fft)
        normalize_tree(T, sigma)
        t1 = _re(); samp["LDL (ffldl_fft+norm)"].append(t1-t0)

        t0 = _rs(); h = div_zq(g, f); t1 = _re()
        samp["NTT (div_zq h=g/f)"].append(t1-t0)

        samp["f,g sample+norm+other"].append(0)
        print(f"  trial {i+1}/{trials}  ntru_gen={samp['NTRU (ntru_gen)'][-1]/CPU_HZ*1e3:.0f} ms")

    return print_table(f"KEYGEN Falcon-{n} ({trials} trials)", KG, samp)

# ── SIGN benchmark ────────────────────────────────────────────────────────────
def bench_sign(n, trials, sk, sigma, sigmin, sig_bound, sig_bytelen):
    (f, g, F, G, B0_fft, T_fft) = sk
    [[a,b],[c,d]] = B0_fft
    payload = sig_bytelen - HEAD_LEN - SALT_LEN

    SG = ["BaseSampler (RCDT)", "berexp+approxexp",
          "FFT ops (split/merge/l10)", "hash+target+compress",
          "Gram (pre-computed)", "LDL build (pre-computed)"]
    samp = defaultdict(list)
    done = 0; attempts = 0

    while done < trials:
        attempts += 1
        salt = urandom(SALT_LEN)

        # hash + target setup
        t0 = _rs()
        hashed = hash_to_point(msg, salt, n)
        point_fft = fft(hashed)
        t0_fft = [point_fft[i]*d[i]/q for i in range(n)]
        t1_fft = [-point_fft[i]*b[i]/q for i in range(n)]
        t1 = _re(); hash_tgt = t1 - t0

        # instrumented ffsampling
        s_acc = SamplerAcc(); poly_acc = defaultdict(int)
        z_fft = _ffsamp_timed([t0_fft, t1_fft], T_fft, sigmin, urandom, s_acc, poly_acc)

        # lattice point + norm + compress
        t0 = _rs()
        v0_fft = add_fft([z_fft[0][i]*a[i] for i in range(n)],
                         [z_fft[1][i]*c[i] for i in range(n)])
        v1_fft = add_fft([z_fft[0][i]*b[i] for i in range(n)],
                         [z_fft[1][i]*d[i] for i in range(n)])

        from fft import sub, neg as fneg
        v0 = [int(round(x)) for x in ifft(v0_fft)]
        v1 = [int(round(x)) for x in ifft(v1_fft)]
        s0 = sub(hashed, v0)
        s1 = fneg(v1)
        nsq = sum(x*x for x in s0) + sum(x*x for x in s1)
        t1 = _re(); post = t1 - t0

        if nsq > sig_bound: continue
        t0 = _rs(); enc = compress(s1, payload); t1 = _re()
        if enc is False: continue
        comp = t1 - t0

        fft_cyc = poly_acc["split"] + poly_acc["merge"] + poly_acc["l10_mulsub"]
        samp["BaseSampler (RCDT)"].append(s_acc.bs)
        samp["berexp+approxexp"].append(s_acc.be)
        samp["FFT ops (split/merge/l10)"].append(fft_cyc)
        samp["hash+target+compress"].append(hash_tgt + post + comp)
        samp["Gram (pre-computed)"].append(0)
        samp["LDL build (pre-computed)"].append(0)
        done += 1

    if attempts > trials:
        print(f"\n  (rejection rate {(attempts-trials)/attempts*100:.1f}% — {attempts} attempts)")

    return print_table(f"SIGN Falcon-{n} ({trials} trials) [tree mode]", SG, samp)

# ── VERIFY benchmark ──────────────────────────────────────────────────────────
def bench_verify(n, trials, vk, sigs, sig_bound, sig_bytelen):
    def _deserialize(bs, n):
        BITS = 14; buf = int.from_bytes(bs, 'little'); mask = (1 << BITS) - 1
        return [(buf >> (i*BITS)) & mask for i in range(n)]
    h = _deserialize(vk, n)
    payload = sig_bytelen - HEAD_LEN - SALT_LEN

    VF = ["NTT poly_mul (mul_zq s1·h)", "Norm check", "decode+hash+other"]
    samp = defaultdict(list)

    for sig in sigs:
        salt_v = sig[HEAD_LEN:HEAD_LEN+SALT_LEN]
        enc_s  = sig[HEAD_LEN+SALT_LEN:]

        # decode + hash
        t0 = _rs()
        s1 = decompress(enc_s, payload, n)
        hashed = hash_to_point(msg, salt_v, n)
        t1 = _re(); dh = t1 - t0

        # NTT poly mul: s0 = hashed - s1*h mod q
        t0 = _rs()
        s0 = sub_zq(hashed, mul_zq(s1, h))
        s0 = [(x + (q>>1)) % q - (q>>1) for x in s0]
        t1 = _re(); ntt = t1 - t0

        # norm check
        t0 = _rs()
        ns = sum(x*x for x in s0) + sum(x*x for x in s1)
        _ = ns <= sig_bound
        t1 = _re(); nc = t1 - t0

        samp["NTT poly_mul (mul_zq s1·h)"].append(ntt)
        samp["Norm check"].append(nc)
        samp["decode+hash+other"].append(dh)

    return print_table(f"VERIFY Falcon-{n} ({trials} trials)", VF, samp)

# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n",             type=int, default=512, choices=[512,1024])
    ap.add_argument("--trials",        type=int, default=20)
    ap.add_argument("--keygen-trials", type=int, default=3)
    args = ap.parse_args()
    n, ST, KT = args.n, args.trials, args.keygen_trials

    # Inline param table (avoids Python 3.8 compat issues in falcon.py)
    _params = {
        512:  dict(sigma=165.7366171829776,  sigmin=1.2778336969128337,
                   sig_bound=34034726,  sig_bytelen=666),
        1024: dict(sigma=168.38857144654395, sigmin=1.298280334344292,
                   sig_bound=70265242,  sig_bytelen=1280),
    }
    p           = _params[n]
    sigma       = p["sigma"]
    sigmin      = p["sigmin"]
    sig_bound   = p["sig_bound"]
    sig_bytelen = p["sig_bytelen"]

    print(f"\nFalcon-{n} domain profiler  [falcon.py-master 3]")
    print(f"  Timer  : {RDTSC_SRC}")
    print(f"  CPU    : {CPU_GHZ:.3f} GHz")
    print(f"  keygen : {KT} trials   sign : {ST} valid sigs   verify : {ST} trials")

    # ── Keygen ──
    kg_tot = bench_keygen(n, KT, sigma)

    # ── Build one key for sign/verify ──
    print(f"\nBuilding key for sign/verify …", end="", flush=True)
    from fft import neg
    f, g, F, G = ntru_gen(n)
    B0     = [[g, neg(f)], [G, neg(F)]]
    G0     = gram(B0)
    B0_fft = [[fft(e) for e in row] for row in B0]
    G0_fft = [[fft(e) for e in row] for row in G0]
    T_fft  = ffldl_fft(G0_fft)
    normalize_tree(T_fft, sigma)
    h      = div_zq(g, f)
    def _serialize_h(poly):
        BITS = 14; buf = 0
        for i, c in enumerate(poly): buf ^= c << (i*BITS)
        return buf.to_bytes((len(poly)*BITS+7)>>3, 'little')
    sk = (f, g, F, G, B0_fft, T_fft)
    vk = _serialize_h(h)
    print(" done.")

    msg = b"bench-domains-test-message-xxxx"

    # ── Sign ──
    sg_tot = bench_sign(n, ST, sk, sigma, sigmin, sig_bound, sig_bytelen)

    # ── Collect sigs for verify (inline sign pipeline, no Falcon class) ──
    print(f"\nCollecting {ST} sigs for verify …", end="", flush=True)
    from fft import sub, neg as fneg
    from ffsampling import ffsampling_fft
    int_header = (0x30 + {512:9, 1024:10}[n]).to_bytes(1, "little")
    payload    = sig_bytelen - HEAD_LEN - SALT_LEN
    [[a2,b2],[c2,d2]] = sk[4]
    sigs = []
    for _ in range(ST):
        while True:
            salt2   = urandom(SALT_LEN)
            hashed2 = hash_to_point(msg, salt2, n)
            pf2     = fft(hashed2)
            t0f     = [pf2[i]*d2[i]/q for i in range(n)]
            t1f     = [-pf2[i]*b2[i]/q for i in range(n)]
            z2      = ffsampling_fft([t0f,t1f], sk[5], sigmin, urandom)
            v0f2    = add_fft([z2[0][i]*a2[i] for i in range(n)],
                              [z2[1][i]*c2[i] for i in range(n)])
            v1f2    = add_fft([z2[0][i]*b2[i] for i in range(n)],
                              [z2[1][i]*d2[i] for i in range(n)])
            v02 = [int(round(x)) for x in ifft(v0f2)]
            v12 = [int(round(x)) for x in ifft(v1f2)]
            s02 = sub(hashed2, v02); s12 = fneg(v12)
            ns2 = sum(x*x for x in s02) + sum(x*x for x in s12)
            if ns2 > sig_bound: continue
            enc2 = compress(s12, payload)
            if enc2 is False: continue
            sigs.append(int_header + salt2 + enc2); break
    print(" done.")

    # ── Verify ──
    vf_tot = bench_verify(n, ST, vk, sigs, sig_bound, sig_bytelen)

    # ── Throughput ──
    print(f"\n{'='*72}")
    print(f"  THROUGHPUT SUMMARY  Falcon-{n}  @ {CPU_GHZ:.3f} GHz")
    print(f"{'='*72}")
    print(f"  {'Op':<10}  {'Cycles(med)':>12}  {'Time':>10}  {'ops/sec':>12}")
    print(f"  {'-'*10}  {'-'*12}  {'-'*10}  {'-'*12}")
    for label, cyc in [("keygen", kg_tot), ("sign", sg_tot), ("verify", vf_tot)]:
        print(f"  {label:<10}  {fc(cyc)}  {ft(cyc)}  {CPU_HZ/cyc:>12.1f}")
    print()
