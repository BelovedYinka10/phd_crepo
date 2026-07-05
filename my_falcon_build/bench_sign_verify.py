"""
bench_sign_verify.py — RDTSC phase breakdown for sign and verify
=================================================================
Builds one Falcon-512/1024 key, then profiles every phase of sign()
and verify() with hardware RDTSC (mfence·rdtsc / rdtscp·mfence).
but
Usage:
  python bench_sign_verify.py [--n 512|1024] [--trials 20]
"""

import argparse, ctypes, os, subprocess, sys, tempfile, statistics
from collections import defaultdict

# ── RDTSC ─────────────────────────────────────────────────────────────────────
_RDTSC_C = """\
#include <stdint.h>
uint64_t rs(void){uint64_t a,d;__asm__ __volatile__("mfence\\n\\trdtsc":"=a"(a),"=d"(d));return(d<<32)|a;}
uint64_t re(void){uint64_t a,d;unsigned c;__asm__ __volatile__("rdtscp":"=a"(a),"=d"(d),"=c"(c));__asm__ __volatile__("mfence");return(d<<32)|a;}
"""
def _build():
    td=tempfile.mkdtemp(); src=td+"/r.c"; lib=td+"/r.so"
    open(src,"w").write(_RDTSC_C)
    r=subprocess.run(["cc","-O2","-shared","-fPIC","-o",lib,src],capture_output=True,timeout=10)
    if r.returncode!=0: raise RuntimeError(r.stderr.decode())
    dll=ctypes.CDLL(lib)
    dll.rs.restype=dll.re.restype=ctypes.c_uint64
    return dll.rs, dll.re
try:
    _rs, _re = _build()
    RDTSC_SRC = "hardware RDTSC"
except Exception as e:
    import time
    try: hz=int(subprocess.run(["sysctl","-n","hw.cpufrequency_max"],capture_output=True,text=True,timeout=3).stdout.strip())
    except: hz=2_600_000_000
    ghz=hz/1e9
    print(f"[warn] {e}; fallback perf_counter_ns×{ghz:.2f}")
    def _rs(): return int(time.perf_counter_ns()*ghz)
    _re=_rs; RDTSC_SRC=f"perf_counter_ns×{ghz:.2f}"

sys.path.insert(0, os.path.dirname(__file__))
from params import Q, FALCON_PARAMS, SALT_LEN
from fft import fft, ifft, poly_mul_fft, poly_add_fft
from ffsampling import ffsampling_fft
from sign import hash_to_point, compress, decompress
from poly_arith import poly_mul, poly_center

CPU_HZ=2_600_000_000
try: CPU_HZ=int(subprocess.run(["sysctl","-n","hw.cpufrequency_max"],capture_output=True,text=True,timeout=3).stdout.strip())
except: pass
CPU_GHZ=CPU_HZ/1e9

def fc(c):
    if c>=1e9: return f"{c/1e9:8.3f} Gc"
    if c>=1e6: return f"{c/1e6:8.3f} Mc"
    if c>=1e3: return f"{c/1e3:8.1f} Kc"
    return f"{c:8.0f}  c"
def fu(c):
    us=c/CPU_HZ*1e6
    if us>=1e3: return f"{us/1e3:8.2f} ms"
    return f"{us:8.2f} µs"

def table(title, phases, samples):
    meds={p:statistics.median(samples[p]) for p in phases}
    tot=sum(meds.values())
    W=max(len(p) for p in phases)+2
    print(f"\n  ── {title} ──")
    print(f"  {'Phase':<{W}}  {'Median':>12}  {'Time':>10}  {'Min':>12}  {'Max':>12}  {'%':>6}")
    print(f"  {'-'*W}  {'-'*12}  {'-'*10}  {'-'*12}  {'-'*12}  {'-'*6}")
    for p in phases:
        s=samples[p]; med=meds[p]
        pct=100*med/tot
        print(f"  {p:<{W}}  {fc(med)}  {fu(med)}  {fc(min(s))}  {fc(max(s))}  {pct:5.1f}%")
    print(f"  {'TOTAL':<{W}}  {fc(tot)}  {fu(tot)}")
    worst=max(phases,key=lambda p:meds[p])
    print(f"\n  *** Bottleneck: [{worst}]  {100*meds[worst]/tot:.1f}% ***")
    return tot

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",     type=int,default=512,choices=[512,1024])
    ap.add_argument("--trials",type=int,default=20)
    args=ap.parse_args()
    n=args.n; trials=args.trials
    p=FALCON_PARAMS[n]
    sigmin=p["sigma_min"]; beta_sq=p["beta_sq"]
    payload_len=p["sig_bytelen"]-1-SALT_LEN

    print(f"\nFalcon-{n}  sign/verify phase profiler")
    print(f"  Timer  : {RDTSC_SRC}")
    print(f"  CPU    : {CPU_GHZ:.3f} GHz")
    print(f"  Trials : {trials}")

    print("\nBuilding tree-mode key (ntru_gen running) …", end="", flush=True)
    from keygen import keygen
    sk, pk = keygen(n, mode="tree")
    print(" done.")

    [[a,b],[c_mat,d]]=sk.B_fft
    h=pk.h
    msg=b"rdtsc benchmark message"

    # ── SIGN ──────────────────────────────────────────────────────────────────
    SG_PHASES=["hash_to_point","fft_c","targets (2×poly_mul_fft)","ffsampling",
               "v_fft (4×mul+2×add)","ifft+round (2×)","norm_check","compress"]
    sg=defaultdict(list)

    completed=0; attempts=0
    while completed<trials:
        attempts+=1
        salt=os.urandom(SALT_LEN)

        t0=_rs(); c_poly=hash_to_point(msg,salt,n); t1=_re()
        _h=t1-t0

        t0=_rs(); c_fft=fft([float(x) for x in c_poly]); t1=_re()
        _fc=t1-t0

        t0=_rs()
        tgt0=[x/Q for x in poly_mul_fft(c_fft,d)]
        tgt1=[-x/Q for x in poly_mul_fft(c_fft,b)]
        t1=_re(); _tg=t1-t0

        t0=_rs(); z=ffsampling_fft([tgt0,tgt1],sk.T,sigmin,os.urandom); t1=_re()
        _fs=t1-t0

        t0=_rs()
        v0f=poly_add_fft(poly_mul_fft(z[0],a),poly_mul_fft(z[1],c_mat))
        v1f=poly_add_fft(poly_mul_fft(z[0],b),poly_mul_fft(z[1],d))
        t1=_re(); _vf=t1-t0

        t0=_rs()
        v0=[int(round(x)) for x in ifft(v0f)]
        v1=[int(round(x)) for x in ifft(v1f)]
        t1=_re(); _iv=t1-t0

        s0=[c_poly[i]-v0[i] for i in range(n)]
        s1=[-v1[i] for i in range(n)]

        t0=_rs(); norm_sq=sum(x*x for x in s0)+sum(x*x for x in s1); ok=norm_sq<=beta_sq; t1=_re()
        _nc=t1-t0
        if not ok: continue

        t0=_rs(); enc=compress(s1,payload_len); t1=_re()
        _cp=t1-t0
        if enc is False: continue

        sg["hash_to_point"].append(_h)
        sg["fft_c"].append(_fc)
        sg["targets (2×poly_mul_fft)"].append(_tg)
        sg["ffsampling"].append(_fs)
        sg["v_fft (4×mul+2×add)"].append(_vf)
        sg["ifft+round (2×)"].append(_iv)
        sg["norm_check"].append(_nc)
        sg["compress"].append(_cp)
        completed+=1

    if attempts>trials:
        print(f"\n  (rejection rate {(attempts-trials)/attempts*100:.1f}%  — {attempts} attempts for {trials} valid sigs)")

    sg_total=table(f"SIGN  Falcon-{n}  ({trials} trials)", SG_PHASES, sg)

    # ── generate sigs for verify ──────────────────────────────────────────────
    from sign import sign as fsign
    print(f"\nGenerating {trials} sigs for verify …", end="", flush=True)
    sigs=[fsign(sk,msg) for _ in range(trials)]
    print(" done.")

    # ── VERIFY ────────────────────────────────────────────────────────────────
    VF_PHASES=["decompress","hash_to_point","poly_mul_ntt (s1·h)","center+norm_check"]
    vf=defaultdict(list)

    for sig in sigs:
        salt_v=sig[1:1+SALT_LEN]; enc_s1=sig[1+SALT_LEN:]

        t0=_rs(); s1=decompress(enc_s1,payload_len,n); t1=_re()
        vf["decompress"].append(t1-t0)

        t0=_rs(); c=hash_to_point(msg,salt_v,n); t1=_re()
        vf["hash_to_point"].append(t1-t0)

        t0=_rs()
        s1m=[x%Q for x in s1]; s1h=poly_mul(s1m,h)
        s0=[(c[i]-s1h[i])%Q for i in range(n)]
        t1=_re(); vf["poly_mul_ntt (s1·h)"].append(t1-t0)

        t0=_rs()
        s0c=poly_center(s0)
        ns=sum(x*x for x in s0c)+sum(x*x for x in s1); _=ns<=beta_sq
        t1=_re(); vf["center+norm_check"].append(t1-t0)

    vf_total=table(f"VERIFY  Falcon-{n}  ({trials} trials)", VF_PHASES, vf)

    # ── Summary ────────────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"  THROUGHPUT SUMMARY  Falcon-{n}  @ {CPU_GHZ:.2f} GHz")
    print(f"{'='*60}")
    print(f"  {'Op':<10}  {'Cycles (med)':>14}  {'Time':>10}  {'ops/s':>10}")
    print(f"  {'-'*10}  {'-'*14}  {'-'*10}  {'-'*10}")
    for label,cyc in [("sign",sg_total),("verify",vf_total)]:
        print(f"  {label:<10}  {fc(cyc)}  {fu(cyc)}  {CPU_HZ/cyc:>8.1f}/s")

if __name__=="__main__":
    main()
