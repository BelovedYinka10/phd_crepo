# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Personal from-scratch Python + Verilog implementation of **Falcon**, a NIST-standardized post-quantum lattice-based digital signature scheme. Ring: R = Z[x]/(x^n+1), q = 12289. Security levels: Falcon-512 (Cat 1), Falcon-1024 (Cat 5).

**Status: keygen + sign/verify are complete. Next stage: optimisation, side-channel hardening, or hardware integration.**

No external dependencies (stdlib only), except `gen_keygen_arch.py` which requires `python-docx` to regenerate `Falcon_Keygen_Architecture.docx`. No Makefile for Python.

## Test Commands

```bash
python test_ntt.py           # NTT correctness + cross-backend consistency
python test_poly_arith.py    # polynomial arithmetic
python test_samplerz.py      # discrete Gaussian sampler (statistical tests)
python test_fft.py           # FFT round-trip, mul, split/merge, adj, LDL (11 tests)
python test_ntrugen.py       # NTRU keygen: karatsuba, field_norm, ntru_solve, ntru_gen_512 (10 tests)
python test_ffsampling.py    # FFT sampling: gram, ldl, ffldl, normalize, ffnp, ffsampling, full pipeline (7 tests)
python test_keygen.py        # keygen: poly_div_ntt, key structure, h*f≡g, NTRU eq, LDL tree, B_fft, serialize, dynamic/expand (9 tests)
python test_sign.py          # sign/verify: hash_to_point, compress/decompress, sign, verify, tamper rejection (7 tests; test 8 skipped by default — n=1024, slow)

python check_keyes.py        # quick keygen smoke-test: prints sk.f/sk.g/pk.h (not in test suite)
```

## Benchmarks

```bash
python ntt.py                                                          # NTT + Montgomery backend timing
python bench_karatsuba_vs_toom3_vs_mutualized.py --trials 10          # poly mul backends
python bench_halftree.py                                               # half-tree storage + timing + symmetry check

# Cycle-level profiling (hardware RDTSC; falls back to perf_counter_ns on non-x86)
python bench_rdtsc.py [--n 512|1024] [--trials 20] [--keygen-trials 5]
# → median cycles+µs+% for every Falcon phase: keygen, sign, verify

python bench_sign_verify.py [--n 512|1024] [--trials 20]
# → RDTSC phase table for sign and verify only (finer breakdown + bottleneck callout)

python bench_ffsampling_internal.py [--n 512|1024] [--trials 30]
# → per-operation drill-down inside ffsampling_fft: samplerz, l10_mulsub, poly_split/merge
#   plus per-depth wall-time table (depth 0 = root, log2(n)−1 = leaf)

python bench_python_domains.py [--n 512|1024] [--trials 20] [--keygen-trials 5]
# → domain table matching C bench_domains: BaseSampler, berexp, FFT ops, hash/compress, NTT

cd rtl/
python bench_basesampler.py                                            # standalone basesampler throughput
python bench_signing_share.py [--passes 200] [--n 512|1024]           # BaseSampler % of ffsampling time
# (both import from parent dir via sys.path.insert; run from rtl/)
```

## RTL (Verilator simulation + Vivado synthesis)

```bash
cd rtl/
make test          # generate vectors.txt then simulate with Verilator
make vectors       # (re)generate vectors.txt only (NVEC=2000 default)
make sim           # build + run testbench against existing vectors.txt
make clean

# Vivado synthesis (batch):
cd rtl/vivado/
# Edit PART in run_vivado.tcl to your target device, then:
vivado -mode batch -source run_vivado.tcl
# Reports: fmax.rpt, util.rpt, timing.rpt
```

## Module Dependency Order

```
params.py                          # Q, sigma tables, header bytes, encoding constants
  ├─ ntt.py                        # negacyclic NTT/iNTT over Z_q; two backends
  │    └─ poly_arith.py            # polynomial arithmetic wrapping NTT
  ├─ samplerz.py                   # discrete Gaussian sampler (RCDT + berexp)
  └─ fft.py                        # complex FFT/iFFT in split-halves layout
       └─ ntrugen.py               # NTRU keygen: karatsuba, field_norm, ntru_solve, ntru_gen
            └─ ffsampling.py       # LDL tree construction + Gaussian sampling
                 └─ keygen.py      # ties everything together; outputs SecretKey + PublicKey
```

## Architecture

### Two NTT backends (`ntt.py`)

- **Plain `%`**: `ntt()` / `intt()` — readable, uses Python `%` for reduction
- **Montgomery**: `ntt_montgomery()` / `intt_montgomery()` — shift-based, no division, FPGA-ready

Both produce identical outputs. Montgomery uses R=2^16, Q0I=12287, R2=10952. Primitive root is **g=11** (not 7 — `ord(7)=2048 ≠ q−1`; `ord(11)=12288 = q−1`).

Twiddle tables `GMB`/`IGMB` are precomputed at import in bit-reversed order to match the iterative CT/GS loops. Both n=512 and n=1024 tables are built at import time.

### FFT layout (`fft.py`)

All signing-side FFT arrays use **split-halves layout**:
- `f[0..n/2-1]` = real parts, `f[n/2..n-1]` = imaginary parts
- One flat list of n floats stores n/2 complex values
- Twiddle table `GM_TAB` is bit-reversed; the first butterfly level (twiddle = i) is implicit
- iFFT scales by 2/n

### Keygen pipeline (`keygen.py`)

```
ntru_gen(n) → (f, g, F, G)
B = [[g, -f], [G, -F]]
gram(B) → G_mat → fft() → ffldl_fft() → normalize_tree(T, sigma)
h = poly_div_ntt(g, f)    # public key: g·f⁻¹ mod q
```

`keygen(n, mode)` has two modes:
- `"tree"` — precomputes `B_fft` and `T` at keygen time; `sk.B_fft` and `sk.T` are populated; fast signing
- `"dynamic"` — stores only `(f,g,F,G)`; `sk.B_fft=None`, `sk.T=None`; call `expand_secret_key(sk)` before signing

### LDL tree structure (`ffsampling.py`)

`ffldl_fft()` returns nested lists; shape matters for traversal:
- Internal node: `[L10, T_left, T_right]` — `len == 3`
- Leaf (before normalize): `[L10, D00, D11]` — `len == 3` but `len(D00) == 2`
- After `normalize_tree(T, sigma)`: leaves become `[sigma/sqrt(D), 0]` — `len == 2`

The guard `len(T[1]) == 2` distinguishes the bottom internal node (leaf level) from higher internal nodes in `ffsampling_fft` and `ffnp_fft`.

### Half-tree optimization (`ffldl_fft_half` / `ffsampling_fft_half`)

The NTRU equation `fG − gF = q` implies `D[0][0] · D[1][1] = q²` pointwise in FFT domain. This forces a bit-complement symmetry: `L10(T1) = −L10(T0)` at every tree node.

- `ffldl_fft_half(G)` → returns `(L10_root, T0, T1)` — builds both subtrees; caller normalizes both, then calls `_extract_sigma_tree(T1)` to keep only T1's sigmas and discard its L10s
- `ffsampling_fft_half(t, L10_root, T0, T1_sigma_tree, ...)` → derives T1 on-the-fly via `_negate_l10_tree(T0, T1_sigma_tree)`
- Storage reduction for L10 values: ~40% for n=512
- **Sigma (leaf) values in T1 are NOT derivable from T0** — they differ by ~25%; must be stored via `_extract_sigma_tree(T1)` after normalization

### Signing pipeline (`sign.py`)

```
salt  ← os.urandom(40)
c     = hash_to_point(salt||msg, n)     # SHAKE256 rejection-sampling → Z_q^n
t0    = c·(−F)/q,  t1 = c·f/q          # target in FFT domain via B_fft
z     = ffsampling_fft([t0,t1], T, sigmin, rng)
v     = z·B_fft → ifft → round         # nearest lattice point
(s0, s1) = (c−v0, −v1)                 # short signature vector
check ||s||² ≤ β²; compress(s1) → bytes
```

`verify()` recovers `s0 = c − s1·h mod q` (via NTT poly_mul), checks `||s0||² + ||s1||² ≤ β²`.
Signature format: `header(1) + salt(40) + compressed_s1(payload_len bytes)`.

### Sampler pipeline (`samplerz.py`)

```
basesampler()       # RCDT inversion: 72-bit u → z0 ∈ {0..18} from half-Gaussian(σ_max=1.8205)
approxexp(x,ccs)    # integer approx of 2^63·ccs·exp(-x); degree-12 polynomial via Horner
berexp(x,ccs)       # Bernoulli trial with prob ≈ ccs·exp(-x); bit-by-bit comparison
samplerz(μ,σ,σ_min) # rejection sampling combining the above
```

### RTL (`rtl/`)

`base_sampler.v` implements `basesampler()` as a **combinational** module:
- Input: 72-bit `u`; Output: 5-bit `z0`
- 18 parallel comparisons against RCDT entries → popcount tree
- Bit-for-bit identical to `samplerz.py:basesampler()`

`rtl/vivado/` contains the synthesis kit: `base_sampler_top.v` (registered wrapper), `base_sampler.xdc` (clock constraint), `run_vivado.tcl` (batch flow). Reported Fmax covers only the 18-comparator datapath — the PRNG (SHAKE256/ChaCha20) is not yet implemented in RTL.

Fmax from reports: `Fmax (MHz) = 1000 / (period_ns − WNS_ns)`. To find the true datapath ceiling, lower `PERIOD` in `base_sampler.xdc` (and the matching value in `run_vivado.tcl`) until WNS goes slightly negative. Throughput is 1 sample/cycle by construction, so samples/sec = Fmax. Software baseline for comparison: ~48 M basesampler calls/sec (≈20.8 ns/call) on Apple M-series; use a number from comparable hardware for a fair claim.

## Critical Constants

| Constant | Value | Note |
|---|---|---|
| Q (modulus) | 12289 | prime = 1 + 12·1024 |
| NTT root | g = **11** | ord(11)=12288=q−1; ord(7)=2048≠q−1 |
| ψ for n=512 | 10302 | 11^((q−1)/1024) mod q |
| ψ for n=1024 | 1945 | 11^((q−1)/2048) mod q |
| Montgomery R | 2^16 | Q0I=12287, R2=10952 |
| RCDT precision | 72 bits | 9 bytes/sample; 18 table entries |
| sigma_max | 1.8205 | RCDT table built for this σ |

## Key encoding

Public key: 1 header byte (`0x00 | log2n`) + 14 bits per coefficient (n coefficients, packed).
Private key: header `0x50 | log2n`; stores `f`, `g`, `F` (bit-width from `FG_BITWIDTH[n]`); G is recomputed.
Signature: header `0x30 | log2n`; compressed encoding with 7 low bits binary + unary high bits.
