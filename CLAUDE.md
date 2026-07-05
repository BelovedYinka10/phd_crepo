# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

PhD research on the **Falcon** post-quantum lattice-based digital signature scheme (NIST PQC standard). Ring: R = Z[x]/(x^n+1), q = 12289. Security levels: Falcon-512 (Cat 1), Falcon-1024 (Cat 5).

Subdirectory overview:

| Directory | Language | Purpose |
|---|---|---|
| `falcon-round3/` | C (C99) | Official NIST Round 3 submission (Thomas Pornin / NCC Group) |
| `falcon.py-master 3/` | Python | Educational reference implementation (Thomas Prest / PQShield) |
| `my_falcon_build/` | Python + Verilog/SV | Personal from-scratch implementation (active development) |
| `cshake-core/` | Verilog | cSHAKE/SHAKE-256/128 hardware core (Yale/Jungk) — reference RTL |
| `FalconSign/` | SystemVerilog | High-speed Falcon hardware: SamplerZ + Arithmetic (ZCU104) |
| `XKCP/` | C / Python | Official eXtended Keccak Code Package — authoritative reference |

## Build & Test Commands

### C implementation (`falcon-round3/Extra/c/`)
```bash
cd falcon-round3/Extra/c/
make                # builds test_falcon and speed binaries
./test_falcon       # self-tests
./speed             # benchmarks (Falcon-256, 512, 1024)
```

NIST reference/optimized builds (for KAT generation):
```bash
cd falcon-round3/Reference_Implementation/falcon512/falcon512int/
make                # → build/kat512int
```

### Python reference (`falcon.py-master 3/`)
```bash
pip install pycryptodome beartype numpy   # required deps
cd "falcon.py-master 3/"
make test           # runs python3 test.py (n=64..1024)
```

### cSHAKE core (`cshake-core/`)
Configurable PARALLEL_SLICES (1,2,4,8,16,32) — set in `verilog/keccak_pkg.v`.
```bash
# Step 1 — generate test vectors
cd cshake-core/testvectors_generator
make NUMTESTS=10 HDL="verilog" cshake256   # produces verilog/testvectors.v

# Step 2 — simulate (Vivado XSim)
cd cshake-core/Vivado
make simulate      # run and check pass/fail
make waveform      # open waveform for debugging

# Full regression (all SHAKE/cSHAKE variants × all PARALLEL_SLICES)
cd cshake-core/regression_tests
./run_tests.sh
```

### FalconSign (`FalconSign/`)
SamplerZ + Arithmetic modules for Falcon hardware signing on Zynq UltraScale+ (ZCU104, Vivado 2022.2).
```bash
# Floating-point unit tests
cd FalconSign/test
make               # generates fp_vectors.txt then runs fp_tb.sv via Verilator

# SamplerZ simulation
cd FalconSign/sampler_sim
make               # runs samplerz_tb.sv + statistical test

# Vivado synthesis / simulation
cd FalconSign/vivado
make               # create_project.tcl + run_sim.tcl
```

### XKCP (`XKCP/`)
Authoritative Keccak reference. The most useful files for cross-checking RTL:
- `lib/low/KeccakP-1600/ref-64bits/KeccakP-1600-reference.c` — round-by-round reference
- `Standalone/CompactFIPS202/C/Keccak-readable-and-compact.c` — compact readable FIPS 202
- `Standalone/CompactFIPS202/Python/CompactFIPS202.py` — Python reference
```bash
cd XKCP
make           # builds all targets per Makefile.build
```

### Personal build (`my_falcon_build/`)
No external dependencies (stdlib only). No Makefile for Python.
```bash
cd my_falcon_build/
python test_ntt.py          # NTT correctness + cross-backend consistency
python test_poly_arith.py   # polynomial arithmetic
python test_samplerz.py     # discrete Gaussian sampler (statistical tests)
python test_fft.py          # FFT round-trip, mul, split/merge, adj, LDL (11 tests)
python test_ntrugen.py      # NTRU keygen: karatsuba, field_norm, ntru_solve, ntru_gen_512 (10 tests)
python test_ffsampling.py   # FFT sampling: gram, ldl, ffldl, normalize, ffnp, ffsampling, full pipeline (7 tests)
python test_keygen.py       # keygen: poly_div_ntt, key structure, h*f≡g, NTRU eq, LDL tree, B_fft, serialize, dynamic/expand (9 tests)
python test_sign.py         # sign/verify: hash_to_point, compress/decompress, sign, verify, tamper rejection (7 tests; n=1024 test skipped by default)
python ntt.py               # standalone NTT test + benchmark
python bench_karatsuba_vs_toom3_vs_mutualized.py --trials 10  # multiplication backend benchmark
python bench_halftree.py    # LDL half-tree optimization benchmark (storage + timing + symmetry check)
python bench_rdtsc.py [--n 512|1024] [--trials 20] [--keygen-trials 5]  # median cycles+µs+% for every phase
python bench_sign_verify.py [--n 512|1024] [--trials 20]                # RDTSC sign/verify breakdown
python bench_ffsampling_internal.py [--n 512|1024] [--trials 30]        # per-op drill-down inside ffsampling_fft
python bench_python_domains.py [--n 512|1024] [--trials 20] [--keygen-trials 5]  # domain table matching C bench_domains
```

### RTL (`my_falcon_build/rtl/`)
Modules: `base_sampler.v` (combinational BaseSampler), `shake256.sv` + `keccak_f1600.sv` (SHAKE256 XOF), `hash_to_point.sv` (SHAKE256 → Z_q^n coefficient stream).

```bash
cd my_falcon_build/rtl/
make test          # generate vectors.txt then simulate with Verilator
make vectors       # (re)generate vectors.txt only  (NVEC=2000 default)
make sim           # build + run testbench against existing vectors.txt
make clean

# Synthesis (batch Vivado):
cd vivado/
# Edit PART in run_vivado.tcl to your target device, then:
vivado -mode batch -source run_vivado.tcl
# Reports: fmax.rpt (Fmax/WNS), util.rpt (LUT/FF), timing.rpt
```

## Architecture (`my_falcon_build/`)

Building bottom-up. **Keygen + sign/verify are complete. Next: optimisation, side-channel hardening, or hardware integration.**

> `my_falcon_build/CLAUDE.md` contains a more detailed per-module reference.

### Dependency order (bottom → top)

```
params.py
  ├─ ntt.py          (uses Q, psi tables from params)
  │    └─ poly_arith.py    (wraps NTT)
  ├─ samplerz.py     (uses RCDT table from params)
  └─ fft.py          (uses GM_TAB twiddle table)
       └─ ntrugen.py       (uses fft for poly ops)
            └─ ffsampling.py    (uses fft split/merge/ldl)
                 └─ keygen.py       (uses ntrugen + ffsampling + ntt)
                      └─ sign.py        (uses keygen + ffsampling + fft + poly_arith)
```

### Key data-flow

```
keygen:  ntru_gen(n) → (f,g,F,G)
         B = [[g,-f],[G,-F]]
         gram(B) → G_mat  →  ffldl_fft(G_fft) → normalize_tree(T, sigma)
         h = poly_div_ntt(g, f)   # public key

sign:    salt ← os.urandom(40)
         c    = hash_to_point(salt||msg, n)  # SHAKE256 rejection → Z_q^n
         t0   = c·(−F)/q,  t1 = c·f/q      # via B_fft in FFT domain
         z    = ffsampling_fft([t0,t1], T, sigmin, rng)
         v    = z·B → ifft → round          # nearest lattice point
         (s0,s1) = (c−v0, −v1); check ||s||²≤β²; compress(s1) → sig

verify:  c    = hash_to_point(salt||msg, n)
         s0   = c − s1·h mod q (centered)
         accept iff ||s0||²+||s1||² ≤ β²
```

### LDL tree structure

`ffldl_fft` returns a nested list used by `ffsampling_fft`:
- Internal node: `[L10, T_left, T_right]` (len 3)
- Bottom internal (leaf): `[L10_2, [sigma0, 0], [sigma1, 0]]` (len 3, T[1] has len 2)
- After `normalize_tree`: leaf values are `sigma / sqrt(D)`, not raw D values

`keygen(n, mode)` has two modes:
- `"tree"` — precomputes `B_fft` and `T` at keygen; fast signing
- `"dynamic"` — stores only `(f,g,F,G)`; rebuilds tree per signature via `expand_secret_key`

### Half-tree optimization (`ffldl_fft_half`, `ffsampling_fft_half`)

Implemented in `ffsampling.py` based on the 2022 Springer paper. The algebraic identity `D[0][0] · D[1][1] = q²` (from the NTRU equation `fG - gF = q`) forces a **bit-complement symmetry** across the tree:

> For every node in T0 at path P, the node in T1 at the bit-complement path P̄ satisfies `L10(T1_P̄) = -L10(T0_P)`.

This means T1's L10 values are entirely derivable from T0 by negation + child-swap. `ffldl_fft_half` builds only T0; `ffsampling_fft_half(t, L10_root, T0, ...)` reconstructs T1 on-the-fly via `_negate_l10_tree`. Empirically verified for n=512: L10 reconstruction error < 1e-12. **Note:** sigma (leaf) values in T1 are NOT derivable from T0's sigmas — they differ by ~30% and must be stored or recomputed separately. Storage reduction for L10 values alone is ~40%.

### RTL (`my_falcon_build/rtl/`)

| Module | Type | Description |
|---|---|---|
| `base_sampler.v` | Verilog (combinational) | RCDT inversion: 72-bit `u` → 5-bit `z0`; 18 parallel comparisons + popcount; bit-for-bit identical to `samplerz.py:basesampler()` |
| `keccak_f1600.sv` | SystemVerilog | Keccak-f[1600] permutation (rounds loop) |
| `shake256.sv` | SystemVerilog | SHAKE256 XOF: absorb byte stream → squeeze byte stream |
| `hash_to_point.sv` | SystemVerilog | Wraps `shake256`; rejection-samples 16-bit words to emit n coefficients in Z_q |

`gen_vectors.py` produces `vectors.txt` for the BaseSampler testbench. `bench_signing_share.py` measures BaseSampler's fraction of total ffsampling time.

The `vivado/` subdirectory contains the synthesis kit: registered wrapper (`base_sampler_top.v`), clock constraint (`base_sampler.xdc`), and batch TCL flow (`run_vivado.tcl`). Reported Fmax covers only the 18-comparator datapath; `hash_to_point.sv` and SHAKE256 are not yet included in synthesis.

## Critical Mathematical Constants

- **NTT primitive root: g = 11** (NOT 7). ord(11) = 12288 = q−1. ord(7) = 2048 ≠ q−1.
- ψ = 11^((q−1)/(2n)) mod q: for n=512 → ψ=10302; for n=1024 → ψ=1945
- Montgomery: R=2^16, Q0I=12287, R2=10952
- FFT uses split-halves layout: `f[0..n/2-1]` = real parts, `f[n/2..n-1]` = imaginary parts. The first butterfly level (twiddle = i) is implicit; iFFT scales by 2/n.

## Signature formats (C API)
Three formats: `FALCON_SIG_COMPRESSED` (variable length), `FALCON_SIG_PADDED` (fixed length), `FALCON_SIG_CT` (constant-time, fixed length).
