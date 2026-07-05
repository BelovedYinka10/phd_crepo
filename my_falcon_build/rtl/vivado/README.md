# BaseSampler — Vivado synthesis kit

Get real Fmax + area for the Falcon BaseSampler datapath on your target FPGA.

## Files to bring to Vivado
- `../base_sampler.v`   — the combinational core (RTL under test)
- `base_sampler_top.v`  — registered wrapper (synth/impl THIS as top)
- `base_sampler.xdc`    — clock constraint for Fmax
- `run_vivado.tcl`      — batch flow: synth → impl → reports

The testbench `../base_sampler_tb.v` is **simulation only** (file I/O) — not
needed for synthesis. You can run it in Vivado XSim if you also copy the
`vectors.txt` produced by `../gen_vectors.py`.

## Quick run (batch)
1. Edit `PART` in `run_vivado.tcl` to your device.
2. From this directory:
   ```
   vivado -mode batch -source run_vivado.tcl
   ```
3. Read results:
   - `fmax.rpt`   — Fmax (MHz), WNS
   - `util.rpt`   — LUT / FF utilization
   - `timing.rpt` — full timing summary

## GUI alternative
Create project → add `base_sampler.v` + `base_sampler_top.v` (top) →
add `base_sampler.xdc` → set part → Run Synthesis → Run Implementation →
Reports: Timing Summary (WNS) + Utilization. Fmax = 1000 / (period_ns − WNS_ns).

## To find the TRUE ceiling
Lower `PERIOD` in `base_sampler.xdc` (and the matching value in
`run_vivado.tcl`) until WNS goes slightly negative — that period is the
critical-path limit.

## IMPORTANT honesty notes
- This characterises **only** the 18-comparator + popcount datapath. The real
  BaseSampler needs a PRNG (SHAKE256 / ChaCha20) to generate `u`; that block
  dominates area and sets the true system Fmax. Quote this number as an
  **upper bound on the datapath**, not the whole sampler.
- Throughput is **1 sample/cycle** by construction, so:
  `samples/sec = Fmax`. Compare against the measured scalar-C software rate
  (~48 M gaussian0/sec, i.e. 20.8 ns/call on Apple M-series) for the speedup —
  but use a C number from comparable hardware for a fair claim.
