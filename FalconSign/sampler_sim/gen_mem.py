"""
gen_mem.py — Generate memory.hex for the samplerz_tb testbench.

Memory layout (256-bit words, address = word index):
  Address 0  : mu values  — {mu_r[63:0], mu_l[63:0], 128'b0}
               mu = 0.0 → both doubles are 0x0000000000000000
  Address 2  : isigma values — {isigma1024, isigma512, isigma1024, isigma512}
               pre_samp reads r_data[64*isigma_index+:64] for index 0..3
  Address 130: ChaCha20 PRNG seed (256-bit key) — any non-zero value works
  Address 131: ChaCha20 nonce / counter

All other addresses: 0 (safe default for double = 0.0).

Task wiring (from samplerz.sv):
  src0_addr → mu   at word 0  (address 0)
  src1_addr → isigma at word 2 (address 2)
  dst_addr  → sample output at word 64 (address 64)
"""

import struct, sys

DEPTH      = 4 * 2048   # BANK_NUM * BANK_DEPTH
WORD_BITS  = 256
WORD_BYTES = WORD_BITS // 8

mem = [0] * DEPTH   # 256-bit words, stored as Python ints

def pack_doubles(*vals):
    """Pack 1-4 IEEE 754 doubles into a 256-bit word (little-endian, d[0] at LSB)."""
    result = 0
    for i, v in enumerate(vals):
        bits = struct.unpack('Q', struct.pack('d', v))[0]
        result |= bits << (64 * i)
    return result

# ── Address 0: mu values ──────────────────────────────────────────────────────
# pre_samp reads r_data[127:0] = {fpr_mu_r[63:0], fpr_mu_l[63:0]}
# so double at bits [63:0]   = mu_l
#    double at bits [127:64] = mu_r
MU = 0.0
mem[0] = pack_doubles(MU, MU, MU, MU)   # all four positions = 0.0

# ── Address 2: isigma values ──────────────────────────────────────────────────
# pre_samp reads r_data[64*isigma_index+:64] for isigma_index in {0,1}
# index 0 → sigma for the current sample
# index 1 → sigma for the second sample
# Use the same isigma for both (Falcon-512)
# Use sigma=1.5 — a representative leaf-level sigma (acceptance rate ~80%)
# This lets us clearly distinguish: correct dist (std≈1.5) vs always-accept (std≈1.82) vs always-reject (std≈0.5)
SIGMA_TEST  = 1.5
ISIGMA_TEST = 1.0 / SIGMA_TEST   # 0x3fe5555555555555

mem[2] = pack_doubles(ISIGMA_TEST, ISIGMA_TEST, ISIGMA_TEST, ISIGMA_TEST)

# ── Address 130: ChaCha20 key (init_state[255:0]) ────────────────────────────
# chacha20.sv: if(cnt == 1+SAMPLERZ_READ_DELAY) init_state[255:0] <= mem_rd_chacha20_data
# 256 bits = 8 × 32-bit key words.  Any non-zero key produces good randomness.
# Use a well-known test vector key to get reproducible, verifiable output.
KEY = 0x03020100_07060504_0B0A0908_0F0E0D0C_13121110_17161514_1B1A1918_1F1E1D1C
mem[130] = KEY & ((1 << 256) - 1)

# ── Address 131: nonce + counter (init_state[383:256] and cc) ────────────────
# chacha20.sv: if(cnt == 2+SAMPLERZ_READ_DELAY) init_state[383:256] <= mem_rd_chacha20_data[127:0]
#              cc[i] <= mem_rd_chacha20_data[191:128] + i
# Nonce at [127:0]: 4 × 32-bit nonce words
# Counter at [191:128]: 64-bit initial counter value
NONCE   = 0x09000000_4A000000_00000000_00000000  # RFC 8439 test nonce
COUNTER = 1  # initial counter = 1
mem[131] = (COUNTER << 128) | NONCE

# ── Write memory.hex (one 256-bit word per line, 64 hex digits, big-endian) ──
out = sys.argv[1] if len(sys.argv) > 1 else "memory.hex"
with open(out, "w") as f:
    for word in mem:
        f.write(f"{word:064X}\n")

print(f"Written {DEPTH} words ({DEPTH * WORD_BYTES // 1024} KB) to {out}")
print(f"  addr  0: mu    = {MU}")
print(f"  addr  2: isigma_test = {ISIGMA_TEST:.10f}  (sigma={SIGMA_TEST})")
print(f"  addr 64: dst (output samples will be written here)")
print(f"  addr130: ChaCha20 key  = {hex(KEY)[:18]}...")
print(f"  addr131: nonce+counter = {hex(NONCE)[:18]}...")
