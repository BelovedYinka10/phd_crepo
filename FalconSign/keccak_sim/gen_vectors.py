#!/usr/bin/env python3
"""
gen_vectors.py — Generate SHAKE256 KAT vectors using XKCP Python reference.
Writes shake256_vectors.hex (one byte per line, 2 hex chars).

Layout:
  [4 bytes] num_tests  (big-endian)
  per test:
    [4 bytes] msg_len
    [4 bytes] out_len
    [msg_len bytes] message
    [out_len bytes] SHAKE256(message)[0:out_len]
"""
import sys, os

XKCP = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     '../../XKCP/Standalone/CompactFIPS202/Python')
sys.path.insert(0, XKCP)
from CompactFIPS202 import SHAKE256

tests = [
    # (message,            out_bytes)   what it exercises
    (b'\x61',             32),          # 1 byte  "a"
    (b'\x61\x62\x63',     32),          # 3 bytes "abc"
    (bytes(range(17)),    32),          # 17 bytes (< 1 rate block)
    (bytes(range(136)),   32),          # exactly 1 rate block (136 B)
    (bytes(range(200)),   64),          # multi-block absorb + multi-squeeze
]

data = bytearray()

def w32(n):
    data.extend(n.to_bytes(4, 'big'))

w32(len(tests))
for msg, out_len in tests:
    expected = SHAKE256(msg, out_len)
    w32(len(msg))
    w32(out_len)
    data.extend(msg)
    data.extend(expected)

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'shake256_vectors.hex')
with open(out, 'w') as f:
    for b in data:
        f.write(f'{b:02x}\n')

print(f'Generated {out}  ({len(data)} bytes, {len(tests)} tests)')
for i, (msg, out_len) in enumerate(tests):
    h = SHAKE256(msg, out_len)
    print(f'  test {i}: msg={msg[:4].hex()}... ({len(msg)}B) → {h[:8].hex()}...')
