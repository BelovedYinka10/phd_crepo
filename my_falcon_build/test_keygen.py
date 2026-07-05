"""
test_keygen.py — Tests for Falcon key generation module
========================================================
"""

import math
import time

from keygen import (keygen, SecretKey, PublicKey, poly_div_ntt,
                    expand_secret_key,
                    serialize_public_key, deserialize_public_key)
from ntt import ntt, intt, ntt_mul
from ntrugen import karamul
from fft import fft, ifft
from params import Q, FALCON_PARAMS


# ── Helpers ──────────────────────────────────────────────────────────────────

def _tree_depth(T):
    """Count the depth of an LDL tree (leaves are len-2 arrays)."""
    if len(T) == 2:
        return 0
    return 1 + max(_tree_depth(T[1]), _tree_depth(T[2]))


def _tree_leaf_count(T):
    """Count leaves in an LDL tree."""
    if len(T) == 2:
        return 1
    return _tree_leaf_count(T[1]) + _tree_leaf_count(T[2])


# ── Tests ────────────────────────────────────────────────────────────────────

def test_poly_div_ntt():
    """h = g/f mod q  →  h*f ≡ g mod q."""
    from ntrugen import gen_poly
    n = 512
    f = gen_poly(n)
    g = gen_poly(n)

    # Ensure f is invertible (gen_poly should produce invertible polys)
    f_ntt = ntt([c % Q for c in f])
    if any(v == 0 for v in f_ntt):
        print("[SKIP] poly_div_ntt (f not invertible)")
        return

    h = poly_div_ntt(g, f)

    # Verify: h * f ≡ g (mod q, x^n+1)
    h_ntt = ntt(h)
    f_ntt = ntt([c % Q for c in f])
    product_ntt = ntt_mul(h_ntt, f_ntt)
    product = intt(product_ntt)

    g_mod = [c % Q for c in g]
    for i in range(n):
        assert product[i] == g_mod[i], \
            f"h*f[{i}] = {product[i]} != g[{i}] = {g_mod[i]}"

    print("[OK]  poly_div_ntt")


def test_keygen_structure():
    """keygen returns correct SecretKey and PublicKey structures."""
    print("       keygen structure running...", end="", flush=True)
    t0 = time.time()

    n = 512
    sk, pk = keygen(n)

    # SecretKey fields
    assert isinstance(sk, SecretKey)
    assert sk.n == n
    assert sk.mode == "tree"
    assert len(sk.f) == n
    assert len(sk.g) == n
    assert len(sk.F) == n
    assert len(sk.G) == n
    assert len(sk.B_fft) == 2
    assert len(sk.B_fft[0]) == 2
    assert len(sk.B_fft[0][0]) == n
    assert sk.T is not None

    # PublicKey fields
    assert isinstance(pk, PublicKey)
    assert pk.n == n
    assert len(pk.h) == n
    # h coefficients should be in [0, q)
    assert all(0 <= c < Q for c in pk.h), "h coefficients out of range"

    dt = time.time() - t0
    print(f"\r[OK]  keygen structure (n=512) — {dt:.1f}s")


def test_public_key_h_f_eq_g():
    """Verify h*f ≡ g (mod q, x^n+1)."""
    print("       h*f ≡ g running...", end="", flush=True)
    t0 = time.time()

    n = 512
    sk, pk = keygen(n)

    h_ntt = ntt(pk.h)
    f_ntt = ntt([c % Q for c in sk.f])
    product_ntt = ntt_mul(h_ntt, f_ntt)
    product = intt(product_ntt)

    g_mod = [c % Q for c in sk.g]
    for i in range(n):
        assert product[i] == g_mod[i], \
            f"h*f[{i}] = {product[i]} != g[{i}] = {g_mod[i]}"

    dt = time.time() - t0
    print(f"\r[OK]  h*f ≡ g (mod q) (n=512) — {dt:.1f}s")


def test_ntru_equation():
    """Verify f*G - g*F = q (mod x^n+1)."""
    print("       NTRU equation running...", end="", flush=True)
    t0 = time.time()

    n = 512
    sk, _ = keygen(n)

    fG = karamul(sk.f, sk.G)
    gF = karamul(sk.g, sk.F)

    # f*G - g*F should equal [q, 0, 0, ..., 0]
    diff = [fG[i] - gF[i] for i in range(n)]
    assert diff[0] == Q, f"Constant term = {diff[0]}, expected {Q}"
    for i in range(1, n):
        assert diff[i] == 0, f"diff[{i}] = {diff[i]}, expected 0"

    dt = time.time() - t0
    print(f"\r[OK]  NTRU equation (n=512) — {dt:.1f}s")


def test_tree_structure():
    """LDL tree has correct depth and leaf count."""
    print("       tree structure running...", end="", flush=True)
    t0 = time.time()

    n = 512
    sk, _ = keygen(n)
    T = sk.T

    expected_depth = int(math.log2(n))
    actual_depth = _tree_depth(T)
    assert actual_depth == expected_depth, \
        f"Depth {actual_depth} != expected {expected_depth}"

    expected_leaves = n
    actual_leaves = _tree_leaf_count(T)
    assert actual_leaves == expected_leaves, \
        f"Leaves {actual_leaves} != expected {expected_leaves}"

    dt = time.time() - t0
    print(f"\r[OK]  tree structure (depth={actual_depth}, leaves={actual_leaves}) — {dt:.1f}s")


def test_B_fft_consistency():
    """B_fft matches fft of basis polynomials."""
    print("       B_fft consistency running...", end="", flush=True)
    t0 = time.time()

    n = 512
    sk, _ = keygen(n)

    # Recompute B_fft from (f, g, F, G)
    neg_f = [-c for c in sk.f]
    neg_F = [-c for c in sk.F]
    B = [[sk.g, neg_f], [sk.G, neg_F]]
    B_fft_expected = [[fft([float(c) for c in poly]) for poly in row] for row in B]

    for i in range(2):
        for j in range(2):
            err = max(abs(a - b) for a, b in zip(sk.B_fft[i][j], B_fft_expected[i][j]))
            assert err < 1e-6, f"B_fft[{i}][{j}] mismatch: err={err}"

    dt = time.time() - t0
    print(f"\r[OK]  B_fft consistency — {dt:.1f}s")


def test_serialize_deserialize():
    """Round-trip serialization of public key."""
    print("       serialize/deserialize running...", end="", flush=True)
    t0 = time.time()

    n = 512
    sk, pk = keygen(n)

    data = serialize_public_key(pk.h, n)

    # Check expected size: 1 header + ceil(14*512/8) = 1 + 896 = 897
    expected_len = FALCON_PARAMS[n]["pk_bytelen"]
    assert len(data) == expected_len, \
        f"Serialized length {len(data)} != expected {expected_len}"

    # Header byte check
    assert data[0] == 0x09, f"Header byte {data[0]:02x} != 0x09 (logn=9)"

    # Deserialize and verify
    h_recovered = deserialize_public_key(data, n)
    assert h_recovered == pk.h, "Round-trip failed: h != deserialize(serialize(h))"

    dt = time.time() - t0
    print(f"\r[OK]  serialize/deserialize ({len(data)} bytes) — {dt:.1f}s")


def test_dynamic_mode():
    """Dynamic mode: sk has no B_fft or T until expanded."""
    print("       dynamic mode running...", end="", flush=True)
    t0 = time.time()

    n = 512
    sk, pk = keygen(n, mode="dynamic")

    # Dynamic key has no precomputed data
    assert sk.mode == "dynamic"
    assert sk.B_fft is None
    assert sk.T is None
    assert len(sk.f) == n
    assert len(pk.h) == n

    # h*f ≡ g still holds
    h_ntt = ntt(pk.h)
    f_ntt = ntt([c % Q for c in sk.f])
    product_ntt = ntt_mul(h_ntt, f_ntt)
    product = intt(product_ntt)
    g_mod = [c % Q for c in sk.g]
    for i in range(n):
        assert product[i] == g_mod[i]

    dt = time.time() - t0
    print(f"\r[OK]  dynamic mode (n=512) — {dt:.1f}s")


def test_expand_secret_key():
    """expand_secret_key converts dynamic → tree mode."""
    print("       expand_secret_key running...", end="", flush=True)
    t0 = time.time()

    n = 512
    sk_dyn, _ = keygen(n, mode="dynamic")
    assert sk_dyn.B_fft is None

    # Expand
    sk_tree = expand_secret_key(sk_dyn)
    assert sk_tree.mode == "tree"
    assert sk_tree.B_fft is not None
    assert sk_tree.T is not None
    assert len(sk_tree.B_fft) == 2
    assert len(sk_tree.B_fft[0][0]) == n

    # Tree has correct structure
    expected_depth = int(math.log2(n))
    assert _tree_depth(sk_tree.T) == expected_depth
    assert _tree_leaf_count(sk_tree.T) == n

    # Polynomials unchanged
    assert sk_tree.f == sk_dyn.f
    assert sk_tree.g == sk_dyn.g

    dt = time.time() - t0
    print(f"\r[OK]  expand_secret_key — {dt:.1f}s")


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    test_poly_div_ntt()
    test_keygen_structure()
    test_public_key_h_f_eq_g()
    test_ntru_equation()
    test_tree_structure()
    test_B_fft_consistency()
    test_serialize_deserialize()
    test_dynamic_mode()
    test_expand_secret_key()
    print(f"\nAll tests passed.")
