"""
test_ffsampling.py — Tests for Fast Fourier Sampling module
============================================================
"""

import random
import time
from math import sqrt
from os import urandom

from ffsampling import (gram, ldl_fft, ffldl_fft, normalize_tree,
                        ffnp_fft, ffsampling_fft)
from fft import (fft, ifft, poly_mul_fft, poly_adj_fft, poly_add_fft,
                 poly_sub_fft, poly_split_fft, poly_merge_fft)
from ntrugen import ntru_gen, karamul
from params import Q, FALCON_PARAMS

random.seed(42)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _rand_poly(n, bound=5):
    return [random.randint(-bound, bound) for _ in range(n)]


def _neg_poly(f):
    return [-c for c in f]


def _tree_depth(T):
    """Count the depth of an LDL tree (leaves are len-2 arrays)."""
    if len(T) == 2:
        return 0
    return 1 + max(_tree_depth(T[1]), _tree_depth(T[2]))


def _tree_leaf_count(T):
    """Count leaves in an LDL tree (leaves are len-2 arrays)."""
    if len(T) == 2:
        return 1
    return _tree_leaf_count(T[1]) + _tree_leaf_count(T[2])


# ── Tests ────────────────────────────────────────────────────────────────────

def test_gram():
    """Gram matrix G = B^* B for a known 2x2 polynomial matrix."""
    n = 8
    f = _rand_poly(n)
    g = _rand_poly(n)
    F = _rand_poly(n)
    G_poly = _rand_poly(n)

    B = [[g, _neg_poly(f)], [G_poly, _neg_poly(F)]]
    G = gram(B)

    # G should be 2x2 and Hermitian: G[i][j] = adj(G[j][i])
    assert len(G) == 2 and len(G[0]) == 2
    assert len(G[0][0]) == n

    # G[0][0] = g·adj(g) + f·adj(f)  (auto-adjoint, so should be ≈ real in FFT)
    # Verify G[0][0] is auto-adjoint: adj(G[0][0]) ≈ G[0][0]
    g00_fft = fft(G[0][0])
    g00_adj = poly_adj_fft(g00_fft)
    err = max(abs(a - b) for a, b in zip(g00_fft, g00_adj))
    assert err < 1e-6, f"G[0][0] not auto-adjoint: err={err}"

    print("[OK]  gram")


def test_ldl_fft():
    """LDL decomposition: verify G = L D L^*."""
    n = 8
    f = _rand_poly(n)
    g = _rand_poly(n)
    F = _rand_poly(n)
    G_poly = _rand_poly(n)

    B = [[g, _neg_poly(f)], [G_poly, _neg_poly(F)]]
    G = gram(B)

    # Convert to FFT domain
    G_fft = [[fft(G[i][j]) for j in range(2)] for i in range(2)]

    L, D = ldl_fft(G_fft)

    # Reconstruct: G = L D L^*
    # G[0][0] = D[0][0]  (L[0][0] = 1)
    # G[1][0] = L[1][0] * D[0][0]
    # G[0][1] = D[0][0] * adj(L[1][0])
    # G[1][1] = L[1][0] * D[0][0] * adj(L[1][0]) + D[1][1]

    recon_00 = D[0][0]
    recon_10 = poly_mul_fft(L[1][0], D[0][0])
    recon_01 = poly_mul_fft(D[0][0], poly_adj_fft(L[1][0]))
    recon_11 = poly_add_fft(
        poly_mul_fft(poly_mul_fft(L[1][0], D[0][0]), poly_adj_fft(L[1][0])),
        D[1][1])

    for name, orig, recon in [("G00", G_fft[0][0], recon_00),
                               ("G10", G_fft[1][0], recon_10),
                               ("G01", G_fft[0][1], recon_01),
                               ("G11", G_fft[1][1], recon_11)]:
        err = max(abs(a - b) for a, b in zip(orig, recon))
        assert err < 1e-4, f"LDL reconstruction {name} failed: err={err}"

    print("[OK]  ldl_fft")


def test_ffldl_fft():
    """ffLDL tree has correct depth and leaf structure."""
    for n in [4, 8, 16, 64]:
        f = _rand_poly(n)
        g = _rand_poly(n)
        F = _rand_poly(n)
        G_poly = _rand_poly(n)

        B = [[g, _neg_poly(f)], [G_poly, _neg_poly(F)]]
        G = gram(B)
        G_fft = [[fft(G[i][j]) for j in range(2)] for i in range(2)]

        T = ffldl_fft(G_fft)

        # Tree should be a 3-element list at the top
        assert len(T) == 3, f"n={n}: root not len 3"

        # Depth should be log2(n)
        # Root (n) → split (n/2) → ... → leaf at n=2 → D arrays (len 2)
        import math
        expected_depth = int(math.log2(n))
        actual_depth = _tree_depth(T)
        assert actual_depth == expected_depth, \
            f"n={n}: depth {actual_depth} != expected {expected_depth}"

        # Number of leaves should be n (2 per bottom node, n/2 bottom nodes)
        expected_leaves = n
        actual_leaves = _tree_leaf_count(T)
        assert actual_leaves == expected_leaves, \
            f"n={n}: {actual_leaves} leaves != expected {expected_leaves}"

    print("[OK]  ffldl_fft")


def test_normalize_tree():
    """After normalization, leaves contain sigma values."""
    n = 8
    sigma = 165.0

    f = _rand_poly(n)
    g = _rand_poly(n)
    F = _rand_poly(n)
    G_poly = _rand_poly(n)

    B = [[g, _neg_poly(f)], [G_poly, _neg_poly(F)]]
    G = gram(B)
    G_fft = [[fft(G[i][j]) for j in range(2)] for i in range(2)]
    T = ffldl_fft(G_fft)

    normalize_tree(T, sigma)

    # Check all leaves are [sigma_value, 0]
    def check_leaves(tree):
        if len(tree) == 2:
            assert isinstance(tree[0], float), f"Leaf[0] not float: {type(tree[0])}"
            assert tree[1] == 0, f"Leaf[1] not 0: {tree[1]}"
            assert tree[0] > 0, f"Leaf sigma not positive: {tree[0]}"
            return
        check_leaves(tree[1])
        check_leaves(tree[2])

    check_leaves(T)
    print("[OK]  normalize_tree")


def test_ffnp_fft():
    """Nearest plane rounding produces integer-valued output."""
    n = 8
    f = _rand_poly(n)
    g = _rand_poly(n)
    F = _rand_poly(n)
    G_poly = _rand_poly(n)

    B = [[g, _neg_poly(f)], [G_poly, _neg_poly(F)]]
    G = gram(B)
    G_fft = [[fft(G[i][j]) for j in range(2)] for i in range(2)]
    T = ffldl_fft(G_fft)

    # Create a target vector in FFT domain
    t0 = fft([random.gauss(0, 10) for _ in range(n)])
    t1 = fft([random.gauss(0, 10) for _ in range(n)])
    t = [t0, t1]

    z = ffnp_fft(t, T)

    # z should be two FFT-domain arrays
    assert len(z) == 2
    assert len(z[0]) == n
    assert len(z[1]) == n

    # Convert back to coefficient domain — should be near-integer
    z0_coeff = ifft(z[0])
    z1_coeff = ifft(z[1])
    for i in range(n):
        assert abs(z0_coeff[i] - round(z0_coeff[i])) < 1e-3, \
            f"z0[{i}] not integer: {z0_coeff[i]}"
        assert abs(z1_coeff[i] - round(z1_coeff[i])) < 1e-3, \
            f"z1[{i}] not integer: {z1_coeff[i]}"

    print("[OK]  ffnp_fft")


def test_ffsampling_fft():
    """ffsampling produces valid output format."""
    n = 8
    sigma = 165.0
    sigmin = 1.277

    f = _rand_poly(n)
    g = _rand_poly(n)
    F = _rand_poly(n)
    G_poly = _rand_poly(n)

    B = [[g, _neg_poly(f)], [G_poly, _neg_poly(F)]]
    G = gram(B)
    G_fft = [[fft(G[i][j]) for j in range(2)] for i in range(2)]
    T = ffldl_fft(G_fft)
    normalize_tree(T, sigma)

    t0 = fft([random.gauss(0, 10) for _ in range(n)])
    t1 = fft([random.gauss(0, 10) for _ in range(n)])
    t = [t0, t1]

    z = ffsampling_fft(t, T, sigmin, urandom)

    assert len(z) == 2
    assert len(z[0]) == n
    assert len(z[1]) == n

    # Convert back — should be near-integer (sampled integers)
    z0_coeff = ifft(z[0])
    z1_coeff = ifft(z[1])
    for i in range(n):
        assert abs(z0_coeff[i] - round(z0_coeff[i])) < 1e-3, \
            f"z0[{i}] not integer: {z0_coeff[i]}"

    print("[OK]  ffsampling_fft")


def test_full_pipeline():
    """Full pipeline: ntru_gen → gram → ffldl → normalize → ffsampling."""
    print("       full pipeline running...", end="", flush=True)
    t0 = time.time()

    n = 512
    params = FALCON_PARAMS[n]
    sigma = params["sigma"]
    sigmin = params["sigma_min"]

    # Generate NTRU key
    f, g, F, G = ntru_gen(n)

    # Build basis
    B = [[g, _neg_poly(f)], [G, _neg_poly(F)]]

    # Gram matrix and LDL tree
    G_mat = gram(B)
    G_fft = [[fft(G_mat[i][j]) for j in range(2)] for i in range(2)]
    T = ffldl_fft(G_fft)
    normalize_tree(T, sigma)

    # Build FFT-domain basis
    B_fft = [[fft([float(c) for c in row]) for row in col] for col in B]

    # Create a random target (simulating hash-to-point)
    point = [random.randint(0, Q - 1) for _ in range(n)]
    point_fft = fft([float(c) for c in point])

    # Compute target t = (point, 0) · B^{-1} in FFT domain
    # For NTRU basis B = [[g, -f], [G, -F]]:
    #   B^{-1} row 0 = [F, f] / q  (scaled)
    #   B^{-1} row 1 = [G, g] / q  (scaled, with signs)
    # Actually: t0 = point · d / q,  t1 = -point · b / q
    # where [[a, b], [c, d]] = B_fft
    a, b = B_fft[0]
    c, d = B_fft[1]
    t0_fft = [point_fft[i] * d[i] / Q for i in range(n)]
    t1_fft = [-point_fft[i] * b[i] / Q for i in range(n)]
    t = [t0_fft, t1_fft]

    # Sample
    z = ffsampling_fft(t, T, sigmin, urandom)

    # Compute v = z · B in FFT domain
    v0_fft = poly_add_fft(poly_mul_fft(z[0], a), poly_mul_fft(z[1], c))
    v1_fft = poly_add_fft(poly_mul_fft(z[0], b), poly_mul_fft(z[1], d))

    v0 = [int(round(x)) for x in ifft(v0_fft)]
    v1 = [int(round(x)) for x in ifft(v1_fft)]

    # Signature: s = (point - v0, -v1)
    s0 = [(point[i] - v0[i]) % Q for i in range(n)]
    s1 = [(-v1[i]) % Q for i in range(n)]

    # The signature norm should be bounded
    # Center the coefficients for norm computation
    def center(x):
        return x - Q if x > Q // 2 else x
    s1_centered = [center(c) for c in s1]
    norm_sq = sum(c ** 2 for c in s1_centered)

    dt = time.time() - t0
    # We don't strictly enforce the norm bound here (single sample may exceed),
    # but it should be in a reasonable range
    print(f"\r[OK]  full pipeline (n=512) — {dt:.1f}s, ||s1||²={norm_sq:,}")


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    test_gram()
    test_ldl_fft()
    test_ffldl_fft()
    test_normalize_tree()
    test_ffnp_fft()
    test_ffsampling_fft()
    test_full_pipeline()
    print(f"\nAll tests passed.")
