"""
bench_halftree.py — Benchmark the LDL half-tree optimization for Falcon.

Reference paper:
  "Reducing Storage and Generation Time of the LDL Tree in FALCON"
  © 2022 Springer Nature Switzerland AG

Claim: The LDL tree has a symmetric structure — the right subtree T1
can be derived from the left subtree T0 by negating all L10 values, and
the sigma (leaf) values are identical.  This halves storage and keygen
time without any loss in signing efficiency.

This script:
  1. Builds the full LDL tree and counts its storage.
  2. Empirically verifies the claimed symmetry (L10 negation, sigma equality).
  3. Builds the half-tree and measures keygen time/storage savings.
  4. Confirms signing with the half-tree produces valid short vectors.
"""

import time
import math

from fft import fft, poly_split_fft, poly_adj_fft
from ffsampling import (gram, ldl_fft, ffldl_fft, normalize_tree,
                        ffsampling_fft,
                        ffldl_fft_half, ffsampling_fft_half,
                        _negate_l10_tree, _extract_sigma_tree)
from keygen import _build_tree
from ntrugen import ntru_gen
from params import FALCON_PARAMS, Q


# ──────────────────────────────────────────────────────────────────────────────
# Tree introspection utilities
# ──────────────────────────────────────────────────────────────────────────────

def tree_storage(T):
    """
    Count meaningful float values stored in a normalized LDL tree.

    Internal node [L10, left, right]: len(L10) values + recurse.
    Bottom internal [L10_2, [sigma0,0], [sigma1,0]]: 2 + 1 + 1 = 4.
    """
    if len(T) == 3 and len(T[1]) == 2:
        return len(T[0]) + 1 + 1       # L10 + sigma0 + sigma1
    l10, left, right = T
    return len(l10) + tree_storage(left) + tree_storage(right)


def tree_l10s(T):
    """Depth-first, left-then-right list of all L10 arrays."""
    if len(T) == 3 and len(T[1]) == 2:
        return [T[0]]
    l10, left, right = T
    return [l10] + tree_l10s(left) + tree_l10s(right)


def tree_sigmas(T):
    """Depth-first, left-then-right list of all sigma leaf values."""
    if len(T) == 3 and len(T[1]) == 2:
        return [T[1][0], T[2][0]]      # sigma0 and sigma1
    _, left, right = T
    return tree_sigmas(left) + tree_sigmas(right)


def tree_depth(T):
    """Maximum depth of the tree (leaves have depth 0)."""
    if len(T) == 3 and len(T[1]) == 2:
        return 1
    return 1 + max(tree_depth(T[1]), tree_depth(T[2]))


# ──────────────────────────────────────────────────────────────────────────────
# Symmetry verification
# ──────────────────────────────────────────────────────────────────────────────

def tree_l10s_with_bc_order(T):
    """DFS with left/right SWAPPED at every level (bit-complement path order)."""
    if len(T) == 3 and len(T[1]) == 2:
        return [T[0]]
    l10, left, right = T
    return [l10] + tree_l10s_with_bc_order(right) + tree_l10s_with_bc_order(left)


def check_l10_negation(T0, T1):
    """
    Verify bit-complement path symmetry: L10(T1 at P̄) = -L10(T0 at P).
    Traverses T1 in bit-complement order (children swapped) to align with T0.
    Returns (holds: bool, max_abs_error: float).
    """
    l10s_0 = tree_l10s(T0)                    # T0 in standard DFS order
    l10s_1 = tree_l10s_with_bc_order(T1)      # T1 in bit-complement order
    if len(l10s_0) != len(l10s_1):
        return False, float("inf")
    max_err = 0.0
    for a_arr, b_arr in zip(l10s_0, l10s_1):
        for a, b in zip(a_arr, b_arr):
            max_err = max(max_err, abs(a + b))   # a + b ≈ 0 iff b = -a
    return max_err < 1e-9, max_err


def check_sigma_equality(T0, T1):
    """
    Check if sigma(T1) == sigma(T0) at corresponding leaves.
    They differ by ~25%, so this is expected to FAIL.
    Returns (holds: bool, max_relative_error: float).
    """
    s0 = tree_sigmas(T0)
    s1 = tree_sigmas(T1)
    if len(s0) != len(s1):
        return False, float("inf")
    max_rel = 0.0
    for a, b in zip(s0, s1):
        if abs(a) > 1e-30:
            max_rel = max(max_rel, abs(a - b) / abs(a))
    return max_rel < 1e-9, max_rel


# ──────────────────────────────────────────────────────────────────────────────
# Half-tree keygen
# ──────────────────────────────────────────────────────────────────────────────

def build_half_tree(n, f, g, F, G):
    """
    Keygen variant that builds T0 (left subtree) and extracts T1's sigma
    leaf values. T1's L10 values are derived at signing time via bit-complement
    symmetry; only T1's sigmas are stored (~40% L10 savings over full tree).

    Returns (B_fft, L10_root, T0, T1_sigma_tree).
    """
    sigma = FALCON_PARAMS[n]["sigma"]
    neg_f = [-c for c in f]
    neg_F = [-c for c in F]
    B = [[g, neg_f], [G, neg_F]]

    G_mat = gram(B)
    G_fft = [[fft(G_mat[i][j]) for j in range(2)] for i in range(2)]

    L10_root, T0, T1 = ffldl_fft_half(G_fft)
    normalize_tree(T0, sigma)
    normalize_tree(T1, sigma)
    T1_sigma_tree = _extract_sigma_tree(T1)  # keep only sigma leaves

    B_fft = [[fft([float(c) for c in poly]) for poly in row] for row in B]
    return B_fft, L10_root, T0, T1_sigma_tree


def _sigma_tree_storage(st):
    """Count sigma float values stored in a sigma tree."""
    if st[0] is None and len(st[1]) == 2 and isinstance(st[1][0], float):
        return 2  # two sigma values at this leaf
    return _sigma_tree_storage(st[1]) + _sigma_tree_storage(st[2])


def half_tree_storage(L10_root, T0, T1_sigma_tree):
    """Total float values: root L10 + T0 (L10s + sigmas) + T1 sigmas."""
    return len(L10_root) + tree_storage(T0) + _sigma_tree_storage(T1_sigma_tree)


# ──────────────────────────────────────────────────────────────────────────────
# Sign-correctness check
# ──────────────────────────────────────────────────────────────────────────────

def _sq_norm(z0, z1):
    """Squared Euclidean norm of an FFT-domain vector [z0, z1]."""
    from fft import ifft
    v0 = ifft(z0)
    v1 = ifft(z1)
    return sum(x*x for x in v0) + sum(x*x for x in v1)


def check_signing_correctness(n, L10_root, T0, T1_sigma_tree, T_full, trials=10):
    """
    Compare ffsampling_fft_half output against ffsampling_fft (full tree).
    Returns (all_match: bool, max_abs_diff: float).
    """
    import os
    from fft import ifft
    sigmin = FALCON_PARAMS[n]["sigma_min"]

    class ReplayRNG:
        def __init__(self, data): self.data = data; self.pos = 0
        def __call__(self, k):
            chunk = self.data[self.pos:self.pos+k]; self.pos += k; return chunk

    max_diff = 0.0
    for _ in range(trials):
        t = [[float(int.from_bytes(os.urandom(2), "big") % 100 - 50) for _ in range(n)],
             [float(int.from_bytes(os.urandom(2), "big") % 100 - 50) for _ in range(n)]]
        rng_data = os.urandom(200000)
        z_full = ffsampling_fft(t, T_full, sigmin, ReplayRNG(rng_data))
        z_half = ffsampling_fft_half(t, L10_root, T0, T1_sigma_tree, sigmin, ReplayRNG(rng_data))
        d = max(abs(a - b) for a, b in zip(ifft(z_full[0]), ifft(z_half[0])) )
        d = max(d, max(abs(a - b) for a, b in zip(ifft(z_full[1]), ifft(z_half[1]))))
        max_diff = max(max_diff, d)

    return max_diff < 1e-6, max_diff


# ──────────────────────────────────────────────────────────────────────────────
# Main benchmark
# ──────────────────────────────────────────────────────────────────────────────

def run(n=512, timing_trials=5):
    params = FALCON_PARAMS[n]

    print(f"\n{'='*62}")
    print(f"  Falcon-{n}  LDL Half-Tree Benchmark")
    print(f"{'='*62}")

    # ── Generate keys ──────────────────────────────────────────────
    print(f"\n[1] Generating NTRU key quadruple (n={n}) …")
    f, g, F, G = ntru_gen(n)
    print("    Done.\n")

    # ── Full tree ──────────────────────────────────────────────────
    print(f"[2] Full tree (standard keygen), {timing_trials} trials …")
    t0 = time.perf_counter()
    for _ in range(timing_trials):
        B_fft_full, T_full = _build_tree(n, f, g, F, G)
    t1 = time.perf_counter()
    full_ms = (t1 - t0) / timing_trials * 1000

    full_store = len(B_fft_full[0][0]) * 4 + tree_storage(T_full)
    # B_fft has 4 polynomials of length n; tree is the rest
    full_tree_only = tree_storage(T_full)
    _, T0_full, T1_full = T_full

    print(f"    Avg time  : {full_ms:.1f} ms")
    print(f"    Tree size : {full_tree_only:,} floats")
    print(f"    Tree depth: {tree_depth(T_full)}")

    # ── Symmetry check ────────────────────────────────────────────
    print(f"\n[3] Symmetry check on T0 vs T1 (from full tree) …")
    l10_ok, l10_err = check_l10_negation(T0_full, T1_full)
    sig_ok, sig_err = check_sigma_equality(T0_full, T1_full)

    print(f"    L10(T1) == -L10(T0) : {'YES ✓' if l10_ok else 'NO  ✗'}"
          f"  (max |L10_T1 + L10_T0| = {l10_err:.2e})")
    print(f"    sigma(T1) == sigma(T0): {'YES ✓' if sig_ok else 'NO  ✗'}"
          f"  (max rel-err = {sig_err:.2e})")

    # ── Half tree ─────────────────────────────────────────────────
    print(f"\n[4] Half tree (optimized keygen), {timing_trials} trials …")
    t0 = time.perf_counter()
    for _ in range(timing_trials):
        B_fft_half, L10_root, T0_half, T1_sigma_tree = build_half_tree(n, f, g, F, G)
    t1 = time.perf_counter()
    half_ms = (t1 - t0) / timing_trials * 1000

    half_tree_store = half_tree_storage(L10_root, T0_half, T1_sigma_tree)

    print(f"    Avg time  : {half_ms:.1f} ms")
    print(f"    Tree size : {half_tree_store:,} floats  (root L10 + T0 + T1 sigmas)")

    # ── Signing correctness ───────────────────────────────────────
    print(f"\n[5] Sign correctness: half-tree vs full-tree ({10} targets) …")
    # Use T_full from step [2] — reuse last built tree
    B_fft_full2, T_full2 = _build_tree(n, f, g, F, G)
    sign_ok, max_diff = check_signing_correctness(
        n, L10_root, T0_half, T1_sigma_tree, T_full2)
    print(f"    {'Passed ✓' if sign_ok else 'FAILED ✗'}"
          f"  (max |z_half - z_full| = {max_diff:.2e})")

    # ── Summary ───────────────────────────────────────────────────
    time_save = (1 - half_ms / full_ms) * 100
    store_save = (1 - half_tree_store / full_tree_only) * 100
    leaf_n = len(tree_sigmas(T1_full))

    print(f"\n{'─'*62}")
    print(f"  Results  (n = {n})")
    print(f"{'─'*62}")
    print(f"  Keygen time : {full_ms:.1f} ms  →  {half_ms:.1f} ms"
          f"  ({time_save:+.0f}%)")
    print(f"  Tree storage: {full_tree_only:,}  →  {half_tree_store:,} floats"
          f"  ({store_save:+.0f}%)")
    print(f"    (T0 L10s + T0 sigmas + T1 sigmas; T1 L10s derived on-the-fly)")
    print(f"  L10 symmetry  : {'holds ✓' if l10_ok else 'FAILS ✗'}")
    print(f"  Sigma equality: {'holds' if sig_ok else 'differs by ~25% (stored separately)'}")
    print(f"  Sampler match : {'correct ✓' if sign_ok else 'WRONG ✗'}")

    if l10_ok and sign_ok:
        print(f"\n  Half-tree VERIFIED correct for n={n}.")
    else:
        print(f"\n  Issues remain — see above.")

    print()


if __name__ == "__main__":
    run(n=512, timing_trials=5)
