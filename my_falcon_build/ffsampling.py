"""
ffsampling.py — Fast Fourier Sampling for Falcon
==================================================

Implements the LDL tree construction and Gaussian sampling used in
Falcon's signature generation (Algorithms 8–11 of the specification).

Pipeline
--------
  keygen:  (f,g,F,G) → B = [[g,-f],[G,-F]] → gram(B) → ffldl_fft → normalize_tree → T
  sign:    (t, T) → ffsampling_fft → short vector z

Data structures
---------------
  LDL tree (nested lists):
    Internal node: [L10, T_left, T_right]   — len == 3
    Leaf (after normalize): [sigma, 0]       — len == 2

  All FFT-domain arrays use split-halves layout:
    a[0..n/2-1] = real parts,  a[n/2..n-1] = imaginary parts

Reference
---------
  Falcon spec v1.2, Algorithms 8–11
  Reference Python implementation by Thomas Prest (PQShield)
"""

__all__ = [
    "gram", "ldl_fft", "ffldl_fft", "normalize_tree",
    "ffnp_fft", "ffsampling_fft",
    # Half-tree optimization (paper: Springer 2022)
    "ffldl_fft_half", "ffsampling_fft_half",
    "_negate_l10_tree", "_extract_sigma_tree",
]

from math import sqrt
from fft import (fft, ifft,
                 poly_add_fft, poly_sub_fft, poly_mul_fft,
                 poly_div_fft, poly_adj_fft,
                 poly_split_fft, poly_merge_fft,
                 poly_mul_coeff, poly_adj_coeff, poly_add_coeff)
from samplerz import samplerz


# ──────────────────────────────────────────────────────────────────────────────
# Gram matrix  (coefficient domain)
# ──────────────────────────────────────────────────────────────────────────────

def gram(B):
    """
    Compute the Gram matrix G = B^* · B for a 2×2 polynomial matrix B.

    B is in coefficient representation: B[i][j] is a list of n floats/ints.
    G[i][j] = sum_k  B[i][k] · adj(B[j][k])

    Returns G as a 2×2 matrix of coefficient-domain polynomials.
    """
    dim = len(B)
    ncols = len(B[0])
    deg = len(B[0][0])
    G = [[[0] * deg for _ in range(dim)] for _ in range(dim)]
    for i in range(dim):
        for j in range(dim):
            for k in range(ncols):
                G[i][j] = poly_add_coeff(G[i][j],
                                         poly_mul_coeff(B[i][k],
                                                        poly_adj_coeff(B[j][k])))
    return G


# ──────────────────────────────────────────────────────────────────────────────
# LDL decomposition  (FFT domain)
# ──────────────────────────────────────────────────────────────────────────────

def ldl_fft(G):
    """
    LDL* decomposition of a 2×2 Gram matrix in FFT domain.

    G = L · D · L*  where L is lower-triangular with 1s on diagonal,
    D is diagonal.

    Returns [L, D] as 2×2 matrices (FFT-domain arrays).
    """
    deg = len(G[0][0])

    zero = [0.0] * deg
    one = [0.0] * deg
    hn = deg >> 1
    # Identity in FFT domain: [1, 1, ..., 1, 0, 0, ..., 0]
    # Actually, fft([1, 0, 0, ...]) = all-ones in split-halves.
    # Re parts = 1.0, Im parts = 0.0.
    for i in range(hn):
        one[i] = 1.0

    D00 = G[0][0][:]
    L10 = poly_div_fft(G[1][0], G[0][0])
    D11 = poly_sub_fft(G[1][1],
                       poly_mul_fft(poly_mul_fft(L10, poly_adj_fft(L10)),
                                    G[0][0]))
    L = [[one[:], zero[:]], [L10, one[:]]]
    D = [[D00, zero[:]], [zero[:], D11]]
    return [L, D]


# ──────────────────────────────────────────────────────────────────────────────
# ffLDL tree construction  (Algorithm 9)
# ──────────────────────────────────────────────────────────────────────────────

def ffldl_fft(G):
    """
    Build the ffLDL decomposition tree from Gram matrix G (FFT domain).

    Recursively bisects the diagonal of D via poly_split_fft until
    reaching the leaf level (array length 2 = one complex value).

    Returns:
        Internal node: [L10, T_left, T_right]  (len == 3)
        Leaf:          [L10, D00, D11]          (len == 3, but D00/D11 have len 2)
    """
    n = len(G[0][0])
    L, D = ldl_fft(G)

    if n > 2:
        # Bisect: split the diagonal elements of D
        d00, d01 = poly_split_fft(D[0][0])
        d10, d11 = poly_split_fft(D[1][1])
        # Form sub-Gram matrices (Hermitian structure)
        G0 = [[d00, d01], [poly_adj_fft(d01), d00]]
        G1 = [[d10, d11], [poly_adj_fft(d11), d10]]
        return [L[1][0], ffldl_fft(G0), ffldl_fft(G1)]
    elif n == 2:
        # Leaf: D00 and D11 are 2-element arrays [re, im]
        return [L[1][0], D[0][0], D[1][1]]


# ──────────────────────────────────────────────────────────────────────────────
# Tree normalization
# ──────────────────────────────────────────────────────────────────────────────

def normalize_tree(tree, sigma):
    """
    Convert LDL tree leaves from ||b_i||^2 to sigma / ||b_i|| (in-place).

    After normalization, each leaf stores the per-coefficient sigma
    used by samplerz during signing.

    Internal nodes (len == 3): recurse on tree[1] and tree[2].
    Leaves (len == 2): tree[0] = sigma / sqrt(tree[0]), tree[1] = 0.
    """
    if len(tree) == 3:
        normalize_tree(tree[1], sigma)
        normalize_tree(tree[2], sigma)
    else:
        # Leaf: tree = [re, im] where re = ||b_i||^2, im ≈ 0
        tree[0] = sigma / sqrt(tree[0])
        tree[1] = 0


# ──────────────────────────────────────────────────────────────────────────────
# Fast Fourier Nearest Plane  (Algorithm 10, rounding version)
# ──────────────────────────────────────────────────────────────────────────────

def ffnp_fft(t, T):
    """
    Fast Fourier Nearest Plane reduction of target vector t using tree T.

    Deterministic rounding (no randomness). Used for testing and
    as a simpler precursor to ffsampling_fft.

    Format: FFT (split-halves layout).

    Base case: when T's children are leaves (len==2), we're at the
    bottom internal node. Process t[1] and t[0] directly using L10,
    then round the real parts.
    """
    n = len(t[0])
    z = [None, None]

    if len(T) == 3 and len(T[1]) == 2:
        # Bottom internal node: T = [L10, D00_leaf, D11_leaf]
        # T[1] and T[2] are [sigma, 0] leaves (or pre-normalize [re, im])
        # t[0] and t[1] are 2-element FFT arrays [re, im]
        l10 = T[0]
        # z[1] = round(Re(t[1]))
        z1_re = round(t[1][0])
        z[1] = [float(z1_re), 0.0]
        # t0' = t[0] + (t[1] - z[1]) * l10
        # All are 2-element arrays: [re, im]
        diff_re = t[1][0] - z[1][0]
        diff_im = t[1][1] - z[1][1]
        l_re = l10[0]
        l_im = l10[1]
        adj_re = diff_re * l_re - diff_im * l_im
        adj_im = diff_re * l_im + diff_im * l_re
        t0b_re = t[0][0] + adj_re
        t0b_im = t[0][1] + adj_im
        z0_re = round(t0b_re)
        z[0] = [float(z0_re), 0.0]
        return z
    else:
        l10, T0, T1 = T
        z[1] = poly_merge_fft(*ffnp_fft(poly_split_fft(t[1]), T1))
        t0b = poly_add_fft(t[0],
                           poly_mul_fft(poly_sub_fft(t[1], z[1]), l10))
        z[0] = poly_merge_fft(*ffnp_fft(poly_split_fft(t0b), T0))
        return z


# ──────────────────────────────────────────────────────────────────────────────
# Fast Fourier Sampling  (Algorithm 11)
# ──────────────────────────────────────────────────────────────────────────────

def ffsampling_fft(t, T, sigmin, randombytes):
    """
    Gaussian sampling over the LDL tree (Algorithm 11, ffSampling).

    At each internal node:
      1. Split t[1], recurse into T1 to sample z[1]
      2. Adjust t[0] using L10 and the residual (t[1] - z[1])
      3. Split adjusted t[0], recurse into T0 to sample z[0]
      4. Merge both results

    At leaves: call samplerz() with mu = Re(t[i]) and sigma from tree.

    Args:
        t: target vector [t0, t1], each an FFT-domain array
        T: normalized LDL tree
        sigmin: minimum sigma (from Falcon parameters)
        randombytes: entropy source (os.urandom or ChaCha20)

    Returns:
        z = [z0, z1] sampled vector in FFT domain
    """
    n = len(t[0])
    z = [None, None]

    if len(T) == 3 and len(T[1]) == 2:
        # Bottom internal node: T = [L10, [sigma0, 0], [sigma1, 0]]
        # t[0], t[1] are 2-element FFT arrays [re, im]
        l10 = T[0]
        sigma1 = T[2][0]   # sigma for z[1]
        sigma0 = T[1][0]   # sigma for z[0]

        # Sample z[1]
        z1_re = samplerz(t[1][0], sigma1, sigmin, randombytes)
        z[1] = [float(z1_re), 0.0]

        # Adjust t[0]: t0' = t[0] + (t[1] - z[1]) * L10
        diff_re = t[1][0] - z[1][0]
        diff_im = t[1][1] - z[1][1]
        l_re = l10[0]
        l_im = l10[1]
        adj_re = diff_re * l_re - diff_im * l_im
        adj_im = diff_re * l_im + diff_im * l_re
        t0b_re = t[0][0] + adj_re
        t0b_im = t[0][1] + adj_im

        # Sample z[0]
        z0_re = samplerz(t0b_re, sigma0, sigmin, randombytes)
        z[0] = [float(z0_re), 0.0]
        return z
    else:
        l10, T0, T1 = T
        z[1] = poly_merge_fft(
            *ffsampling_fft(poly_split_fft(t[1]), T1, sigmin, randombytes))
        t0b = poly_add_fft(t[0],
                           poly_mul_fft(poly_sub_fft(t[1], z[1]), l10))
        z[0] = poly_merge_fft(
            *ffsampling_fft(poly_split_fft(t0b), T0, sigmin, randombytes))
        return z


# ──────────────────────────────────────────────────────────────────────────────
# Half-tree optimization  (Springer 2022 paper)
# ──────────────────────────────────────────────────────────────────────────────
#
# Core algebraic fact (from NTRU equation fG − gF = q):
#   det(Gram) = D[0][0] · D[1][1] = q²  (pointwise in FFT domain)
#
# Bit-complement path symmetry: for every node at path P in T0, the node at
# the bit-complement path P̄ (every L↔R flipped) in T1 satisfies:
#   L10(T1 at P̄) = −L10(T0 at P)
#
# Consequence: T1's L10 values are derivable from T0 by negation + child-swap.
# However, T1's leaf sigma values differ from T0's by ~25% and must be stored
# or recomputed separately.  Storage reduction for L10 values alone is ~40%.
# ──────────────────────────────────────────────────────────────────────────────

def ffldl_fft_half(G):
    """
    Build the left subtree T0 and right subtree T1 of the ffLDL tree.

    Returns (L10_root, T0, T1) where:
      - L10_root is the n-element L10 array at the root level
      - T0 is the full left subtree (ffldl_fft of D[0][0])
      - T1 is the full right subtree (ffldl_fft of D[1][1])

    T1 is returned so the caller can extract sigma leaf values via
    _extract_sigma_tree(T1) after normalization, then discard T1's L10s.
    Call normalize_tree(T0, sigma) and normalize_tree(T1, sigma) after this.
    """
    L, D = ldl_fft(G)
    d00, d01 = poly_split_fft(D[0][0])
    d10, d11 = poly_split_fft(D[1][1])
    G0 = [[d00, d01], [poly_adj_fft(d01), d00]]
    G1 = [[d10, d11], [poly_adj_fft(d11), d10]]
    return L[1][0], ffldl_fft(G0), ffldl_fft(G1)


def _extract_sigma_tree(T):
    """
    Extract only the sigma leaf structure from a normalized LDL tree.

    Returns a tree with the same shape as T but with None in place of every
    L10 array.  Only the leaf sigma values ([sigma, 0] pairs) are kept.
    Used to capture T1's sigmas after normalization, before discarding T1.
    """
    if len(T) == 3 and len(T[1]) == 2:
        return [None, T[1][:], T[2][:]]
    _, left, right = T
    return [None, _extract_sigma_tree(left), _extract_sigma_tree(right)]


def _negate_l10_tree(T0, T1_sigma_tree):
    """
    Reconstruct T1 from T0 using the bit-complement path symmetry.

    At every node: negate L10 and swap left/right children (so that
    T1_derived at path P carries −L10(T0 at P̄), matching T1_actual).
    Sigma leaf values come from T1_sigma_tree (T1's normalized sigmas),
    NOT from T0, because T1 and T0 sigmas differ by ~25%.

    Args:
        T0: normalized left subtree
        T1_sigma_tree: output of _extract_sigma_tree(T1) after normalization
    """
    if len(T0) == 3 and len(T0[1]) == 2:
        # Leaf: use T1's sigma values (T1_sigma_tree[1] and [2])
        return [[-x for x in T0[0]], T1_sigma_tree[1], T1_sigma_tree[2]]
    l10, T0_left, T0_right = T0
    T1_sig_left, T1_sig_right = T1_sigma_tree[1], T1_sigma_tree[2]
    # Swap T0's children (bit-complement), but keep T1's sigma tree aligned
    return [[-x for x in l10],
            _negate_l10_tree(T0_right, T1_sig_left),
            _negate_l10_tree(T0_left, T1_sig_right)]


def ffsampling_fft_half(t, L10_root, T0, T1_sigma_tree, sigmin, randombytes):
    """
    FFSampling using the half-tree (L10_root, T0, T1_sigma_tree).

    Reconstructs T1 on-the-fly from T0 via _negate_l10_tree, which applies
    the bit-complement path symmetry (negate + swap children) for L10 values,
    and pulls leaf sigmas from T1_sigma_tree.  Then calls the standard sampler.

    Args:
        t:             target vector [t0, t1] in FFT domain
        L10_root:      root-level L10 (n floats)
        T0:            normalized left subtree
        T1_sigma_tree: _extract_sigma_tree(T1) — T1's sigma leaf structure
        sigmin:        minimum sigma (Falcon parameter)
        randombytes:   entropy source
    """
    T1_derived = _negate_l10_tree(T0, T1_sigma_tree)
    return ffsampling_fft(t, [L10_root, T0, T1_derived], sigmin, randombytes)
