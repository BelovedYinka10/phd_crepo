"""
test_fft.py — Comprehensive tests for the iterative FFT module.

Tests:
  1. Round-trip:  ifft(fft(f)) ≈ f
  2. Polynomial multiplication vs schoolbook (ground truth)
  3. Split/merge round-trip
  4. Adjoint properties
  5. Division: f/g * g ≈ f
  6. muladj consistency
  7. Autoadj operations
  8. invnorm2
  9. LDL decomposition
"""

import random
from fft import (
    fft, ifft,
    poly_add_fft, poly_sub_fft, poly_neg_fft,
    poly_mul_fft, poly_adj_fft, poly_div_fft,
    poly_mulconst_fft, poly_mulselfadj_fft, poly_muladj_fft,
    poly_mul_autoadj_fft, poly_div_autoadj_fft,
    poly_invnorm2_fft,
    poly_split_fft, poly_merge_fft,
    poly_LDL_fft,
)


def max_err(a, b):
    return max(abs(x - y) for x, y in zip(a, b))


def negacyclic_mul(f, g):
    """Schoolbook f·g mod (x^n+1) — ground truth."""
    n = len(f)
    h = [0.0] * n
    for i in range(n):
        for j in range(n):
            idx = (i + j) % n
            sign = -1 if (i + j) >= n else 1
            h[idx] += sign * f[i] * g[j]
    return h


# ── Tests ────────────────────────────────────────────────────────────────────

def test_round_trip(n):
    f = [random.gauss(0, 10) for _ in range(n)]
    err = max_err(f, ifft(fft(f)))
    assert err < 1e-6, f"round-trip n={n}: err={err}"


def test_mul(n):
    """FFT-based multiplication matches schoolbook."""
    f = [random.gauss(0, 10) for _ in range(n)]
    g = [random.gauss(0, 10) for _ in range(n)]
    h_fft = ifft(poly_mul_fft(fft(f), fft(g)))
    h_naive = negacyclic_mul(f, g)
    err = max_err(h_naive, h_fft)
    # Tolerance scales with n (accumulated floating-point error)
    tol = n * 1e-6
    assert err < tol, f"mul n={n}: err={err} (tol={tol})"


def test_split_merge(n):
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    f0, f1 = poly_split_fft(f_fft)
    err = max_err(f_fft, poly_merge_fft(f0, f1))
    assert err < 1e-6, f"split/merge n={n}: err={err}"


def test_adj(n):
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    # adj(adj(f)) = f
    err = max_err(f_fft, poly_adj_fft(poly_adj_fft(f_fft)))
    assert err < 1e-10, f"adj n={n}: err={err}"
    # f * adj(f) is real (imag = 0)
    hn = n >> 1
    sa = poly_mulselfadj_fft(f_fft)
    max_imag = max(abs(sa[i]) for i in range(hn, n))
    assert max_imag < 1e-6, f"mulselfadj imag n={n}: {max_imag}"


def test_div(n):
    """f / g * g ≈ f"""
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    g_fft = fft([random.gauss(0, 10) for _ in range(n)])
    recovered = poly_mul_fft(poly_div_fft(f_fft, g_fft), g_fft)
    err = max_err(f_fft, recovered)
    assert err < 1e-6, f"div n={n}: err={err}"


def test_muladj(n):
    """muladj(a,b) = mul(a, adj(b))"""
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    g_fft = fft([random.gauss(0, 10) for _ in range(n)])
    r1 = poly_muladj_fft(f_fft, g_fft)
    r2 = poly_mul_fft(f_fft, poly_adj_fft(g_fft))
    err = max_err(r1, r2)
    assert err < 1e-6, f"muladj n={n}: err={err}"


def test_autoadj(n):
    """mul then div by auto-adjoint b recovers original."""
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    b = poly_mulselfadj_fft(f_fft)  # auto-adjoint (real)
    recovered = poly_div_autoadj_fft(poly_mul_autoadj_fft(f_fft, b), b)
    err = max_err(f_fft, recovered)
    assert err < 1e-6, f"autoadj n={n}: err={err}"


def test_invnorm2(n):
    """invnorm2(a,b) · (|a|² + |b|²) ≈ 1"""
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    g_fft = fft([random.gauss(0, 10) for _ in range(n)])
    inv = poly_invnorm2_fft(f_fft, g_fft)
    nf = poly_mulselfadj_fft(f_fft)
    ng = poly_mulselfadj_fft(g_fft)
    hn = n >> 1
    for u in range(hn):
        p = inv[u] * (nf[u] + ng[u])
        assert abs(p - 1.0) < 1e-6, f"invnorm2 n={n} u={u}: {p}"


def test_ldl(n):
    """LDL: D11 + |L10|²·D00 ≈ original G11."""
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    g_fft = fft([random.gauss(0, 10) for _ in range(n)])
    g00 = poly_add_fft(poly_mulselfadj_fft(f_fft), poly_mulselfadj_fft(g_fft))
    g01 = poly_muladj_fft(f_fft, g_fft)
    g11 = list(g00)
    d00, l10, d11 = poly_LDL_fft(list(g00), list(g01), list(g11))
    hn = n >> 1
    for u in range(hn):
        l_norm2 = l10[u] ** 2 + l10[u + hn] ** 2
        g11_recon = d11[u] + l_norm2 * d00[u]
        assert abs(g11_recon - g00[u]) < 1e-3, \
            f"LDL n={n} u={u}: err={abs(g11_recon - g00[u])}"


def test_add_sub_neg(n):
    """add/sub/neg consistency."""
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    g_fft = fft([random.gauss(0, 10) for _ in range(n)])
    # f + g - g ≈ f
    err = max_err(f_fft, poly_sub_fft(poly_add_fft(f_fft, g_fft), g_fft))
    assert err < 1e-10, f"add/sub n={n}: err={err}"
    # f + neg(f) ≈ 0
    zero = poly_add_fft(f_fft, poly_neg_fft(f_fft))
    err = max(abs(x) for x in zero)
    assert err < 1e-10, f"neg n={n}: err={err}"


def test_mulconst(n):
    """mulconst(f, c) = c * f."""
    f_fft = fft([random.gauss(0, 10) for _ in range(n)])
    c = 3.14159
    scaled = poly_mulconst_fft(f_fft, c)
    expected = [x * c for x in f_fft]
    err = max_err(scaled, expected)
    assert err < 1e-10, f"mulconst n={n}: err={err}"


# ── Main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    random.seed(42)

    tests = [
        ("round-trip",   test_round_trip),
        ("mul",          test_mul),
        ("split/merge",  test_split_merge),
        ("adj",          test_adj),
        ("div",          test_div),
        ("muladj",       test_muladj),
        ("autoadj",      test_autoadj),
        ("invnorm2",     test_invnorm2),
        ("LDL",          test_ldl),
        ("add/sub/neg",  test_add_sub_neg),
        ("mulconst",     test_mulconst),
    ]

    sizes = [8, 16, 64, 512, 1024]

    for name, fn in tests:
        for n in sizes:
            if name == "split/merge" and n < 8:
                continue
            fn(n)
        print(f"[OK]  {name}")

    print(f"\nAll {len(tests)} tests passed for n ∈ {sizes}.")
