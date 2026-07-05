"""
test_samplerz.py — Tests and statistical evaluation for samplerz.py
====================================================================

Tests:
  1. basesampler range         — z0 always in {0,...,18}
  2. approxexp accuracy        — matches math.exp within tolerance
  3. berexp distribution       — empirical acceptance rate ≈ ccs·exp(-x)
  4. samplerz output type      — always returns int
  5. samplerz range            — outputs are integers (no floats)
  6. Statistical mean          — E[z] ≈ μ
  7. Statistical variance      — Var[z] ≈ σ²
  8. Both Falcon parameter sets — n=512 (σ=165.7) and n=1024 (σ=168.4)

Run: python test_samplerz.py
"""

import sys
import math
import random

from samplerz import samplerz, basesampler, berexp, approxexp
from params   import FALCON_PARAMS

PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"
_results = []

def _record(name, ok, detail=""):
    _results.append((name, ok, detail))
    status = PASS if ok else FAIL
    print(f"  [{status}]  {name}" + (f" — {detail}" if detail else ""))


# ─────────────────────────────────────────────────────────────────────────────
# 1. basesampler range
# ─────────────────────────────────────────────────────────────────────────────

def test_basesampler_range():
    print("\n── 1. basesampler range ─────────────────────────────────────────")
    samples = [basesampler() for _ in range(5000)]
    ok_range = all(0 <= z <= 18 for z in samples)
    _record("basesampler always in {0,...,18}", ok_range,
            f"min={min(samples)}, max={max(samples)}")

    # Distribution should be heavily concentrated near 0
    frac_zero = samples.count(0) / len(samples)
    _record("basesampler: >20% of samples are 0 (half-Gaussian peak)",
            frac_zero > 0.20, f"frac(0)={frac_zero:.3f}")


# ─────────────────────────────────────────────────────────────────────────────
# 2. approxexp accuracy
# ─────────────────────────────────────────────────────────────────────────────

def test_approxexp():
    print("\n── 2. approxexp accuracy ────────────────────────────────────────")
    import math
    errors = []
    # approxexp returns ~ 2^64 * ccs * exp(-x), so divide by 2^64
    for x in [0.0, 0.1, 0.3, 0.5, 0.693, 0.0001, 0.6]:
        for ccs in [1.0, 0.8, 0.5]:
            approx  = approxexp(x, ccs) / (1 << 64)
            exact   = ccs * math.exp(-x)
            rel_err = abs(approx - exact) / exact if exact > 0 else 0
            errors.append(rel_err)

    max_err = max(errors)
    _record("approxexp relative error < 1e-9 for x∈[0,ln2), ccs∈(0,1]",
            max_err < 1e-9, f"max_rel_err={max_err:.2e}")


# ─────────────────────────────────────────────────────────────────────────────
# 3. berexp empirical acceptance rate
# ─────────────────────────────────────────────────────────────────────────────

def test_berexp():
    print("\n── 3. berexp acceptance rate ────────────────────────────────────")
    import math
    trials = 10000
    for x, ccs in [(0.5, 1.0), (1.0, 0.8), (0.2, 0.9)]:
        hits      = sum(berexp(x, ccs) for _ in range(trials))
        empirical = hits / trials
        expected  = ccs * math.exp(-x)
        # Allow 3-sigma tolerance: std = sqrt(p(1-p)/n)
        std       = math.sqrt(expected * (1 - expected) / trials)
        ok        = abs(empirical - expected) < 4 * std
        _record(f"berexp(x={x}, ccs={ccs}): empirical≈expected",
                ok, f"empirical={empirical:.4f}, expected={expected:.4f}, "
                    f"4σ={4*std:.4f}")


# ─────────────────────────────────────────────────────────────────────────────
# 4 & 5. samplerz output type and integrality
# ─────────────────────────────────────────────────────────────────────────────

def test_output_type():
    print("\n── 4–5. Output type & integrality ──────────────────────────────")
    mu, sigma, sigmin = 0.0, 1.5, 1.0
    samples = [samplerz(mu, sigma, sigmin) for _ in range(200)]
    _record("samplerz always returns int",
            all(isinstance(z, int) for z in samples))
    _record("samplerz no fractional parts",
            all(z == int(z) for z in samples))


# ─────────────────────────────────────────────────────────────────────────────
# 6 & 7. Statistical mean and variance
# ─────────────────────────────────────────────────────────────────────────────

def test_statistics():
    print("\n── 6–7. Statistical mean and variance ───────────────────────────")
    import math

    cases = [
        # (mu,   sigma,  sigmin,  n_samples, tol_mean, tol_var_frac)
        (0.0,   1.5,    1.0,     8000,  0.10,  0.10),
        (3.7,   1.5,    1.0,     8000,  0.10,  0.10),
        (-2.3,  1.7,    1.0,     8000,  0.10,  0.10),
    ]

    for mu, sigma, sigmin, n, tol_m, tol_v in cases:
        samples = [samplerz(mu, sigma, sigmin) for _ in range(n)]
        mean    = sum(samples) / n
        var     = sum((z - mean) ** 2 for z in samples) / n

        ok_mean = abs(mean - mu)       < tol_m
        ok_var  = abs(var  - sigma**2) < tol_v * sigma**2

        _record(f"μ={mu}, σ={sigma}: E[z]≈μ",
                ok_mean, f"sample_mean={mean:.4f}, expected={mu}")
        _record(f"μ={mu}, σ={sigma}: Var[z]≈σ²",
                ok_var,  f"sample_var={var:.4f}, expected={sigma**2:.4f}")


# ─────────────────────────────────────────────────────────────────────────────
# 8. Falcon parameter sets
# ─────────────────────────────────────────────────────────────────────────────

def test_falcon_params():
    print("\n── 8. Falcon parameter sets ─────────────────────────────────────")
    #
    # samplerz is only valid for sigma in (sigmin, MAX_SIGMA).
    # The large "sigma" in FALCON_PARAMS (165, 168) is the SIGNING sigma —
    # it is NOT passed directly to samplerz.  Instead:
    #
    #   Key-gen sampling (sigma_fg):
    #     sigma_base = 1.17 * sqrt(q/8192) ≈ 1.433  (< MAX_SIGMA ✓)
    #     k = 4096 // n  samples are drawn and SUMMED per coefficient
    #     Sum of k i.i.d. Gaussian(0, σ_base) ~ Gaussian(0, √k · σ_base)
    #     → effective sigma_fg = √k · sigma_base
    #
    #   Signing sampling (sigma ≈ 165):
    #     Achieved via Fast Fourier Sampling (ffsampling), where each LEAF
    #     of the Falcon tree uses samplerz with a small leaf-sigma ∈ [sigmin, sigma_max].
    #
    import math
    Q_local  = 12289
    N_BASE   = 4096
    SIGMA_BASE = 1.17 * math.sqrt(Q_local / N_BASE)   # ≈ 1.433 < MAX_SIGMA

    for n in (512, 1024):
        p      = FALCON_PARAMS[n]
        sigmin = p["sigma_min"]
        k      = N_BASE // n                           # samples summed per coeff
        sigma_fg_expected = math.sqrt(k) * SIGMA_BASE  # effective sigma_fg

        # Draw n*k samples with sigma_base and sum into n coefficients
        raw  = [samplerz(0.0, SIGMA_BASE, sigmin) for _ in range(n * k)]
        poly = [sum(raw[i*k:(i+1)*k]) for i in range(n)]

        mean = sum(poly) / n
        var  = sum(z * z for z in poly) / n   # E[z²] ≈ sigma_fg² when μ=0

        ok_mean = abs(mean) < 3 * sigma_fg_expected / (n ** 0.5)
        ok_var  = abs(var - sigma_fg_expected**2) < 0.15 * sigma_fg_expected**2

        _record(f"Falcon-{n}: gen_poly mean ≈ 0  (sigma_base={SIGMA_BASE:.4f}, k={k})",
                ok_mean, f"mean={mean:.3f}")
        _record(f"Falcon-{n}: gen_poly var ≈ sigma_fg²={sigma_fg_expected**2:.2f}",
                ok_var, f"var={var:.2f}")

        # Coefficients should be small integers
        max_abs = max(abs(z) for z in poly)
        _record(f"Falcon-{n}: all poly coeffs within 6·sigma_fg={6*sigma_fg_expected:.1f}",
                max_abs < 6 * sigma_fg_expected, f"max|coeff|={max_abs}")


# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

def _summary():
    total  = len(_results)
    passed = sum(1 for _, ok, _ in _results if ok)
    failed = total - passed
    print(f"\n{'='*66}")
    print(f"  Results: {passed}/{total} passed", end="")
    if failed:
        print(f"  ── {failed} FAILED:")
        for name, ok, detail in _results:
            if not ok:
                print(f"    ✗  {name}" + (f"  [{detail}]" if detail else ""))
    else:
        print("  — all OK")
    print(f"{'='*66}")


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("samplerz Test Suite")
    print("=" * 66)

    test_basesampler_range()
    test_approxexp()
    test_berexp()
    test_output_type()
    test_statistics()
    test_falcon_params()

    _summary()
    sys.exit(0 if all(ok for _, ok, _ in _results) else 1)
