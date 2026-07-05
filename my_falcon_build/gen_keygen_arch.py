"""Generate Falcon Keygen Architecture diagram as a Word document."""

from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH

doc = Document()

# ── Title ──
title = doc.add_heading("Falcon Key Generation Architecture", level=0)
title.alignment = WD_ALIGN_PARAGRAPH.CENTER

doc.add_paragraph("")

# ── Section 1: Pipeline Overview ──
doc.add_heading("1. Keygen Pipeline Overview", level=1)
code = doc.add_paragraph()
code.style = doc.styles['Normal']
run = code.add_run("""
                      keygen(n=512)
                           │
             ┌─────────────┼─────────────┐
             │             │             │
             ▼             ▼             ▼
        ntru_gen()    poly_div_ntt()  _build_tree()
        [ntrugen.py]  [keygen.py]    [keygen.py]
             │             │             │
             ▼             ▼             ▼
        (f, g, F, G)      h         (B_fft, T)
             │             │             │
             └─────────────┼─────────────┘
                           ▼
                  sk = (f, g, F, G, B_fft, T)
                  pk = (h)
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)


# ── Section 2: Module Dependency Map ──
doc.add_heading("2. Module Dependency Map", level=1)
code = doc.add_paragraph()
run = code.add_run("""
      ┌──────────┐
      │ params.py│  Q=12289, sigma, sigma_min
      └────┬─────┘
           │
      ┌────▼──���──┐     ┌──────────┐
      │  ntt.py  │     │samplerz.py│  Discrete Gaussian sampler
      └────┬─────┘     └─────┬─────┘
           │                 │
      ┌���───▼─────┐     ┌─────▼─────┐
      │  fft.py  │────▶│ntrugen.py │  NTRU key generation
      └────┬─────┘     └─────┬─────┘
           │                 │
      ┌────▼─────────────────▼──────┐
      │      ffsampling.py          │  Gram matrix + LDL tree
      └─────────────┬───────────────┘
                    │
      ┌─────────────▼───────────────┐
      │        keygen.py            │  Wires everything together
      └─────────���───────────────────┘
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)


# ── Section 3: Step-by-Step ──
doc.add_heading("3. Step-by-Step Function Calls", level=1)

doc.add_heading("Step 1: Generate NTRU Key Quadruple", level=2)
code = doc.add_paragraph()
run = code.add_run("""ntru_gen(512)
    ├── gen_poly(512) ──────────────▶ f  (small Gaussian polynomial)
    ├── gen_poly(512) ──────────────▶ g  (small Gaussian polynomial)
    ├─��� gs_norm(f, g, q) ──────────▶ Check norm ≤ 1.17² · q
    ├── ntt(f) ────────────────────▶ Check f invertible mod q
    └── ntru_solve(f, g) ──────────▶ (F, G) solving f·G - g·F = q
            ├── field_norm(f), field_norm(g)   Reduce degree n → n/2
            ├── ntru_solve(fp, gp)             Recurse until n=1
            │       └── xgcd(f0, g0)           Base case: extended GCD
            ├── lift(Fp), lift(Gp)             Lift back to degree n
            ├── karamul(lifted, conjugate)      Multiply over Z (exact)
            └── reduce(f, g, F, G)             Babai reduction

OUTPUT: (f, g, F, G) — four integer polynomials in Z[x]/(xⁿ+1)
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)

doc.add_heading("Step 2: Compute Public Key", level=2)
code = doc.add_paragraph()
run = code.add_run("""poly_div_ntt(g, f)
    ├── ntt(g mod q) ──────────────▶ g_ntt
    ├── ntt(f mod q) ──────────────▶ f_ntt
    ├── h_ntt[i] = g[i] · f[i]⁻¹ mod q   (Fermat inverse: pow(f, q-2, q))
    └── intt(h_ntt) ───────────────▶ h

OUTPUT: h — public key polynomial in Z_q[x]/(xⁿ+1)
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)

doc.add_heading("Step 3: Build LDL Tree (tree mode only)", level=2)
code = doc.add_paragraph()
run = code.add_run("""_build_tree(512, f, g, F, G)
    │
    ├── B = [[g, -f], [G, -F]]            Basis matrix (coeff domain)
    │
    ├── gram(B)                            Gram matrix G = B* · B
    │       └── G[i][j] = Σ_k B[i][k] · adj(B[j][k])
    │               Uses: poly_mul_coeff, poly_adj_coeff, poly_add_coeff
    │
    ├── fft(G_mat) ────────────────▶ G_fft  (2×2 matrix in FFT domain)
    │
    ├── ffldl_fft(G_fft)                   Build LDL decomposition tree
    │       ├── ldl_fft(G)                 2×2 LDL: G = L · D · L*
    │       │       D00 = G[0][0]
    │       │       L10 = G[1][0] / G[0][0]
    │       │       D11 = G[1][1] - L10·adj(L10)·G[0][0]
    │       ├── poly_split_fft(D00)        Bisect diagonal entries
    │       ├── poly_split_fft(D11)        Bisect diagonal entries
    │       ├── ffldl_fft(G_left)          Recurse on left half
    │       ├── ffldl_fft(G_right)         Recurse on right half
    │       └── ... until leaf (n=2)
    │
    ├── normalize_tree(T, sigma)           Convert leaves for sampling
    │       Internal node (len==3): recurse on children
    │       Leaf (len==2): value = sigma / sqrt(value)
    │
    └── B_fft = fft(each poly in B)        FFT-domain basis for signing
            [[fft(g), fft(-f)],
             [fft(G), fft(-F)]]

OUTPUT: (B_fft, T) — precomputed signing data
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)

doc.add_heading("Step 4: Package Keys", level=2)
code = doc.add_paragraph()
run = code.add_run("""sk = SecretKey(n, mode, f, g, F, G, B_fft, T)
pk = PublicKey(n, h)
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)


# ── Section 4: Two Signing Modes ──
doc.add_heading("4. Two Signing Modes", level=1)
code = doc.add_paragraph()
run = code.add_run("""
  TREE MODE (mode="tree")              DYNAMIC MODE (mode="dynamic")
  ┌───────────────────────┐            ┌───────────────────────┐
  │ keygen(512, "tree")   │            │keygen(512, "dynamic") │
  │                       │            │                       │
  │ sk contains:          │            │ sk contains:          │
  │   f, g, F, G          │            │   f, g, F, G          │
  │   B_fft  ✓            │            │   B_fft = None        │
  │   T      ✓            │            │   T     = None        │
  │                       │            │                       │
  │ Signing: FAST         │            │ Before signing call:  │
  │ (tree ready to use)   │            │ expand_secret_key(sk) │
  │                       │            │ rebuilds B_fft + T    │
  │ Key size: LARGER      │            │                       │
  └───────────────────────┘            │ Signing: SLOWER       │
                                       │ Key size: SMALLER     │
                                       └───────────────────────┘
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)


# ── Section 5: Data Flow ��─
doc.add_heading("5. Data Flow Diagram", level=1)
code = doc.add_paragraph()
run = code.add_run("""
      Gaussian         Gaussian
      Sampler          Sampler
         │                │
         ▼                ▼
         f                g          ← Small secret polynomials (~6 bits)
         │                │
         ├────────────────┤
         │   ntru_solve   │
         │  f·G - g·F = q │
         ▼                ▼
         F                G          ← Larger polynomials (~8 bits)
         │                │
         ├────┬───────────┤
         │    │           │
         │    ▼           │
         │  B=[[g,-f],    │
         │    [G,-F]]     │
         │    │           │
         │    ▼           │
         │  gram(B)       │
         │    │           │
         │    ▼           │
         │  ffldl_fft     │
         │    │           │
         │    ▼           │
         │  normalize     │
         │    │           │
         │    ▼           │
         │   T (tree) ────┤──────▶  SECRET KEY
         │                │         (f, g, F, G, B_fft, T)
         │   B_fft  ──────┤
         │                │
         ▼                ▼
      h = g/f mod q  ────────────▶  PUBLIC KEY
      (NTT division)                (h)
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)


# ── Section 6: Tree Structure ──
doc.add_heading("6. LDL Tree Structure (after ffldl_fft + normalize)", level=1)
code = doc.add_paragraph()
run = code.add_run("""
              [L10, T_left, T_right]              depth 0  (n=512)
             /                      \\
    [L10, T_l, T_r]          [L10, T_l, T_r]     depth 1  (n=256)
   /              \\         /              \\
  ...             ...     ...             ...
 /                  \\
[sigma, 0]      [sigma, 0]                        depth 9  (leaf)

Internal node: [L10_fft_array, T_left, T_right]   len == 3
Leaf:          [sigma_value, 0]                    len == 2

Total depth:  log2(512) = 9
Total leaves: 512
Each leaf sigma = sigma_signing / sqrt(||b_i||²)
""")
run.font.name = 'Courier New'
run.font.size = Pt(9)


# ── Save ──
path = "/Users/mac/Desktop/project/phd_code/my_falcon_build/Falcon_Keygen_Architecture.docx"
doc.save(path)
print(f"Saved to {path}")
