/*
 * LDL Tree construction for Falcon.
 *
 * This file contains:
 *   - ffLDL_fft()             recursive LDL decomposition of Gram matrix
 *   - ffLDL_binary_normalize() normalize tree leaves with sigma
 *   - expand_privkey()        top-level: builds B0 matrix + LDL tree
 *                             from raw private key (f, g, F, G)
 *
 * Output is the expanded_key[] array consumed by gaussian_sampling.c
 *
 * ==========================(LICENSE BEGIN)============================
 *
 * Copyright (c) 2017-2019  Falcon Project
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   Thomas Pornin <thomas.pornin@nccgroup.com>
 */

#include "inner.h"

#define MKN(logn)   ((size_t)1 << (logn))

/* ====================================================================== */
/*
 * LDL tree size.
 *
 * For logn = 0 (polynomials are constant), the "tree" is a single element.
 * Otherwise, the tree node has size 2^logn, and has two child trees for
 * size logn-1 each:
 *   s(0)     = 1
 *   s(logn)  = 2^logn + 2*s(logn-1)
 */
static inline unsigned
ffLDL_treesize(unsigned logn)
{
	return (logn + 1) << logn;
}

/* ====================================================================== */
/*
 * Expanded private key layout offsets.
 * The expanded key stores B0 (four polynomials) followed by the LDL tree.
 */
static inline size_t skoff_b00(unsigned logn) { (void)logn; return 0; }
static inline size_t skoff_b01(unsigned logn) { return MKN(logn); }
static inline size_t skoff_b10(unsigned logn) { return 2 * MKN(logn); }
static inline size_t skoff_b11(unsigned logn) { return 3 * MKN(logn); }
static inline size_t skoff_tree(unsigned logn) { return 4 * MKN(logn); }

/* ====================================================================== */
/*
 * Convert a small-integer polynomial (int8_t) to fpr (IEEE-754 double).
 */
static void
smallints_to_fpr(fpr *r, const int8_t *t, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		r[u] = fpr_of(t[u]);
	}
}

/* ====================================================================== */
/*
 * Inner recursive function for ffLDL_fft().
 * Expects the matrix to be auto-adjoint and quasicyclic.
 * Uses g0 and g1 as modifiable temporaries.
 * tmp[] must have room for at least one polynomial.
 */
static void
ffLDL_fft_inner(fpr *restrict tree,
	fpr *restrict g0, fpr *restrict g1, unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;

	n = MKN(logn);
	if (n == 1) {
		tree[0] = g0[0];
		return;
	}
	hn = n >> 1;

	/*
	 * LDL decomposition: yields L (written into tree) and d11 (into tmp).
	 * d00 stays in g0.
	 */
	Zf(poly_LDLmv_fft)(tmp, tree, g0, g1, g0, logn);

	/*
	 * Split d00 (in g0) and d11 (in tmp) into half-size polynomials.
	 * Reuse g0 and g1 as temporary storage:
	 *   d00 splits into g1, g1+hn
	 *   d11 splits into g0, g0+hn
	 */
	Zf(poly_split_fft)(g1, g1 + hn, g0, logn);
	Zf(poly_split_fft)(g0, g0 + hn, tmp, logn);

	/*
	 * Recurse on both halves.
	 */
	ffLDL_fft_inner(tree + n,
		g1, g1 + hn, logn - 1, tmp);
	ffLDL_fft_inner(tree + n + ffLDL_treesize(logn - 1),
		g0, g0 + hn, logn - 1, tmp);
}

/*
 * Compute the ffLDL tree of an auto-adjoint matrix G.
 * G is provided as three polynomials in FFT representation (upper triangle).
 *
 * Output tree[] has size (logn+1)*(2^logn) elements.
 * tmp[] must have room for at least three polynomials.
 */
static void
ffLDL_fft(fpr *restrict tree, const fpr *restrict g00,
	const fpr *restrict g01, const fpr *restrict g11,
	unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	fpr *d00, *d11;

	n = MKN(logn);
	if (n == 1) {
		tree[0] = g00[0];
		return;
	}
	hn = n >> 1;
	d00 = tmp;
	d11 = tmp + n;
	tmp += n << 1;

	memcpy(d00, g00, n * sizeof *g00);
	Zf(poly_LDLmv_fft)(d11, tree, g00, g01, g11, logn);

	Zf(poly_split_fft)(tmp, tmp + hn, d00, logn);
	Zf(poly_split_fft)(d00, d00 + hn, d11, logn);
	memcpy(d11, tmp, n * sizeof *tmp);
	ffLDL_fft_inner(tree + n,
		d11, d11 + hn, logn - 1, tmp);
	ffLDL_fft_inner(tree + n + ffLDL_treesize(logn - 1),
		d00, d00 + hn, logn - 1, tmp);
}

/*
 * Normalize the ffLDL tree: each leaf value x is replaced with
 * sqrt(x) * fpr_inv_sigma[orig_logn].
 *
 * This saves a division in the sampler (precomputes 1/sigma into leaves).
 */
static void
ffLDL_binary_normalize(fpr *tree, unsigned orig_logn, unsigned logn)
{
	size_t n;

	n = MKN(logn);
	if (n == 1) {
		tree[0] = fpr_mul(fpr_sqrt(tree[0]), fpr_inv_sigma[orig_logn]);
	} else {
		ffLDL_binary_normalize(tree + n, orig_logn, logn - 1);
		ffLDL_binary_normalize(tree + n + ffLDL_treesize(logn - 1),
			orig_logn, logn - 1);
	}
}

/* ====================================================================== */
/* see inner.h */

/*
 * Expand private key (f, g, F, G) into the expanded_key[] array.
 *
 * expanded_key[] layout:
 *   [0         .. n-1  ]  b00 = g   (FFT)
 *   [n         .. 2n-1 ]  b01 = -f  (FFT, negated)
 *   [2n        .. 3n-1 ]  b10 = G   (FFT)
 *   [3n        .. 4n-1 ]  b11 = -F  (FFT, negated)
 *   [4n        .. end  ]  LDL tree  (size = ffLDL_treesize(logn))
 */
void
Zf(expand_privkey)(fpr *restrict expanded_key,
	const int8_t *f, const int8_t *g,
	const int8_t *F, const int8_t *G,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n;
	fpr *rf, *rg, *rF, *rG;
	fpr *b00, *b01, *b10, *b11;
	fpr *g00, *g01, *g11, *gxx;
	fpr *tree;

	n = MKN(logn);
	b00 = expanded_key + skoff_b00(logn);
	b01 = expanded_key + skoff_b01(logn);
	b10 = expanded_key + skoff_b10(logn);
	b11 = expanded_key + skoff_b11(logn);
	tree = expanded_key + skoff_tree(logn);

	/*
	 * Load private key into B0 = [[g, -f], [G, -F]].
	 */
	rf = b01;
	rg = b00;
	rF = b11;
	rG = b10;

	smallints_to_fpr(rf, f, logn);
	smallints_to_fpr(rg, g, logn);
	smallints_to_fpr(rF, F, logn);
	smallints_to_fpr(rG, G, logn);

	/*
	 * FFT-transform all four key elements, then negate f and F.
	 */
	Zf(FFT)(rf, logn);
	Zf(FFT)(rg, logn);
	Zf(FFT)(rF, logn);
	Zf(FFT)(rG, logn);
	Zf(poly_neg)(rf, logn);
	Zf(poly_neg)(rF, logn);

	/*
	 * Compute Gram matrix G = B * B^†  (upper triangle only):
	 *   g00 = b00*adj(b00) + b01*adj(b01)
	 *   g01 = b00*adj(b10) + b01*adj(b11)
	 *   g11 = b10*adj(b10) + b11*adj(b11)
	 */
	g00 = (fpr *)tmp;
	g01 = g00 + n;
	g11 = g01 + n;
	gxx = g11 + n;

	memcpy(g00, b00, n * sizeof *b00);
	Zf(poly_mulselfadj_fft)(g00, logn);
	memcpy(gxx, b01, n * sizeof *b01);
	Zf(poly_mulselfadj_fft)(gxx, logn);
	Zf(poly_add)(g00, gxx, logn);

	memcpy(g01, b00, n * sizeof *b00);
	Zf(poly_muladj_fft)(g01, b10, logn);
	memcpy(gxx, b01, n * sizeof *b01);
	Zf(poly_muladj_fft)(gxx, b11, logn);
	Zf(poly_add)(g01, gxx, logn);

	memcpy(g11, b10, n * sizeof *b10);
	Zf(poly_mulselfadj_fft)(g11, logn);
	memcpy(gxx, b11, n * sizeof *b11);
	Zf(poly_mulselfadj_fft)(gxx, logn);
	Zf(poly_add)(g11, gxx, logn);

	/*
	 * Build the LDL tree from the Gram matrix.
	 */
	ffLDL_fft(tree, g00, g01, g11, logn, gxx);

	/*
	 * Normalize tree leaves: replace x with sqrt(x) * inv_sigma.
	 */
	ffLDL_binary_normalize(tree, logn, logn);
}
