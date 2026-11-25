# Re-derive the exponent A in R << (MT)^{A(1-sigma)} from the estimates listed in Stucky's paper
# (arXiv 2008.11325, sections 4-7), as a function of the pointwise bound exponent kappa for the long
# smooth polynomials (|L(1/2+it, xi)| << (m^2+t^2)^kappa: Weyl 1/6, convexity 1/4, none),
# with and without the "trivial" estimate D = MT + N min(M,T) (Ricci Thm C) and Huxley (Coleman 7.3).
# L = |log(M/T)|/log(MT); in the full-norm-range case theta = 1 one has T <= x^eps, so L -> 1.
# A sector width x^{-th} at full norm range is reached iff A < 1/th (MT <= x^{th+eps}).
import numpy as np

betas = np.concatenate([np.linspace(0.05, 3, 1180), np.linspace(3.005, 12, 300)])[:, None]
sigmas = np.linspace(0.5004, 0.9996, 1249)[None, :]

def worstA(L, kappa, use_C=True, use_hux=True, gmax=130):
    b, s = betas, sigmas
    best = np.ones((b.shape[0], s.shape[1]))          # R <= MT
    for g in range(1, gmax + 1):
        e1 = np.maximum(1 + L + g*b*(1 - 2*s), 2*g*b*(1 - s))
        best = np.minimum(best, e1)
        if use_C:
            e2 = np.maximum(1 + g*b*(1 - 2*s), (1 - L)/2 + 2*g*b*(1 - s))
            best = np.minimum(best, e2)
        if use_hux:
            e3 = np.where(s >= 0.75, np.maximum(1 + L + g*b*(4 - 6*s), 2*g*b*(1 - s)), np.inf)
            best = np.minimum(best, e3)
    A = best/(1 - s)
    if kappa is not None:
        excluded = (b > 2/3) & (b*(s - 0.5) > kappa*(1 + L) + 1e-12)
        A = np.where(excluded, 0.0, A)
    i, j = np.unravel_index(np.argmax(A), A.shape)
    return A[i, j], float(betas[i, 0]), float(sigmas[0, j])


import sys
cases = [(L, k, kn, C, H) for L in [1.0, 0.6, 0.0] for (k, kn) in [(1/6, "Weyl 1/6"), (1/4, "convexity 1/4"), (None, "none")]
         for C in [True, False] for H in [True, False]]
for L, kappa, kn, use_C, use_H in cases:
    A, bw, sw = worstA(L, kappa, use_C, use_H)
    print(f"L={L:.2f}  pointwise={kn:14s} RicciC={str(use_C):5s} Huxley={str(use_H):5s}  A={A:.4f}  1/A={1/A:.4f}  worst beta={bw:.3f} sigma={sw:.4f}", flush=True)
