"""Referee check of the exponent A(kappa) = 2 + 8 kappa of rh_proof.md Prop. 7.4 (own code).

A factor of length N = K^beta (K = size of the k-family) has R large values at level sigma with
  R <= K^{e},  e = min( 1 (trivial, Prop 7.3),  min_g [ max(1, g beta) + g beta (1 - 2 sigma) ] (Prop 7.1) ),
and, for a long smooth factor, the exclusion (Prop 7.2): no large values if sigma >= 1/2 + 2 kappa / beta.
The method needs e <= A (1 - sigma); A(kappa) = sup e/(1-sigma).
Exact rational arithmetic is not needed: a fine grid plus a local refinement around the maximiser.
"""
import numpy as np


def ratio(beta, sigma, kappa, gmax=80, exclusion=True):
    e = np.ones_like(beta)
    for g in range(1, gmax + 1):
        e = np.minimum(e, np.maximum(1, g * beta) + g * beta * (1 - 2 * sigma))
    r = e / (1 - sigma)
    if exclusion:
        r = np.where(sigma >= 0.5 + 2 * kappa / beta, 0.0, r)
    return r


def scan(kappa, exclusion=True, bmin=0.02, bmax=8.0, gmax=80):
    B, S = np.meshgrid(np.linspace(bmin, bmax, 2400), np.linspace(0.5, 0.999, 900), indexing='ij')
    R = ratio(B, S, kappa, gmax=gmax, exclusion=exclusion)
    i = np.unravel_index(np.argmax(R), R.shape)
    b0, s0 = B[i], S[i]
    # refine
    B2, S2 = np.meshgrid(np.linspace(max(bmin, b0 - 0.01), b0 + 0.01, 801),
                         np.linspace(max(0.5, s0 - 0.002), min(0.9995, s0 + 0.002), 801), indexing='ij')
    R2 = ratio(B2, S2, kappa, gmax=gmax, exclusion=exclusion)
    j = np.unravel_index(np.argmax(R2), R2.shape)
    return R2[j], B2[j], S2[j]


if __name__ == '__main__':
    print("kappa   sup e/(1-sigma)   2+8kappa   at beta   sigma   predicted beta=1+4k, sigma=1/2+2k/(1+4k)")
    for kappa in (1 / 8, 1 / 6, 0.2, 0.2012, 0.24, 0.249, 0.25):
        r, b, s = scan(kappa)
        print(f"{kappa:.4f}  {r:14.4f}   {2 + 8 * kappa:8.4f}   {b:7.4f}  {s:7.4f}   {1 + 4 * kappa:.4f} {0.5 + 2 * kappa / (1 + 4 * kappa):.4f}")
    r, b, s = scan(0.2, exclusion=False)
    print(f"long factors without the exclusion: sup = {r:.4f} at beta {b:.4f}, sigma {s:.4f} (= 1/(1-sigma) at the grid edge: unbounded)")
    r, b, s = scan(0.2, exclusion=False, bmin=0.005, bmax=0.05, gmax=400)
    print(f"short factors beta in [0.005, 0.05], no exclusion: sup = {r:.4f} (2(1+beta) bound)")
    # what the Weyl exponent and MV-type exponents give for theta = 1/A
    for name, k in (("Weyl kappa=1/6", 1 / 6), ("convexity kappa=1/4", 0.25),
                    ("Burgess-type GL1 kappa = 1/4 - (1-2*7/64)/16", 0.25 - (1 - 14 / 64) / 16)):
        print(f"{name:45s} theta < {1 / (2 + 8 * k):.5f}")
