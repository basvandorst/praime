# Exponent A (R << (MT)^{A(1-sigma)}) in Stucky's scheme when the false "Ricci C" estimate D = MT + N min(M,T)
# is replaced by the elementary large sieve D = N + MT + max(M,T) sqrt(N) (valid for arbitrary coefficients,
# proof: Beurling-Selberg majorant + Schur test + lattice points in a polar box).
import numpy as np
betas = np.concatenate([np.linspace(0.05, 3, 1180), np.linspace(3.005, 12, 300)])[:, None]
sigmas = np.linspace(0.5004, 0.9996, 1249)[None, :]

def worstA(L, kappa, use_C=False, use_elem=True, use_hux=True, gmax=130):
    b, s = betas, sigmas
    best = np.ones((b.shape[0], s.shape[1]))
    for g in range(1, gmax + 1):
        best = np.minimum(best, np.maximum(1 + L + g*b*(1 - 2*s), 2*g*b*(1 - s)))          # Coleman LS
        if use_C:
            best = np.minimum(best, np.maximum(1 + g*b*(1 - 2*s), (1 - L)/2 + 2*g*b*(1 - s)))
        if use_elem:
            e = np.maximum(np.maximum(1 + g*b*(1 - 2*s), 2*g*b*(1 - s)), (1 + L)/2 + g*b*(1.5 - 2*s))
            best = np.minimum(best, e)
        if use_hux:
            best = np.minimum(best, np.where(s >= 0.75, np.maximum(1 + L + g*b*(4 - 6*s), 2*g*b*(1 - s)), np.inf))
    A = best/(1 - s)
    if kappa is not None:
        A = np.where((b > 2/3) & (b*(s - 0.5) > kappa*(1 + L) + 1e-12), 0.0, A)
    i, j = np.unravel_index(np.argmax(A), A.shape)
    return A[i, j], float(betas[i, 0]), float(sigmas[0, j])

for L in [1.0, 0.9, 0.6, 0.3, 0.0]:
    for kappa, kn in [(1/6, "Weyl 1/6"), (0.2, "0.20"), (0.24, "0.24"), (1/4, "convexity 1/4")]:
        A, bw, sw = worstA(L, kappa)
        print(f"elementary LS, no RicciC: L={L:.2f} kappa={kn:14s} A={A:.4f} 1/A={1/A:.4f} worst beta={bw:.3f} sigma={sw:.4f}", flush=True)
