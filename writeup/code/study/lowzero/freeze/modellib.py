"""modellib.py: random models of h -> ln|L(1/2 + i(t + h), chi)| on one interval h in [0, theta), theta = 2 pi,
and the per-interval observables (same definitions as obslib.py).

Prime probe (exact, the Euler factor at p with a uniform phase; chi(p) = +-1 is absorbed in the phase):
    Y_p(h; phi) = -ln|1 - p^(-1/2) exp(i(phi + h ln p))|.
Models of the field (S = the set of primes dividing the conductor):
  'hyb<k0>'  CUE(N) on the full circle mapped onto the interval (FHK's matching: N zeros per 2 pi <-> L zeros
             per 2 pi in t), with its Fourier modes k < k0 removed and replaced by the explicit Euler factors of the
             primes p with ln p < k0 - 1/2, p not in S:
               W(h) = ln|Lambda_N(e^{ih})| + sum_{k<k0} Re(e^{ikh} T_k)/k + sum_{p: ln p < k0-1/2, p not in S} Y_p(h),
             T_k = sum_j e^{-ik theta_j}  (ln|Lambda| = -Re sum_k e^{ikh} T_k / k).
  'cue'      CUE(N) with the probes added on top (no mode removal; double-counts the low modes).
  'ep'       randomized Euler product: sum_{p <= X, p not in S} Y_p(h), ln X = e^{-gamma} L (GHK: E|P_X|^2 matches
             the second moment's leading term); no zeros.
Grid: K = 126 points h_k = k theta/K (step 0.0499); Z_beta = mean_k e^{2 beta W}; M = max with a parabola.
"""
import numpy as np

THETA = 2 * np.pi
K = 126
HGRID = np.arange(K) * THETA / K
BETAS = (0.5, 1.0, 1.5, 2.0, 3.0)
EULER_GAMMA = 0.5772156649015329


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def cue_angles(rng, B, N):
    Z = (rng.standard_normal((B, N, N)) + 1j * rng.standard_normal((B, N, N))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    d = np.diagonal(R, axis1=1, axis2=2)
    Q = Q * (d / np.abs(d))[:, None, :]
    return np.angle(np.linalg.eigvals(Q))


def cue_field(rng, B, N, k0=1, h=HGRID):
    """ln|Lambda_N(e^{ih})| with modes k < k0 removed; (B, len(h))."""
    th = cue_angles(rng, B, N)
    W = np.zeros((B, len(h)))
    for j in range(N):
        W += np.log(np.abs(2 * np.sin((h[None, :] - th[:, j:j + 1]) / 2)) + 1e-300)
    for k in range(1, k0):
        Tk = np.exp(-1j * k * th).sum(axis=1)
        W += np.real(np.exp(1j * k * h)[None, :] * Tk[:, None]) / k
    return W


def probe(rng, B, p, h=HGRID, phi=None):
    if phi is None:
        phi = rng.uniform(0, 2 * np.pi, B)
    return -np.log(np.abs(1 - p ** -0.5 * np.exp(1j * (phi[:, None] + h[None, :] * np.log(p)))))


def ep_field(rng, B, L, skip=(), h=HGRID, chunk=200):
    X = int(np.exp(np.exp(-EULER_GAMMA) * L))
    ps = [p for p in primes_upto(max(X, 2)) if p not in skip]
    W = np.zeros((B, len(h)))
    for p in ps:
        W += probe(rng, B, p, h)
    return W


def stats(W, betas=BETAS):
    """per-row observables: dict of arrays (B,)."""
    B, Kp = W.shape
    out = {}
    km = np.argmax(W, axis=1)
    r = np.arange(B)
    b = W[r, km]
    kl = np.clip(km - 1, 0, Kp - 1); kr = np.clip(km + 1, 0, Kp - 1)
    a = W[r, kl]; c = W[r, kr]
    den = a - 2 * b + c
    inner = (km > 0) & (km < Kp - 1) & (den < 0)
    d = np.where(inner, 0.5 * (a - c) / np.where(den < 0, den, -1), 0.0)
    d = np.clip(d, -1, 1)
    out['M'] = np.where(inner, b - 0.25 * (a - c) * d, b)
    out['Mgrid'] = b
    mx = b[:, None]
    for be in betas:
        out['lnZ_%g' % be] = 2 * be * b + np.log(np.mean(np.exp(2 * be * (W - mx)), axis=1))
    out['mean_ln'] = W.mean(axis=1)
    return out
