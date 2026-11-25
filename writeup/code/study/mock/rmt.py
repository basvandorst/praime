"""Random-matrix (U(N), N -> infinity) cumulants of the linear statistic Z_f = sum_j f(N theta_j / 2 pi).

kappa_m^U(f) = int_{u_1+..+u_m=0} prod_j fhat(u_j) Q_m(u) du_1..du_{m-1}
Q_m(u) = sum over ordered set partitions (B_1..B_n) of {1..m} of (-1)^(n+1)/n * (1 - R)_+,
         R = range of the partial sums S_0=0, S_j = w_1+..+w_j (j < n), w_i = sum_{l in B_i} u_l.
(derivation.md section 2: log det(I + (e^{tf}-1) K_sine), K_sine has Fourier transform 1_[-1/2,1/2].)
Cross-check: toeplitz_cumulants(N) = exact finite-N cumulants from log det T_N(e^{tF}) (Heine).
Self-contained: numpy only.
"""
import itertools
import numpy as np


def fhat_tri(u, s):
    """Fejer test function: fhat(u) = (1 - |u|/s)_+, f(x) = s (sin(pi s x)/(pi s x))^2."""
    return np.clip(1.0 - np.abs(u) / s, 0.0, None)


def ordered_set_partitions(m):
    """All ordered set partitions of range(m) as lists of tuples."""
    out = []
    for labels in itertools.product(range(m), repeat=m):
        n = max(labels) + 1
        if set(labels) != set(range(n)):
            continue
        blocks = [tuple(i for i in range(m) if labels[i] == b) for b in range(n)]
        out.append(blocks)
    return out


_OSP = {m: ordered_set_partitions(m) for m in (2, 3, 4, 5)}


def Q(us):
    """us: array (m, K) of points on the hyperplane sum = 0. Returns Q_m at each point."""
    m = us.shape[0]
    tot = np.zeros(us.shape[1])
    for blocks in _OSP[m]:
        n = len(blocks)
        S = np.zeros(us.shape[1])
        mx = np.zeros(us.shape[1])
        mn = np.zeros(us.shape[1])
        for b in blocks[:-1]:
            S = S + us[list(b)].sum(axis=0)
            mx = np.maximum(mx, S)
            mn = np.minimum(mn, S)
        tot += ((-1) ** (n + 1) / n) * np.clip(1.0 - (mx - mn), 0.0, None)
    return tot


def kappa_U(m, fhat, s, h=None, chunk=2_000_000):
    """N=infinity cumulant by the midpoint rule on [-s, s]^(m-1); u_m = -sum."""
    if m == 2:
        u = np.linspace(-s, s, 200001)
        y = np.minimum(np.abs(u), 1.0) * fhat(u, s) ** 2
        return np.trapezoid(y, u) if hasattr(np, 'trapezoid') else np.trapz(y, u)
    if h is None:
        h = s / (600 if m == 3 else 90)
    g = -s + (np.arange(int(round(2 * s / h))) + 0.5) * h
    total = 0.0
    if m == 3:
        U1, U2 = np.meshgrid(g, g, indexing='ij')
        U1 = U1.ravel(); U2 = U2.ravel(); U3 = -U1 - U2
        keep = np.abs(U3) < s
        U1, U2, U3 = U1[keep], U2[keep], U3[keep]
        us = np.vstack([U1, U2, U3])
        w = fhat(U1, s) * fhat(U2, s) * fhat(U3, s)
        total = np.sum(w * Q(us)) * h * h
        return total
    if m == 4:
        G1, G2 = np.meshgrid(g, g, indexing='ij')
        G1 = G1.ravel(); G2 = G2.ravel()
        for u3 in g:
            U4 = -G1 - G2 - u3
            keep = np.abs(U4) < s
            if not keep.any():
                continue
            a, b, d = G1[keep], G2[keep], U4[keep]
            c = np.full(a.shape, u3)
            w = fhat(a, s) * fhat(b, s) * fhat(c, s) * fhat(d, s)
            nz = w > 0
            if not nz.any():
                continue
            us = np.vstack([a[nz], b[nz], c[nz], d[nz]])
            total += np.sum(w[nz] * Q(us))
        return total * h ** 3
    raise ValueError(m)


def toeplitz_cumulants(N, fhat, s):
    """Exact cumulants 1..4 of Z = sum_n (1/N) fhat(n/N) Tr U^n over U(N) (Heine: E e^{tZ} = det T_N(e^{tF}))."""
    K = int(np.floor(s * N))
    n = np.arange(-K, K + 1)
    c = fhat(n / N, s) / N                       # Fourier coefficients of F
    # powers F^k as coefficient arrays (convolution)
    pw = [None, c]
    for k in range(2, 5):
        pw.append(np.convolve(pw[-1], c))
    def T(coef):
        L = (len(coef) - 1) // 2
        idx = np.arange(N)[:, None] - np.arange(N)[None, :]
        M = np.zeros((N, N))
        ok = np.abs(idx) <= L
        M[ok] = coef[idx[ok] + L]
        return M
    A = [None] + [T(pw[k]) for k in range(1, 5)]
    tr = np.trace
    A1, A2, A3, A4 = A[1], A[2], A[3], A[4]
    A11 = A1 @ A1
    k1 = tr(A1)
    k2 = tr(A2) - tr(A11)
    k3 = tr(A3) - 3 * tr(A1 @ A2) + 2 * tr(A11 @ A1)
    k4 = tr(A4) - 4 * tr(A1 @ A3) - 3 * tr(A2 @ A2) + 12 * tr(A11 @ A2) - 6 * tr(A11 @ A11)
    return k1, k2, k3, k4


if __name__ == '__main__':
    import sys
    for s in (0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0):
        k2 = kappa_U(2, fhat_tri, s)
        k3 = kappa_U(3, fhat_tri, s)
        k4 = kappa_U(4, fhat_tri, s, h=s / 60)
        print(f"s={s:4.2f} k2U={k2:.6f} k3U={k3:+.6f} k4U={k4:+.6f}  k3/k2^1.5={k3/k2**1.5:+.4f} k4/k2^2={k4/k2**2:+.4f}")
        sys.stdout.flush()
