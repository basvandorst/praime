"""cue_check.py: the Palm (two- and three-point) formula for Var log|Lambda'(theta_n)| in CUE(N), against the exact value.

Exact (Hughes-Keating-O'Connell, finite N): K''(0) = sum_{j=1}^{N-1} [psi'(j+2) - psi'(j+1)/2].
Palm formula: with a zero at 0, the other eigenangles have the one- and two-point densities R2(0,u)/rho1, R3(0,u,v)/rho1,
rho1 = N/2pi, K(x) = (1/2pi) sum_{m=0}^{N-1} e^{imx}; f(u) = ln|2 sin(u/2)|, Y = sum_{j != n} f(theta_j - theta_n):
  Var Y = int f^2 R2/rho1 + int int f(u) f(v) [R3(0,u,v)/rho1 - R2(0,u) R2(0,v)/rho1^2].
(1) exact evaluation by Fourier series (f_k = -1/(2|k|)): closed finite sums (derivation.md sect. 5);
(2) direct 2-D quadrature of the determinants (Gauss panels graded at the log singularities u = 0, v = 0), the same
    scheme as the L-function integrals; (3) the pieces A, B, V3, mu^2 of derivation.md.
usage: cue_check.py
"""
import math, numpy as np
from mpmath import psi, mpf, harmonic


def exact_var(N): return float(sum(psi(1, j + 2) - psi(1, j + 1) / 2 for j in range(1, N)))


def palm_fourier(N):
    H = [float(harmonic(k)) for k in range(N + 1)]
    B = N * math.pi ** 2 / 12 - sum((N - k) / (2 * k * k) for k in range(1, N))
    A = -math.pi ** 2 / 12 - sum((N - k) * (H[k - 1] + H[k]) / k for k in range(1, N)) / N
    V3 = sum((H[b] + H[N - 1 - b]) ** 2 for b in range(N)) / (2 * N)
    mu = H[N] - 1
    return dict(A=A, B=B, V3=V3, mu=mu, var=A + B + V3 - mu * mu)


def nodes_axis(N, npan_per=4, ng=10):
    """nodes on [0, pi]: graded panels at 0 (log singularity), then panels of width ~ pi/(npan_per*N)."""
    x, w = np.polynomial.legendre.leggauss(ng)
    edges = [0.0]
    h = math.pi / (npan_per * N)
    for k in range(6, 0, -1): edges.append(h * 4.0 ** (-k))
    e = h
    while e < math.pi - 1e-12: edges.append(e); e += h
    edges.append(math.pi); edges = np.unique(np.array(edges))
    a, b = edges[:-1], edges[1:]
    X = ((b - a)[:, None] / 2 * (x[None, :] + 1) + a[:, None]).ravel(); W = ((b - a)[:, None] / 2 * w[None, :]).ravel()
    return X, W


def palm_quad(N):
    X, W = nodes_axis(N); X = np.concatenate([-X[::-1], X]); W = np.concatenate([W[::-1], W])
    rho = N / (2 * math.pi)
    def K(x):  # (1/2pi) sum e^{imx}
        m = np.arange(N); return np.exp(1j * np.multiply.outer(x, m)).sum(-1) / (2 * math.pi)
    f = np.log(np.abs(2 * np.sin(X / 2)))
    Ku = K(X); R2 = rho * rho - np.abs(Ku) ** 2
    I1 = np.sum(W * f * f * R2) / rho
    mu = np.sum(W * f * R2) / rho
    I2 = 0.0
    for i in range(len(X)):
        u = X[i]; v = X
        Kuv = K(u - v)
        R3 = rho ** 3 - rho * (np.abs(Ku[i]) ** 2 + np.abs(Ku) ** 2 + np.abs(Kuv) ** 2) + 2 * np.real(np.conj(Ku[i]) * Kuv * Ku)
        I2 += W[i] * f[i] * np.sum(W * f * (R3 / rho - R2[i] * R2 / rho ** 2))
    return dict(I1=I1, I2=I2, mu=mu, var=I1 + I2)


if __name__ == '__main__':
    print('N    exact K\'\'(0)      Palm (Fourier sums)  diff      Palm (2-D quadrature) diff     mean H_N-1 / quad')
    for N in (2, 10, 20, 50):
        ex = exact_var(N); pf = palm_fourier(N)
        pq = palm_quad(N) if N >= 10 else None
        print('%-4d %.12f  %.12f  %.1e' % (N, ex, pf['var'], pf['var'] - ex) +
              ('   %.10f  %.1e   %.8f / %.8f' % (pq['var'], pq['var'] - ex, pf['mu'], pq['mu']) if pq else ''))
        print('     pieces (Fourier): A %.6f  B %.6f  V3 %.6f  mu^2 %.6f   (B - [ln N/2 + (1+gamma)/2] = %.2e)' % (
            pf['A'], pf['B'], pf['V3'], pf['mu'] ** 2, pf['B'] - (0.5 * math.log(N) + (1 + 0.5772156649015329) / 2)))
    for N in (10, 20, 50, 200, 1000):
        print('N=%5d  Var - (ln N/2 - 0.678785) = %+.6f' % (N, exact_var(N) - (0.5 * math.log(N) + (0.5772156649015329 + 1) / 2 + 1 - math.pi ** 2 / 4)))
