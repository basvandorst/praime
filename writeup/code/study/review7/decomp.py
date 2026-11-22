"""Decompose the g = 2 long-range term at train into pair sums: A = sum w (rho_YY - 1), A+- = shifted-pair sums,
and compare with SXX, SXY; tests whether the first-order cancellation (total = 2A P^2) holds."""
import numpy as np, sys
sys.path.insert(0, '/home/bas/scratch/study-review7')
from rho_mine import primes_upto, logrho, Q, WIN
L0, W = WIN['train']; n = W // Q; g = 2
PR = primes_upto(Q * n + 64).astype(np.int64)
r = (Q / 48) / np.log(L0 + Q * np.arange(n) + Q / 2)
N2 = 1 << int(np.ceil(np.log2(2 * n))); F = np.fft.rfft(r, N2); w = np.fft.irfft(F * np.conj(F), N2)[:n]
m = r.sum()
rYY = np.exp(logrho(n, 0, 'YY', PR)); rXX = np.exp(logrho(n, g, 'XX', PR))
rXp = np.exp(logrho(n, g, 'XYp', PR)); rXm = np.exp(logrho(n, g, 'XYm', PR))
# shifted pair {0, d -+ g} with both classes fixed: rho_pm = rXp / (pair factor?) -> compute directly as XY-type 2-point
# pair (g, d): ratio S({g,d})/(S{g}S{d}) restricted p>=11: generic (1-2/p)/(1-1/p)^2, p | 210k - g: 1/(1-1/p)
A = np.sum(w[1:] * (rYY[1:] - 1))
SXX = 2 * np.sum(w[1:] * (rXX[1:] - 1)); SXY = np.sum(w[1:] * (rXp[1:] - 1 + rXm[1:] - 1))
print('per m: 2A %.4f  SXX %.4f  SXY %.4f  SXX-2SXY+2A %.4f' % (2 * A / m, SXX / m, SXY / m, (SXX - 2 * SXY + 2 * A) / m))
# product-of-pairs approximation for rXX: rYY^2 * rho(d+g) * rho(d-g)
def pair_shift(s):
    # rho for the pair {0, 210k + s g}: generic YY factor, p | 210k + s g -> 1/(1-1/p) ; reuse XY 'lin' structure
    P = PR[PR >= 11]; pf = P.astype(float)
    gen = (1 - 2 / pf) / (1 - 1 / pf) ** 2; lin = 1 / (1 - 1 / pf)
    L = np.full(n, np.log(gen).sum())
    from rho_mine import modinv_vec
    inv = modinv_vec(np.full_like(P, Q), P); k0 = ((-s * g) % P * inv) % P
    cl = np.log(lin / gen); sm = P < n
    for i in np.nonzero(sm)[0]:
        L[int(k0[i])::int(P[i])] += cl[i]
    big = (~sm) & (k0 < n) & (k0 >= 1); np.add.at(L, k0[big], cl[big])
    return np.exp(L)
rp = pair_shift(1); rm = pair_shift(-1)
Ap = np.sum(w[1:] * (rp[1:] - 1)); Am = np.sum(w[1:] * (rm[1:] - 1))
print('per m: 2A+ %.4f 2A- %.4f' % (2 * Ap / m, 2 * Am / m))
prod = rYY ** 2 * rp * rm
print('max |rXX/prod - 1| %.3e, mean %.3e' % (np.max(np.abs(rXX[1:] / prod[1:] - 1)), np.mean(rXX[1:] / prod[1:] - 1)))
print('first-order prediction of SXX: 2(2A + A+ + A-)/m = %.4f; SXY: (2A + A+ + A-)/m = %.4f' % (2 * (2 * A + Ap + Am) / m, (2 * A + Ap + Am) / m))
print('mean over k of rho-1: YY %.2e XX %.2e XYp %.2e shifted %.2e' % (np.mean(rYY[1:] - 1), np.mean(rXX[1:] - 1), np.mean(rXp[1:] - 1), np.mean(rp[1:] - 1)))
