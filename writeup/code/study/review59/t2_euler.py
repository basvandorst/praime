"""Task 2(a, c): the shifted-prime constant C(h) = prod_{l odd} L_l, L_l = 1 + rho_l/(l - 1) (l !| h), 1 - rho_l (l | h),
rho_l = (lambda(l)^2 - 1)/l (split l: lambda = 2 cos 4theta; inert: lambda = 0 so rho = -1/l); the factor at 2 is 1 (ramified).
lambda from own lattice enumeration for l <= LMAX.  Then the separation, in registered errors, of nearby wrong constants
from the holdout readings (readings_holdout.json, A4|shift2 and A4|shift6)."""
import json, math
import numpy as np
LMAX = 10 ** 7

def isqrt_arr(x):
    r = np.floor(np.sqrt(x.astype(np.float64))).astype(np.int64)
    r[r * r > x] -= 1; r[(r + 1) * (r + 1) <= x] += 1
    return r

def lam_upto(M, seg=2 * 10 ** 6):
    out = np.zeros(M)
    for L in range(0, M, seg):
        U = min(L + seg, M); amax = math.isqrt(U - 1); a = np.arange(amax + 1, dtype=np.int64)
        hi = isqrt_arr(U - 1 - a * a)
        lo2 = L - a * a; lo = np.where(lo2 <= 0, 0, isqrt_arr(np.maximum(lo2 - 1, 0)) + 1)
        cnt = np.maximum(hi - lo + 1, 0); tot = int(cnt.sum())
        A = np.repeat(a, cnt); off = np.repeat(np.cumsum(cnt) - cnt, cnt); B = np.repeat(lo, cnt) + (np.arange(tot) - off)
        af = A.astype(float); bf = B.astype(float)
        v = ((af * af - bf * bf) ** 2 - 4 * af * af * bf * bf) * np.where(A > 0, 2.0, 1.0) * np.where(B > 0, 2.0, 1.0)
        out[L:U] += np.bincount(A * A + B * B - L, weights=v, minlength=U - L)
    n = np.arange(M, dtype=float); n[0] = 1
    return out / (4 * n * n)

lam = lam_upto(LMAX + 1)
print('checks: lambda(5) %.4f (2cos4theta = -0.56), lambda(9) %.4f (1), lambda(25) %.4f (lambda5^2-1 = -0.6864), lambda(15) %.1e, lambda(2) %.3f'
      % (lam[5], lam[9], lam[25], lam[15], lam[2]))
s = np.ones(LMAX + 1, bool); s[:2] = False
for q in range(2, math.isqrt(LMAX) + 1):
    if s[q]: s[q * q::q] = False
P = np.nonzero(s)[0]; P = P[P > 2]
rho = (lam[P] ** 2 - 1) / P
def C(h, sign=1.0, drop=(), lmax=LMAX):
    m = P <= lmax; l = P[m]; r = sign * rho[m]
    f = np.where(h % l == 0, 1 - r, 1 + r / (l - 1.0))
    for d in drop: f[l == d] = 1.0
    return float(np.exp(np.sum(np.log(f))))
c2, c6 = C(2), C(6)
print('C(2) = %.5f (registered 0.7791)  C(6) = %.5f (registered 1.2466)  ratio %.4f;  C(2) to l <= 1e5 / 1e6: %.5f %.5f'
      % (c2, c6, c6 / c2, C(2, lmax=10 ** 5), C(2, lmax=10 ** 6)))
R = json.load(open('/home/bas/praime/writeup/data/study/r156/readings_holdout.json'))['stats']
meas = {h: (R['A4|shift%d' % h]['reading'][0], R['A4|shift%d' % h]['sigma_pred']) for h in (2, 6)}
own = json.load(open('t2_shift_out.json'))   # own sub-range sums
print('holdout readings: h=2 %.4f +- %.4f, h=6 %.4f +- %.4f' % (meas[2] + meas[6]))
alts = {'registered': (c2, c6), 'no factor at l=3': (C(2, drop=(3,)), C(6, drop=(3,))),
        'rho sign flipped': (C(2, -1), C(6, -1)), 'h=2 and h=6 swapped (l | h rule off at 3)': (c6, c2),
        'l | h rule ignored (1 + rho/(l-1) everywhere)': (c2, C(1)), 'only l <= 3': (C(2, lmax=3), C(6, lmax=3)),
        'only l <= 100': (C(2, lmax=100), C(6, lmax=100)), 'no factor at l=5': (C(2, drop=(5,)), C(6, drop=(5,))),
        'no factor at l=7': (C(2, drop=(7,)), C(6, drop=(7,))), 'no factors l=5..13': (C(2, drop=(5, 7, 11, 13)), C(6, drop=(5, 7, 11, 13))),
        'all integers (1)': (1.0, 1.0), 'primes (c = 2.0264)': (2.026446, 2.026446)}
print('%-46s %8s %8s %8s %8s %7s' % ('alternative constant', 'C(2)', 'C(6)', 'sep h=2', 'sep h=6', 'chi2'))
for k, (a2, a6) in alts.items():
    s2 = (meas[2][0] - a2 / c2) / meas[2][1]; s6 = (meas[6][0] - a6 / c6) / meas[6][1]
    print('%-46s %8.4f %8.4f %8.2f %8.2f %7.1f' % (k, a2, a6, s2, s6, s2 * s2 + s6 * s6))
