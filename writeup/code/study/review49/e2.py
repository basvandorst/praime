# Explicit formula in the time domain: S_X(t) = -(1/pi) sum_{p^k <= X} b(p^k) p^{-k/2} sin(k t ln p)/k, no zeros used.
# Sigma2_X(l) = var over the registered window grid of S_X(t(x + l)) - S_X(t(x)); compared with the zeros' Sigma2
# (= var of S(t(x+l)) - S(t(x)) with the true S) and with the diagonal (orthogonal-line) sum D_X.
import sys, json, math, numpy as np, cypari2
from scipy.special import loggamma
pari = cypari2.Pari(); pari.allocatemem(300 * 10 ** 6)
Z = np.load(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/data/study/r144/fresh_zeros.npz')
FS = json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']
LS = np.arange(0.5, 10.01, 0.5); STEP = 0.01; LO, HI = 20.0, 400.0
XS = [10, 100, 1000, 10 ** 4, 10 ** 5]


def theta(N, t):
    t = np.asarray(t, float); return 0.5 * t * math.log(N) - t * math.log(2 * math.pi) + np.imag(loggamma(1 + 1j * t))


def run(c):
    N = c['N']; g = np.unique(np.sort(Z[c['id']])); g = g[g > 1e-6]
    b = g[(g >= LO) & (g <= HI)]; ub = theta(N, b) / math.pi
    u0, u1 = ub[0], ub[-1]
    xg = u0 + STEP * np.arange(int((u1 - u0) / STEP) + 2)
    tt = np.linspace(LO - 1, HI + 1, 400001); uu = theta(N, tt) / math.pi
    tg = np.interp(xg, uu, tt)
    E = pari.ellinit(c['ainv']); PMX = max(XS)
    pr = [int(p) for p in pari.primes(pari.primepi(PMX))]
    ap = np.array([int(pari.ellap(E, p)) for p in pr]); pr = np.array(pr)
    bad = np.array([N % p == 0 for p in pr]); lam = ap / np.sqrt(pr)
    terms = []                                             # (n, coefficient, frequency)
    for p, l, bd in zip(pr, lam, bad):
        b0, b1 = 2.0, l; k = 1
        while p ** k <= PMX:
            bk = (l ** k) if bd else b1
            terms.append((p ** k, -bk / (math.pi * k * p ** (k / 2)), k * math.log(p)))
            if not bd: b0, b1 = b1, l * b1 - b0
            k += 1
    terms.sort(); terms = np.array(terms)
    # true S on the grid from all zeros (difference only, so the constant r/2 does not matter)
    uall = theta(N, g) / math.pi
    Strue = np.searchsorted(uall, xg, side='right') - xg
    out = {'id': c['id'], 'N': N, 'g2': c['g2'], 'g3': c['g3'], 'X': {}, 'meas': [], 'check': []}
    def s2(S):
        v = []
        for l in LS:
            m = int(round(l / STEP)); nwin = len(np.arange(u0, u1 - l, STEP))
            d = S[m:m + nwin] - S[:nwin]; v.append(float(d.var()))
        return v
    # consistency: var of the true S difference = counts variance of the registered estimator
    out['meas'] = s2(Strue)
    v = []
    for l in LS:
        x = np.arange(u0, u1 - l, STEP); cnt = np.searchsorted(ub, x + l) - np.searchsorted(ub, x); v.append(float(cnt.var()))
    out['check'] = v
    SX = np.zeros_like(xg); i = 0
    for X in XS:
        j = int(np.searchsorted(terms[:, 0], X, side='right'))
        for a in range(i, j, 64):
            T = terms[a:min(a + 64, j)]
            SX += (T[:, 1][None, :] * np.sin(np.outer(tg, T[:, 2]))).sum(1)
        i = j
        R = Strue - SX
        s_x = s2(SX); s_r = s2(R)
        out['X'][str(X)] = {'s2X': s_x, 's2R': s_r, 'nterms': j}
    return out


if __name__ == '__main__':
    ids = sys.argv[2].split(',')
    res = [run(c) for c in FS if c['id'] in ids]
    json.dump(res, open(sys.argv[1], 'w'))
    for r in res: print(r['id'], 'meas-check', max(abs(a - b) for a, b in zip(r['meas'], r['check'])))
