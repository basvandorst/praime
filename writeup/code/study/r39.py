"""Round 39 (preregister/round39.txt; not blind): the under-dispersion of the residue tables of round 13 derived from
the long-range rigidity of the Hardy-Littlewood correlations. For row a (mod 210) and gap g:
  X = gap-g events (p, p + g consecutive) with p = a, Y = primes = a, P = the law's share P_a(g),
  E[(X - P Y)^2] = Var X - 2 P Cov(X, Y) + P^2 Var Y,
with Var Y from pairs {0, d}, Var X from {0, g, d, d + g}, Cov from {0, g, +-d}, d = 210 k (factors p | 210 cancel).
For g >= 6 the emptiness of two interiors 210 k apart is taken as independent (heuristic).
usage: r39.py [sub]"""
import os, json, math, pickle
import numpy as np
from scipy.signal import fftconvolve
import common as C
import theory_r3 as R3
import r7

Q = 210
PR = C.small_primes(20_000_000).astype(np.int64)
P11 = PR[PR >= 11]


def logconst(f_none, g):
    """sum over p >= 11 of log f(p, none), with p | g using its own generic factor; tail beyond 2e7 below 1e-8."""
    p = P11.astype(np.float64)
    v = np.log(f_none(p, (g % P11) == 0))
    return float(v.sum())


def rho(n, g, kind):
    """rho(k), k = 0..n-1 (index 0 unused), for kind 'YY', 'XX', 'XYp' (d = +210k), 'XYm' (d = -210k)."""
    def f(p, case, pg):
        # local factor at p >= 11 of S(set)/prod S(parts); case: 'none', 'k' (p | k), 'lin' (p | 210k -+ g)
        if kind == 'YY':
            return {'none': (1 - 2 / p) / (1 - 1 / p) ** 2, 'k': 1 / (1 - 1 / p)}[case]
        if kind == 'XX':
            if pg:
                return {'none': (1 - 2 / p) / (1 - 1 / p) ** 2, 'k': 1 / (1 - 1 / p)}[case]
            return {'none': (1 - 4 / p) / (1 - 2 / p) ** 2, 'k': 1 / (1 - 2 / p), 'lin': (1 - 3 / p) / (1 - 2 / p) ** 2}[case]
        if pg:                                                   # XY with p | g
            return {'none': (1 - 2 / p) / (1 - 1 / p) ** 2, 'k': 1 / (1 - 1 / p)}[case]
        return {'none': (1 - 3 / p) / ((1 - 2 / p) * (1 - 1 / p)), 'k': 1 / (1 - 1 / p), 'lin': 1 / (1 - 1 / p)}[case]

    pgv = lambda p: (g % p) == 0
    L = np.full(n, logconst(lambda p, pg: np.where(pg, f(p, 'none', True) if kind != 'YY' else f(p, 'none', False),
                                                   f(p, 'none', False)), g if kind != 'YY' else 1))
    # p | k, all primes p >= 11 below n
    for p in P11[P11 < n]:
        p = int(p)
        pg = kind != 'YY' and g % p == 0
        L[p::p] += math.log(f(p, 'k', pg) / f(p, 'none', pg))
    # p | 210 k - g (XX, XYp) and p | 210 k + g (XX, XYm), p does not divide g
    sgns = {'YY': (), 'XX': (-1, 1), 'XYp': (-1,), 'XYm': (1,)}[kind]
    k = np.arange(n, dtype=np.int64)
    for s in sgns:
        cof = np.abs(Q * k + s * g)
        for q in (2, 3, 5, 7):
            for _ in range(40):
                m = cof % q == 0
                if not m.any():
                    break
                cof[m] //= q
        lim = int(math.isqrt(int(Q * n + g))) + 1
        for p in P11[P11 <= lim]:
            p = int(p)
            if g % p == 0:
                k0 = 0                                           # p | 210k +- g iff p | k: already in the p | k term
                idx = np.arange(0, n, p)
            else:
                k0 = (-s * g * pow(Q, -1, p)) % p
                idx = np.arange(k0, n, p)
                L[idx] += math.log(f(p, 'lin', False) / f(p, 'none', False))
            for _ in range(40):
                m = cof[idx] % p == 0
                if not m.any():
                    break
                cof[idx[m]] //= p
        big = (cof > 1)
        pb = cof[big].astype(np.float64)
        L[big] += np.log(f(pb, 'lin', False) / f(pb, 'none', False))
    return np.exp(L)


def window(L0, W, gaps, Pa_rows):
    n = W // Q
    # density profile of one class: r_j = (210/48)/ln(m_j)
    r = (Q / 48) / np.log(L0 + Q * np.arange(n) + Q / 2)
    w = fftconvolve(r, r[::-1])[n - 1:]                     # w[k] = sum_j r_j r_{j+k}
    m = float(r.sum()); s1 = float(np.sum(r * (1 - r)))
    rYY = rho(n, 1, 'YY')
    varY = s1 + 2 * float(np.sum(w[1:] * (rYY[1:] - 1)))
    out = {'n': n, 'm': m, 'vY': varY / m}
    for g in gaps:
        rXX = rho(n, g, 'XX'); rXp = rho(n, g, 'XYp'); rXm = rho(n, g, 'XYm')
        SXX = 2 * float(np.sum(w[1:] * (rXX[1:] - 1)))
        SXY = float(np.sum(w[1:] * ((rXp[1:] - 1) + (rXm[1:] - 1))))
        vals, nulls = [], []
        for P in Pa_rows[g]:
            if P <= 0:
                continue
            varX = P * m - P * P * float(np.sum(r * r)) + P * P * SXX
            cov = P * s1 + P * SXY
            ER2 = varX - 2 * P * cov + P * P * varY
            vals.append(ER2 / (P * m)); nulls.append(1 - P)
        out[str(g)] = {'pred': float(np.mean(vals)), 'null': float(np.mean(nulls)), 'rows': len(vals),
                       'per_row': vals, 'P': [float(x) for x in Pa_rows[g] if x > 0]}
        print(g, 'pred %.3f null %.3f rows %d' % (out[str(g)]['pred'], out[str(g)]['null'], len(vals)), flush=True)
    return out


def main():
    cache = pickle.load(open(r7.CACHE7, 'rb'))
    U = r7.units(Q)
    res = {}
    for win in ('train', 'A', 'B', 'C'):
        L0, W = R3.WIN[win]
        T = r7.trans(cache, Q, L0, W); Pa = T.sum(2).T; Pa = Pa / Pa.sum(1, keepdims=True)
        gaps = list(range(2, 39, 2))
        Pa_rows = {g: [float(Pa[i, g // 2 - 1]) if math.gcd(int(a) + g, Q) == 1 else 0.0 for i, a in enumerate(U)] for g in gaps}
        print(win, flush=True)
        res[win] = window(L0, W, gaps, Pa_rows)
        res[win]['P_pooled40'] = [float(x) for x in Pa[:, 19]]
        print(win, 'vY %.4f' % res[win]['vY'], flush=True)
    json.dump(res, open(os.path.join(R3.OUT, 'r39.json'), 'w'), indent=1)


def implied(win_res, P40):
    num = den = 0.0
    for g in range(2, 39, 2):
        d = win_res[str(g)]
        num += sum(d['per_row']); den += sum(1 - p for p in d['P'])
    P40 = [p for p in P40 if p > 0]
    num += sum(1 - p for p in P40); den += sum(1 - p for p in P40)
    return num / den


def sub_windows():
    """Addendum (measured values known from round 13): the sub-windows of train, the first 9e6 or 9e5 numbers of each
    tenth of [1e7, 1e8), as in r13c.py."""
    import contextlib, io
    cache = pickle.load(open(r7.CACHE7, 'rb'))
    U = r7.units(Q); L, W = R3.WIN['train']
    out = json.load(open(os.path.join(R3.OUT, 'r39.json')))
    for win in ('train', 'A', 'B', 'C'):
        out[win]['implied_chi2_dof'] = implied(out[win], out[win]['P_pooled40'])
    for Wsub in (9 * 10 ** 6, 9 * 10 ** 5):
        vals = []
        for k in range(10):
            Ls = L + k * 9 * 10 ** 6
            T = r7.trans(cache, Q, Ls, Wsub); Pa = T.sum(2).T; Pa = Pa / Pa.sum(1, keepdims=True)
            Pa_rows = {g: [float(Pa[i, g // 2 - 1]) if math.gcd(int(a) + g, Q) == 1 else 0.0 for i, a in enumerate(U)]
                       for g in range(2, 39, 2)}
            with contextlib.redirect_stdout(io.StringIO()):
                wr = window(Ls, Wsub, list(range(2, 39, 2)), Pa_rows)
            vals.append(implied(wr, [float(x) for x in Pa[:, 19]]))
        out['train_sub_W=%d' % Wsub] = {'implied_chi2_dof': float(np.mean(vals)), 'values': vals}
        print('W', Wsub, round(float(np.mean(vals)), 3), flush=True)
    json.dump(out, open(os.path.join(R3.OUT, 'r39.json'), 'w'), indent=1)


if __name__ == '__main__':
    import sys
    sub_windows() if sys.argv[1:] == ['sub'] else main()
