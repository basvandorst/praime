"""dev_closure.py: why surgery (removing p's Euler factor from a character with p good) gives a smaller contrast
than the real contrast.  Development characters only.
(a) Fourier coefficient of ln|L(1/2+it)| at frequency ln p over [100, 1000], against the theorem value chi(p) p^{-1/2}/2;
(b) within-interval variance of ln|L| (clipped at ln 1e-6) for p | D against p not | D, and of R = ln|L| - Y_p;
(c) surgery split: removal (p good) and addition (p bad) separately, for lnZ_0.5, lnZ_1, M."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import obslib
from dev_analysis import load, kron, char_stats

HERE = os.path.dirname(os.path.abspath(__file__))
PR = (2, 3, 5, 7)


def within_var(x):
    out = []
    for idx, w in obslib.weights():
        xi = x[np.clip(idx, 0, obslib.NPT - 1)]
        m = np.dot(w, xi) / obslib.THETA
        out.append(np.dot(w, (xi - m) ** 2) / obslib.THETA)
    return float(np.mean(out))


def main():
    dev8 = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
    Ds = dev8 + json.load(open(os.path.join(HERE, 'dev_extra.json')))['dev_extra']
    g = obslib.grid()
    wts = np.full(len(g), obslib.H); wts[0] = wts[-1] = obslib.H / 2
    T = g[-1] - g[0]
    res = {p: {'four': [], 'vw_good': [], 'vw_R': [], 'vw_bad': [], 'rem': {}, 'add': {}} for p in PR}
    for D in Ds:
        v = load(D)
        lnL = np.log(np.maximum(v, 1e-6))
        vw = within_var(lnL)
        s0, _ = char_stats(v)
        for p in PR:
            c = kron(D, p)
            Yp = -np.log(np.abs(1 - (c if c else 1) * p ** -0.5 * np.exp(-1j * g * np.log(p))))
            if c != 0:
                f = np.sum(wts * lnL * np.exp(1j * g * np.log(p))) / T
                res[p]['four'].append(float(np.real(f) / (c * p ** -0.5 / 2)))
                res[p]['vw_good'].append(vw)
                res[p]['vw_R'].append(within_var(lnL - Yp))
                sr, _ = char_stats(v * np.exp(-Yp))
                for k in ('lnZ_0.5', 'lnZ_1', 'M'):
                    res[p]['rem'].setdefault(k, []).append(sr[k] - s0[k])
            else:
                res[p]['vw_bad'].append(vw)
                for e in (1, -1):
                    Ye = -np.log(np.abs(1 - e * p ** -0.5 * np.exp(-1j * g * np.log(p))))
                    sa, _ = char_stats(v * np.exp(Ye))
                    for k in ('lnZ_0.5', 'lnZ_1', 'M'):
                        res[p]['add'].setdefault(k, []).append(s0[k] - sa[k])
    for p in PR:
        r = res[p]
        Yw = None
        print('p=%d  Fourier coef / theorem value: mean %.3f sd %.3f (n=%d)' % (p, np.mean(r['four']), np.std(r['four']),
              len(r['four'])))
        print('      within-interval var of ln|L|: p good %.4f, p bad %.4f, diff %.4f;  R = ln|L| - Y_p (p good) %.4f'
              % (np.mean(r['vw_good']), np.mean(r['vw_bad']), np.mean(r['vw_good']) - np.mean(r['vw_bad']),
                 np.mean(r['vw_R'])))
        for k in ('lnZ_0.5', 'lnZ_1', 'M'):
            print('      %-8s removal %+.4f (%.4f, n=%d)   addition %+.4f (%.4f, n=%d)' % (
                k, np.mean(r['rem'][k]), np.std(r['rem'][k]) / np.sqrt(len(r['rem'][k])), len(r['rem'][k]),
                np.mean(r['add'][k]), np.std(r['add'][k]) / np.sqrt(len(r['add'][k])), len(r['add'][k])))


if __name__ == '__main__':
    main()
