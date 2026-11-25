"""dev_analysis.py: development characters only (the 8 of round 149's dev list and the 24 of dev_extra.json).
(1) per-character observables (means over the 143 intervals of length 2 pi; and the variance of M over them);
(2) the indicator regression on p = 2, 3, 5, 7 with intercept and ln q, residual sd per character (the noise);
(3) the within-character standard error of each mean (from the interval values, with lag correlations up to 3);
(4) surgery (the data-driven closure): remove p's Euler factor from a character where p is good, add it (both signs)
    where p is bad; the mean contrast 'p missing minus p present' per observable;
(5) zeta (L about 4.3) for reference.
Writes dev_obs.json."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import obslib

HERE = os.path.dirname(os.path.abspath(__file__))
DEV8 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', 'data', 'study', 'r149', 'dev_values')
DEVX = os.path.join(HERE, 'dev_extra_values')
PR = (2, 3, 5, 7)
KEYS = ['M', 'lnZ_0.5', 'lnZ_1', 'lnZ_1.5', 'lnZ_2', 'lnZ_3', 'mean_ln']


def kron(D, p):
    if D % p == 0:
        return 0
    if p == 2:
        return 1 if D % 8 in (1, 7) else -1
    r = pow(D % p, (p - 1) // 2, p)
    return 1 if r == 1 else -1


def load(D):
    for d in (DEV8, DEVX):
        f = os.path.join(d, 'absL_%d.npy' % D)
        if os.path.exists(f):
            return np.load(f)
    raise FileNotFoundError(D)


def char_stats(v):
    o = obslib.observables(v)
    s = {k: float(o[k].mean()) for k in KEYS}
    s['varM'] = float(o['M'].var(ddof=1))
    s['Mgrid'] = float(o['Mgrid'].mean())
    # within-character se of the means: var/n * (1 + 2 sum_{l<=3} (1 - l/n) r_l)
    for k in KEYS:
        x = o[k] - o[k].mean(); n = len(x)
        r = [np.dot(x[:-l], x[l:]) / np.dot(x, x) for l in (1, 2, 3)]
        f = 1 + 2 * sum((1 - l / n) * r[l - 1] for l in (1, 2, 3))
        s['se_' + k] = float(np.sqrt(max(f, 0.05) * x.var(ddof=1) / n))
    return s, o


def ols(y, X):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = res @ res / dof
    cov = s2 * np.linalg.inv(X.T @ X)
    return beta, np.sqrt(np.diag(cov)), np.sqrt(s2)


def main():
    dev8 = [1001, -1003, 1005, -1048, 1068, -1195, 1253, -1267]
    devx = json.load(open(os.path.join(HERE, 'dev_extra.json')))['dev_extra']
    Ds = dev8 + devx
    t = 100 + obslib.THETA * (np.arange(obslib.NINT) + 0.5)
    rows, surg = [], {p: {k: [] for k in KEYS + ['varM']} for p in PR}
    grid = obslib.grid()
    for D in Ds:
        v = load(D)
        s, o = char_stats(v)
        q = abs(D)
        s['D'] = D; s['lnq'] = float(np.log(q)); s['Lbar'] = float(np.mean(np.log(q * t / (2 * np.pi))))
        s['bad'] = [p for p in PR if D % p == 0]
        rows.append(s)
        for p in PR:
            c = kron(D, p)
            if c != 0:      # remove p: |L| * |1 - chi(p) p^{-1/2 - it}|
                vr = v * np.abs(1 - c * p ** -0.5 * np.exp(-1j * grid * np.log(p)))
                sr, _ = char_stats(vr)
                for k in KEYS + ['varM']:
                    surg[p][k].append(sr[k] - s[k])
            else:           # add p with chi(p) = +1 and -1
                for e in (1, -1):
                    va = v / np.abs(1 - e * p ** -0.5 * np.exp(-1j * grid * np.log(p)))
                    sa, _ = char_stats(va)
                    for k in KEYS + ['varM']:
                        surg[p][k].append((s[k] - sa[k]) / 2 if e == 1 else (s[k] - sa[k]) / 2)
        print('D=%6d L=%.2f bad=%-8s ' % (D, s['Lbar'], s['bad']) +
              ' '.join('%s %.3f(%.3f)' % (k, s[k], s['se_' + k]) for k in ['M', 'lnZ_1', 'lnZ_3']) +
              ' varM %.3f' % s['varM'], flush=True)
    # regression
    y = {k: np.array([r[k] for r in rows]) for k in KEYS + ['varM']}
    X = np.column_stack([np.ones(len(rows)), [r['lnq'] - np.mean([r['lnq'] for r in rows]) for r in rows]] +
                        [[1.0 if p in r['bad'] else 0.0 for r in rows] for p in PR])
    print('\ngroup sizes', {p: int(sum(1 for r in rows if p in r['bad'])) for p in PR}, 'n =', len(rows))
    reg = {}
    for k in KEYS + ['varM']:
        b, se, sres = ols(y[k], X)
        within = np.sqrt(np.mean([r['se_' + k] ** 2 for r in rows])) if k != 'varM' else float('nan')
        reg[k] = {'coef': b[2:].tolist(), 'se': se[2:].tolist(), 'resid_sd': float(sres), 'within_se': float(within)}
        print('%-8s resid sd %.4f (within-character se %.4f)  ' % (k, sres, within) +
              ' '.join('c%d %+.4f(%.4f)' % (p, b[2 + i], se[2 + i]) for i, p in enumerate(PR)))
    # surgery means
    print('\nsurgery (p missing minus present; mean over characters, se):')
    sg = {}
    for p in PR:
        sg[p] = {k: [float(np.mean(surg[p][k])), float(np.std(surg[p][k], ddof=1) / np.sqrt(len(surg[p][k])))]
                 for k in KEYS + ['varM']}
        print('p=%d n=%d ' % (p, len(surg[p]['M'])) + ' '.join('%s %+.4f(%.4f)' % (k, *sg[p][k]) for k in KEYS + ['varM']))
    # zeta
    sz, _ = char_stats(np.load(os.path.join(DEV8, 'absL_1.npy')))
    print('\nzeta (L 2.8-5.0): ' + ' '.join('%s %.3f' % (k, sz[k]) for k in KEYS) + ' varM %.3f' % sz['varM'])
    json.dump({'rows': rows, 'reg': reg, 'surgery': sg, 'zeta': sz}, open(os.path.join(HERE, 'dev_obs.json'), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
