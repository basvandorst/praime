"""Independent recomputation of the round-146 variance reading from the stored per-zero values.
Own code: reads the npz/json data only."""
import os, numpy as np, json, sys

R = os.path.dirname(os.path.abspath(__file__))+'/../../../../'
PRED = json.load(open(R + 'writeup/code/study/lowzero/varlp/predictions_varlp.json'))
SETS = {'C': (R + 'writeup/data/study/r146/vals_C.npz', R + 'writeup/data/riemann/lzeros/r134_C.npz'),
        'Bp': (R + 'writeup/data/study/r146/vals_Bp.npz', R + 'writeup/data/riemann/lzeros/r138_Bp.npz')}
PR = [2, 3, 5, 7, 11, 13]


def load(lo=600., hi=1000., detrend='L'):
    rows = []
    for s, (fv, fz) in SETS.items():
        V = np.load(fv, allow_pickle=True)
        Z = np.load(fz, allow_pickle=True)
        preds = {e['D']: e for e in PRED[s]}
        for D, e in preds.items():
            g = V['D%d__gamma' % D]; la = V['D%d__logabs' % D]; aL = V['D%d__absL' % D]
            z = Z['D%d' % D]; z = z[(z >= 600) & (z <= 1000)]
            same = (len(z) == len(g)) and np.allclose(z, g, atol=1e-9, rtol=0)
            q = abs(D)
            Lj = np.log(q * g / (2 * np.pi))
            m = (g >= lo) & (g <= hi)
            rows.append(dict(set=s, D=D, q=q, g=g[m], la=la[m], Lj=Lj[m], ratio=np.max(aL / np.exp(la)),
                             same=same, pred=e))
    return rows


def stat(r, detrend='L', trim=0.0):
    y = r['la'] - (r['Lj'] if detrend == 'L' else np.log(r['Lj']))
    if trim > 0:
        c = np.quantile(r['la'], trim)
        y = y[r['la'] > c]
    return np.var(y, ddof=1), np.mean(y)


def design(rows, primes=PR, extra=None):
    sets = sorted(set(r['set'] for r in rows))
    n = len(rows)
    cols = [np.ones(n)]
    for s in sets[1:]:
        cols.append(np.array([1.0 if r['set'] == s else 0.0 for r in rows]))
    cols.append(np.array([np.mean(np.log(r['Lj'])) for r in rows]))
    if extra is not None:
        cols += extra
    nn = len(cols)
    for p in primes:
        cols.append(np.array([1.0 if r['q'] % p == 0 else 0.0 for r in rows]))
    X = np.array(cols).T
    X[:, 1:] -= X[:, 1:].mean(0)
    return X, nn


def ols(X, y):
    return np.linalg.lstsq(X, y, rcond=None)[0]


if __name__ == '__main__':
    rows = load()
    bad = [r['D'] for r in rows if not r['same']]
    print('L-functions', len(rows), 'height mismatch', bad, 'max |L|/|L\'|', max(r['ratio'] for r in rows))
    for excl in [(197, -395), ()]:
        kept = [r for r in rows if r['D'] not in excl]
        for detrend in ['L', 'lnL']:
            v = np.array([stat(r, detrend)[0] for r in kept])
            X, nn = design(kept)
            b = ols(X, v)
            pr = {k: ols(X, np.array([r['pred'][k] for r in kept])) for k in ['cum2_3pt', 'cum2_ldm', 'cum2_lead']}
            rng = np.random.default_rng(20261126)
            n = len(kept); B = []
            for _ in range(2000):
                ix = rng.integers(0, n, n)
                Xb = X[ix].copy(); Xb[:, 1:] -= Xb[:, 1:].mean(0)
                B.append(ols(Xb, v[ix]))
            B = np.array(B); se = B.std(0, ddof=1)
            print('\nexclude', excl, 'detrend', detrend, 'n', len(kept), 'zeros', sum(len(r['g']) for r in kept),
                  'groups', [sum(1 for r in kept if r['q'] % p == 0) for i, p in enumerate(PR)])
            print('fit ', ' '.join('%+.5f' % x for x in b[nn:]))
            print('se  ', ' '.join('%.5f' % x for x in se[nn:]))
            for k, c in pr.items():
                print('%-9s' % k, ' '.join('%+.5f' % x for x in c[nn:]), ' z', ' '.join('%+.1f' % x for x in (b[nn:] - c[nn:]) / se[nn:]))
            l3 = pr['cum2_3pt'][nn:nn + 3]; ll = pr['cum2_lead'][nn:nn + 3]; ld = pr['cum2_ldm'][nn:nn + 3]
            def amp(bb, z0):
                return np.sum((bb - z0) * (l3 - z0)) / np.sum((l3 - z0) ** 2)
            ga = amp(b[nn:nn + 3], ll); gs = np.std([amp(x[nn:nn + 3], ll) for x in B], ddof=1)
            gl = amp(b[nn:nn + 3], ld); gls = np.std([amp(x[nn:nn + 3], ld) for x in B], ddof=1)
            print('gamma %.4f +- %.4f   gamma_L %.4f +- %.4f' % (ga, gs, gl, gls))
            # joint chi2 over p=2,3,5,7 with bootstrap covariance
            d = b[nn:nn + 4] - pr['cum2_3pt'][nn:nn + 4]
            Cv = np.cov(B[:, nn:nn + 4].T)
            chi2 = d @ np.linalg.solve(Cv, d)
            from scipy.stats import chi2 as C2
            print('joint chi2(4) vs 3pt = %.2f  p = %.3f ; corr matrix:' % (chi2, C2.sf(chi2, 4)))
            print(np.round(np.corrcoef(B[:, nn:nn + 4].T), 2))
            # common-scale: fit - 3pt = k * 3pt ? (fractional deficit)
            w = np.linalg.solve(Cv, np.ones(4))
            frac = d / pr['cum2_3pt'][nn:nn + 4]
            print('fractional deviations', np.round(frac, 3))
