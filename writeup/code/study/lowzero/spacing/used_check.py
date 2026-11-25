"""Task 2, retrospective (the sets were USED for spacing statistics in rounds 115-116; the numbers below were not known in
this form): the 181 real L-functions of rounds 117/120 (all fundamental discriminants 100 <= |D| <= 400), zeros
100 <= gamma <= 600.  Per L-function Var(s) and <r~>; OLS on [1, I(2|D), I(3|D), I(5|D), I(7|D), I(11|D), I(13|D)];
the same OLS on each variant's predictions (regression-matched).  Also: the closure computed on the REAL base, by re-mapping
the zeros of every L-function with p | D (p = 2, 3, 5) with p's wave (8 phases): its response against the sine-base value."""
import sys, os, math, json
H = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, H + '/lib'); sys.path.insert(0, H)
import numpy as np
import predlib
LZ = os.path.join(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..')), 'data', 'riemann', 'lzeros') + '/'
PR = (2, 3, 5, 7, 11, 13)

def load_R120():
    Z = {}
    for f in ('r117_part0.npz', 'r117_part1.npz', 'r120_part0.npz', 'r120_part1.npz'):
        z = np.load(LZ + f); Z.update({k: np.sort(z[k]) for k in z.files})
    return Z

def unfold(g, q):
    return g / (2 * math.pi) * np.log(q * g / (2 * math.pi * math.e))

def stats(u):
    s = np.diff(u); r = np.minimum(s[1:], s[:-1]) / np.maximum(s[1:], s[:-1])
    return float(np.var(s)), float(np.mean(r)), len(s)

def ols(y, X):
    c, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X @ c; n, k = X.shape
    s2 = res @ res / (n - k); cov = s2 * np.linalg.inv(X.T @ X)
    return c, np.sqrt(np.diag(cov)), math.sqrt(s2)

if __name__ == '__main__':
    Z = load_R120(); P = predlib.Predictor()
    rows = []
    for key in sorted(Z, key=lambda k: (abs(int(k[1:])), int(k[1:]))):
        D = int(key[1:]); q = abs(D); g = Z[key]; g = g[(g >= 100) & (g <= 600)]
        v, r, n = stats(unfold(g, q)); pr = P.predict(D, 100, 600, law2=True)
        rows.append(dict(D=D, var=v, r=r, n=n, pred=pr))
    print('%d L-functions, mean spacings %.0f, mean L %.3f' % (len(rows), np.mean([x['n'] for x in rows]), np.mean([x['pred']['meanL'] for x in rows])))
    X = np.array([[1.0] + [float(abs(x['D']) % p == 0) for p in PR] for x in rows])
    VARS = ('law', 'law2', 'slow', 'bblm', 'ignored'); out = {'rows': rows}
    for stat in ('var', 'r'):
        y = np.array([x[stat] for x in rows]); c, se, sres = ols(y, X)
        print('\n[%s] data: level (mean) %.5f; OLS residual sd %.5f' % (stat, y.mean(), sres))
        print('  coefficient     const     p=2      p=3      p=5      p=7      p=11     p=13')
        print('  data          ' + ' '.join('%+.5f' % v for v in c))
        print('  se            ' + ' '.join(' %.5f' % v for v in se))
        res = {'data': c.tolist(), 'se': se.tolist()}
        for vn in VARS:
            yp = np.array([x['pred'][vn][stat] for x in rows]); cp, _, _ = ols(yp, X)
            res[vn] = cp.tolist()
            z = (c[1:4] - cp[1:4]) / se[1:4]
            print('  %-12s  ' % vn + ' '.join('%+.5f' % v for v in cp) + '   z(2,3,5) = ' + ' '.join('%+.1f' % v for v in z) + '   mean(data-pred) %+.5f' % (y - yp).mean())
        out[stat] = res
    # closure on the real base: re-map the zeros of L-functions with p | D by p's wave
    print('\nclosure on the real base (zeros with p missing, re-mapped by p wave, 8 phases): change "missing minus present"')
    Lg = []; emp = {}
    for p in (2, 3, 5):
        dvs, drs, dvp, drp = [], [], [], []
        for x in rows:
            D = x['D']; q = abs(D)
            if q % p: continue
            g = Z['D%d' % D]; g = g[(g >= 90) & (g <= 610)]
            tt = np.linspace(85, 615, 20001); uu = unfold(tt, q)
            u = unfold(g, q); acc = []
            for phi in 2 * math.pi * np.arange(8) / 8:
                xm = u.copy()
                for _ in range(30):
                    tx = np.interp(xm, uu, tt)
                    dN = np.angle(1 - p ** -0.5 * np.exp(1j * (tx * math.log(p) + phi))) / math.pi
                    xm = u - dN
                tx = np.interp(xm, uu, tt); sel = (tx >= 100) & (tx <= 600)
                acc.append(stats(xm[sel])[:2])
            acc = np.mean(acc, 0)
            dvs.append(x['var'] - acc[0]); drs.append(x['r'] - acc[1])
            # sine-base closure at this L-function's mean L, single prime
            L = x['pred']['meanL']; dvp.append(float(P._s('dv', p)(L))); drp.append(float(P._s('dr', p)(L)))
        dvs, drs, dvp, drp = map(np.array, (dvs, drs, dvp, drp))
        emp[p] = dict(n=len(dvs), dvar_real=dvs.mean(), dvar_sine=dvp.mean(), dr_real=drs.mean(), dr_sine=drp.mean(), dvar_real_se=dvs.std() / math.sqrt(len(dvs)), dr_real_se=drs.std() / math.sqrt(len(drs)))
        print('  p=%d (n=%d): dVar real base %+.5f +- %.5f   sine base %+.5f  (ratio %.3f) | d<r~> real %+.5f +- %.5f  sine %+.5f (ratio %.3f)'
              % (p, len(dvs), dvs.mean(), dvs.std() / math.sqrt(len(dvs)), dvp.mean(), dvs.mean() / dvp.mean(), drs.mean(), drs.std() / math.sqrt(len(drs)), drp.mean(), drs.mean() / drp.mean()))
    out['realbase'] = emp
    json.dump(out, open('used_check.json', 'w'), indent=1, default=float)
