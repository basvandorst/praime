"""consec_bks.py: the block-dispersion (column) law for the non-abelian and abelian contrasts, and its comparison.

For a primitive L-function of degree d and conductor q, the variance of sum_{x < n <= x + H} Lambda_F(n), divided by its
diagonal H ln x, is (Montgomery-Soundararajan for d = 1; Bui-Keating-Smith 2016, Theorem C1, for d >= 2, integrated form
converted to a local one by d/dX):
    R(x, H) = [d ln(x/H) + ln q + d B] / ln x,   B = 1 - gamma - ln 2 pi,     for H above the transition (d >= 2),
    R(x, H) = 1 (H0: no off-diagonal term at fixed shifts; the diagonal alone)   below it,
and the model used here is R = min(1, ...) (a kink; the true crossover is smooth).  For d = 1 the first line holds for all
H >= ln x (no regime change).  A function of the label that mixes several L-functions gets the |c|^2-weighted mean.
Blocks are pooled as in consec_analyze (sum T^2 / sum V), so the prediction for a height bin is the V-weighted mean over
its blocks, V ~ H / ln x.

Components (d, q, weight): conductor of the Artin / Hecke L-function behind each function (see consec.md 2.4).
usage: consec_bks.py READINGS.json   (prints measured R against the law per height bin; --json writes predictions)"""
import sys, json, math
import numpy as np
GAM = 0.5772156649015329; B = 1 - GAM - math.log(2 * math.pi)
C1 = (8 / math.pi) / 3                      # E[sign(cos t) * 2 cos t] under Sato-Tate = 8/(3 pi)
COMP = {   # scheme -> list per function (f0, f1, ...) of [(d, q, weight)], then the abelian sign 'chi'
    'C23': [[(2, 23, 1.0)]], 'C31': [[(2, 31, 1.0)]], 'C47': [[(2, 47, 1.0)], [(2, 47, 1.0)]], 'C39': [[(2, 39, 1.0)]],
    'G': [[(2, 4, 1.0)], [(2, 4, 1.0)]],
    'E37s': [[(2, 37, C1 ** 2), (4, 37 ** 3, 1 - C1 ** 2)], [(2, 37, 1.0)], [(3, 37 ** 2, 1.0)]],
    'E37m2': [[(2, 148, 1.0)]],
}
CHIQ = {'C23': 23, 'C31': 31, 'C47': 47, 'C39': 39, 'G': 4, 'E37m2': 37, 'E37m3': 3}

def r_one(x, H, d, q):
    v = (d * math.log(x / H) + math.log(q) + d * B) / math.log(x)
    return v if d == 1 else min(1.0, v)

def predict(comp, H, b0, b1, nsub=400):
    """V-weighted mean of R over blocks of length H in [b0, b1) (integration over x with weight 1/ln x)."""
    b0 = max(b0, 2 * H, 1000); xs = np.linspace(b0, b1, nsub + 1); xm = 0.5 * (xs[1:] + xs[:-1]); w = 1 / np.log(xm)
    r = np.array([sum(wt * r_one(x, H, d, q) for (d, q, wt) in comp) for x in xm])
    return float((r * w).sum() / w.sum())

if __name__ == '__main__':
    R = json.load(open(sys.argv[1])); out = {}
    for nm, fl in COMP.items():
        if nm not in R: continue
        for key, ent in R[nm]['disp'].items():
            H, b0, b1 = key.split('_'); H = int(H[1:]); b0 = int(b0); b1 = int(b1)
            row = []
            for i, comp in enumerate(fl + [[(1, CHIQ.get(nm, 1), 1.0)]] if nm in CHIQ else fl):
                f = 'f%d' % i if i < len(fl) else 'chi'
                if f not in ent or ent[f]['R'] != ent[f]['R']: continue
                p = predict(comp, H, b0, b1); m, e = ent[f]['R'], ent[f]['err']
                row.append('%s %.4f+-%.4f law %.4f z %+.1f' % (f, m, e, p, (m - p) / e if e > 0 else 0))
                out['%s|%s|%s' % (nm, key, f)] = dict(measured=m, err=e, law=p)
            print('%-6s %-26s ' % (nm, key) + '  '.join(row))
    if '--json' in sys.argv: json.dump(out, open(sys.argv[sys.argv.index('--json') + 1], 'w'), indent=0)
