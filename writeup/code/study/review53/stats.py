"""Readings recomputed from my own integrals (meas.npz) and the frozen predictions (numbers only: Pfull etc.).
Own jackknife (delete-one-character), own OLS."""
import json, math, numpy as np, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
m = np.load(os.path.join(HERE, 'meas.npz'))
P = json.load(open(os.path.dirname(os.path.abspath(__file__))+'/../../../../writeup/code/study/lowzero/moment3/predictions_triple.json'))
S = P['sets']; R = S['reserved']; tri = [tuple(t) for t in S['triples']]; pairs = [tuple(p) for p in S['pairs']]
TT = {tuple(r['key']): r for r in P['triples']}; PP = {tuple(r['key']): r for r in P['pairs']}; SS = {r['key'][0]: r for r in P['singles']}
M1 = m['single_trap']; M2 = dict(zip(pairs, m['pair_trap'])); M3 = m['triple_trap']; M4 = m['single4_trap']


def jk(stat, tuples, chars=R):
    full = stat(np.arange(len(tuples)))
    reps = []
    for c in chars:
        keep = np.array([i for i, tp in enumerate(tuples) if c not in tp])
        reps.append(stat(keep))
    reps = np.array(reps); n = len(chars)
    return full, math.sqrt((n - 1) / n * np.sum((reps - reps.mean(0)) ** 2))


def slope(x, y):
    x = x - x.mean(); y = y - y.mean(); return float((x * y).sum() / (x * x).sum())


if __name__ == '__main__':
    r1 = M1 / np.array([SS[D]['P1'] for D in R]); r4 = M4 / np.array([SS[D]['P4'] for D in R])
    print('singles r1 mean %.4f se %.4f | fourth %.4f se %.4f' % (r1.mean(), r1.std(ddof=1) / math.sqrt(56), r4.mean(), r4.std(ddof=1) / math.sqrt(56)))
    r2 = np.array([M2[p] / PP[p]['full'] for p in pairs])
    print('pairs mean %.4f +- %.4f' % jk(lambda k: r2[k].mean(), pairs))
    pf = np.array([TT[t]['full'] for t in tri]); r3 = M3 / pf
    rp = {p: M2[p] / PP[p]['full'] for p in pairs}
    cpair = np.array([rp[(a, b)] + rp[(a, c)] + rp[(b, c)] - 3 for a, b, c in tri])
    zlog = np.array([math.log(rp[(a, b)]) + math.log(rp[(a, c)]) + math.log(rp[(b, c)]) for a, b, c in tri])
    print('triples raw %.4f +- %.4f' % jk(lambda k: r3[k].mean(), tri))
    print('triples CV  %.4f +- %.4f' % jk(lambda k: (r3[k] - cpair[k]).mean(), tri))
    print('mean cpair over triples %.4f +- %.4f' % jk(lambda k: cpair[k].mean(), tri))
    PC = {c: np.array([TT[t][c] for t in tri]) for c in ['full', 'diag', 'kirk', 'no3', 'ind', 'lead']}
    lm = np.log(M3)
    for c in ['no3', 'kirk', 'diag', 'ind', 'lead']:
        xf = np.log(PC['full']); xc = np.log(PC[c])
        g = jk(lambda k: slope(xf[k] - xc[k], lm[k] - zlog[k] - xc[k]), tri)
        g0 = jk(lambda k: slope(xf[k] - xc[k], lm[k] - xc[k]), tri)
        print('gamma_%-4s CV %.4f +- %.4f raw %.4f +- %.4f' % (c, g[0], g[1], g0[0], g0[1]))
    print('slope raw %.4f +- %.4f' % jk(lambda k: slope(np.log(pf[k]), lm[k]), tri))
    print('slope CV  %.4f +- %.4f' % jk(lambda k: slope(np.log(pf[k]), lm[k] - zlog[k]), tri))
