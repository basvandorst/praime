# Review analysis of round 144 with own Sigma2 (s2_own.json) and own law (lawout/*.json)
import json, glob, math, sys, numpy as np, cypari2
from scipy import stats
pari = cypari2.Pari()
D = __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/'
S2 = json.load(open(D + 's2_own.json'))
FS = {c['id']: c for c in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_set.json'))['curves']}
P = {c['id']: c for c in json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_predictions.json'))['curves']}
GUE = np.array(json.load(open(__import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/../../../../writeup/code/study/lowzero/deg2/fresh_predictions.json'))['gue'])
LAW = {}
for f in glob.glob(D + 'lawout/*.json'):
    r = json.load(open(f)); LAW[r['id']] = r
EXCL = ('N2357_0_m1_1_m3_1', 'N427_1_0_1_m8_7')
LS = np.arange(0.5, 10.01, 0.5); I = lambda l: int(round(l / 0.5)) - 1
rng = np.random.default_rng(1)


def arr(ids, band, key):
    return np.array([LAW[i][band][key] for i in ids])


def meas(ids, band='full'):
    return np.array([S2[i][band] for i in ids])


def slope(y, x, w):
    y = y - y.mean(0); x = x - x.mean(0)
    return float((w * (x * y).sum(0)).sum() / (w * (x * x).sum(0)).sum())


def gam(S, T, A, nb=2000):
    """slope of (S - A) on (T - A), between-curve, weights 1/var(S - T) per l; bootstrap over curves"""
    w = 1 / (S - T).var(0, ddof=1)
    g0 = slope(S - A, T - A, w); n = len(S); bs = []
    for _ in range(nb):
        ix = rng.integers(0, n, n); wb = 1 / (S[ix] - T[ix]).var(0, ddof=1)
        bs.append(slope(S[ix] - A[ix], T[ix] - A[ix], wb))
    return g0, float(np.std(bs, ddof=1))


def level(S, T, nb=2000):
    m = (S - T).mean(0); n = len(S)
    bs = np.array([(S[ix] - T[ix]).mean(0) for ix in (rng.integers(0, n, n) for _ in range(nb))])
    return m, bs.std(0, ddof=1)


def fmt(v, f='%+.4f'):
    return ' '.join(f % x for x in v)


def main():
    ids_all = sorted(i for i in LAW if i in S2)
    ids = [i for i in ids_all if i not in EXCL]
    print('curves with own law:', len(ids_all), ' kept (registered rule):', len(ids))
    # 1. own law vs frozen
    dl = np.array([np.array(LAW[i]['full']['full']) - np.array(P[i]['full']) for i in ids_all])
    print('\n1. own law (24 nodes) - frozen law (10 nodes): mean', fmt(dl.mean(0)[[1, 3, 7, 11, 15, 17, 19]], '%+.5f'))
    print('   max|.| per l', fmt(np.abs(dl).max(0)[[1, 3, 7, 11, 15, 17, 19]], '%.5f'), ' (l = 1,2,4,6,8,9,10)')
    print('   rms per l   ', fmt(np.sqrt((dl ** 2).mean(0))[[1, 3, 7, 11, 15, 17, 19]], '%.5f'))
    S = meas(ids); T = arr(ids, 'full', 'full'); Tf = np.array([P[i]['full'] for i in ids])
    for name, TT in (('own law', T), ('frozen law', Tf)):
        m, se = level(S, TT)
        print('\n2. level data - %s (%d curves):' % (name, len(ids))); print('   ' + fmt(m)); print('   se ' + fmt(se))
        Dd = S - TT; n, p = Dd.shape; mm = Dd.mean(0); C = np.cov(Dd.T); T2 = n * mm @ np.linalg.solve(C, mm)
        F = T2 * (n - p) / ((n - 1) * p); print('   Hotelling T2 %.1f, F-p %.3f' % (T2, stats.f.sf(F, p, n - p)))
    # 3. own regressions
    print('\n3. own between-curve amplitudes (own law, own templates):')
    A = {'gamma_A': T * 0, 'gamma_deg1': arr(ids, 'full', 'deg1'), 'gamma_nosym (Sym2 ratio -> 1)': arr(ids, 'full', 'nosym'),
         'gamma_k1p2 (doubled orbit)': T + arr(ids, 'full', 'dk1p2'), 'gamma_w1p2 (first traversal of 2)': T + arr(ids, 'full', 'dw1p2'),
         'gamma_k1 (all repetitions)': T + arr(ids, 'full', 'dk1'), 'gamma_k1p3': T + arr(ids, 'full', 'dk1p3'),
         'gamma_w1p3': T + arr(ids, 'full', 'dw1p3')}
    G = {}
    for k, alt in A.items():
        if k == 'gamma_A':
            G[k] = gam(S, T, np.zeros_like(T))
        else:
            G[k] = gam(S, T, alt)
        print('   %-36s %.4f +- %.4f' % (k, *G[k]))
    # the same readings with the 138 curves (duplicates removed)
    S138 = meas(ids_all); T138 = arr(ids_all, 'full', 'full')
    m, se = level(S138, T138)
    print('\n   138 curves (two duplicated pairs removed): level', fmt(m[[1, 3, 7, 11, 15, 19]]), ' se', fmt(se[[1, 3, 7, 11, 15, 19]], '%.4f'))
    for k in ('gamma_A', 'gamma_deg1', 'gamma_k1p2 (doubled orbit)', 'gamma_w1p2 (first traversal of 2)', 'gamma_nosym (Sym2 ratio -> 1)'):
        if k == 'gamma_A': alt = np.zeros_like(T138)
        elif k == 'gamma_deg1': alt = arr(ids_all, 'full', 'deg1')
        elif k.startswith('gamma_k1p2'): alt = T138 + arr(ids_all, 'full', 'dk1p2')
        elif k.startswith('gamma_w1p2'): alt = T138 + arr(ids_all, 'full', 'dw1p2')
        else: alt = arr(ids_all, 'full', 'nosym')
        print('   138: %-34s %.4f +- %.4f' % (k, *gam(S138, T138, alt)))
    for i in EXCL:
        if i in LAW:
            print('   %s: data - law at l = 1,2,4,6,8,10:' % i, fmt((np.array(S2[i]['full']) - np.array(LAW[i]['full']['full']))[[1, 3, 7, 11, 15, 19]]))
    return ids, ids_all, S, T


if __name__ == '__main__':
    main()
