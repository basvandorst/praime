# All curves: Sigma2 = var(dS_X) [primes <= X, time domain, no zeros] + var(dR_X) [remainder from zeros] + 2 cov.
# (a) cross term; (b) remainder test var(dR_X) vs law - D_X (the conjectural part, with the prime-line scatter removed);
# (c) registered readings with the law corrected by the exact finite-band explicit formula (law - D_X + var(dS_X)).
import json, glob, numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)) + '/')
import analyse as A
from scipy import stats
E = {}
for f in glob.glob(A.D + 'e2all/*.json'):
    r = json.load(open(f)); E[r['id']] = r
ids = sorted(i for i in E if i not in A.EXCL); ids_all = sorted(E)
print('curves', len(ids), '(+ %d dedup)' % (len(ids_all) - len(ids)))
pick = [0, 1, 3, 7, 11, 15, 19]
S = A.meas(ids); T = A.arr(ids, 'full', 'full'); n = len(ids)
print('max |var(true dS) - registered S2|:', max(np.abs(np.array(E[i]['meas']) - np.array(A.S2[i]['full'])).max() for i in ids))
for X in ('1000', '10000', '100000'):
    sX = np.array([E[i]['X'][X]['s2X'] for i in ids]); sR = np.array([E[i]['X'][X]['s2R'] for i in ids]); DX = A.arr(ids, 'full', 'D' + X)
    cr = S - sX - sR
    print('\nX = %s' % X)
    print('  2cov(dS_X, dR_X): mean %s  max|.| %s' % (A.fmt(cr.mean(0)[pick], '%+.5f'), A.fmt(np.abs(cr).max(0)[pick], '%.4f')))
    print('  var(dS_X) - D_X : mean %s  sd %s' % (A.fmt((sX - DX).mean(0)[pick]), A.fmt((sX - DX).std(0)[pick], '%.4f')))
    print('  var(dR_X) (zeros remainder) mean %s ; law - D_X %s' % (A.fmt(sR.mean(0)[pick], '%.4f'), A.fmt((T - DX).mean(0)[pick], '%.4f')))
    q = sR - (T - DX); m, se = A.level(sR, T - DX)
    C = np.cov(q.T); T2 = n * m @ np.linalg.solve(C, m); F = T2 * (n - 20) / ((n - 1) * 20)
    print('  remainder test var(dR_X) - (law - D_X): m %s' % A.fmt(m[pick]))
    print('                                       se %s  max|z| %.1f  Hotelling p %.3g' % (A.fmt(se[pick], '%.4f'), np.abs(m / se).max(), stats.f.sf(F, 20, n - 20)))
    print('  per-curve sd: data - law %s -> remainder test %s' % (A.fmt((S - T).std(0)[pick], '%.4f'), A.fmt(q.std(0)[pick], '%.4f')))
    for Y in ('1000', '10000', 'H'):
        EFr = A.arr(ids, 'full', 'D' + Y) + A.arr(ids, 'full', 'G' + Y) - DX
        mm, ss = A.level(sR, EFr, 500)
        print('  remainder vs EF_%-5s remainder (lines to %s + GUE beyond): m %s  max|z| %.1f' % (Y, Y, A.fmt(mm[pick]), np.abs(mm / ss).max()))
    # registered readings with the corrected law
    Tc = T - DX + sX
    print('  corrected law (law - D_X + var(dS_X)):')
    for k, alt in (('gamma_A', np.zeros_like(T)), ('gamma_deg1', A.arr(ids, 'full', 'deg1') - DX + sX),
                   ('gamma_k1p2', Tc + A.arr(ids, 'full', 'dk1p2')), ('gamma_w1p2', Tc + A.arr(ids, 'full', 'dw1p2')),
                   ('gamma_k1', Tc + A.arr(ids, 'full', 'dk1'))):
        print('     %-11s %.4f +- %.4f' % (k, *A.gam(S, Tc, alt, 1000)))
