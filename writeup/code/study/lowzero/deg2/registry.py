"""Registered numbers from a predictions JSON alone (no zeros): design, mean levels, regression-matched per-prime
contrasts of the law split into traversals (k = 1, 2, >= 3 of p = 2 and 3) and the rest, and the same under the
alternatives. usage: registry.py PRED.json"""
import sys, json
import numpy as np
P = json.load(open(sys.argv[1])); C = P['curves']; LS = P['l']; G = np.array(P['gue'])
show = [LS.index(l) for l in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0)]
def row(v): return ' '.join('%+.4f' % v[i] for i in show)
print('curves %d; band %s; l shown: %s' % (len(C), P['band'], ' '.join('%7.1f' % LS[i] for i in show)))
for g in ('g2', 'g3'):
    print('  %s groups:' % g, {k: sum(c[g] == k for c in C) for k in sorted(set(c[g] for c in C))})
print('  ranks/root numbers: w=+1 %d, w=-1 %d; N from %d to %d; mean L %.2f (%.2f .. %.2f); expected zeros %d (%.0f per curve)' % (
    sum(c['w'] == 1 for c in C), sum(c['w'] == -1 for c in C), min(c['N'] for c in C), max(c['N'] for c in C),
    np.mean([c['meanL'] for c in C]), min(c['meanL'] for c in C), max(c['meanL'] for c in C), sum(c['nz_expected'] for c in C), np.mean([c['nz_expected'] for c in C])))
M = {k: np.array([c[k] for c in C]) for k in ('full', 'k1', 'deg1', 'offzeta', 'k1p2', 'k1p3', 'w1p2', 'w1p3')}
print('mean level, l shown:')
print('  GUE      ' + row(G))
for k in ('full', 'deg1', 'k1', 'offzeta'):
    print('  %-8s ' % k + row(M[k].mean(0)))
print('  sd(full) ' + ' '.join('%.4f' % M['full'].std(0, ddof=1)[i] for i in show))
print('  max(full - GUE) ' + row((M['full'] - G).max(0)) + '   curves with full - GUE > 0.02 at l = 6: %d' % int(((M['full'] - G)[:, LS.index(6.0)] > 0.02).sum()))
for k in ('k1', 'deg1', 'offzeta', 'k1p2', 'w1p2'):
    d = M['full'] - M[k]
    print('  template full - %-7s between-curve sd: ' % k + ' '.join('%.4f' % d.std(0, ddof=1)[i] for i in show))
G2 = ['0', '2', 'b']; G3 = ['0', '2', '3', 'b']
X = np.array([[1.0] + [float(c['g2'] == g) for g in G2] + [float(c['g3'] == g) for g in G3] for c in C]); Xi = np.linalg.pinv(X.T @ X)
def coef(Y): return (Xi @ X.T @ Y)[1:]
pieces = {}
for p in (2, 3):
    for k in ('1', '2', '3'):
        pieces['d%d_%s' % (p, k)] = np.array([c['d%d_%s' % (p, k)] for c in C])
rest = M['full'] - sum(pieces.values())
print('regression-matched per-prime contrasts (baseline |a_2| = 1, |a_3| = 1), law split by traversal:')
names = ['|a2|=0', '|a2|=2', '2|N', '|a3|=0', '|a3|=2', '|a3|=3', '3|N']
for j, nm in enumerate(names):
    p = 2 if j < 3 else 3
    print('  %-7s law   ' % nm + row(coef(M['full'])[j]))
    for k in ('1', '2', '3'):
        print('  %-7s  k=%s  ' % ('', k if k != '3' else '3+') + row(coef(pieces['d%d_%s' % (p, k)])[j]))
    print('  %-7s  rest ' % '' + row(coef(rest)[j]) + '   (off-diagonal, other primes, L)')
    for a in ('k1', 'deg1'):
        print('  %-7s %-5s ' % ('', a) + row(coef(M[a])[j]))
